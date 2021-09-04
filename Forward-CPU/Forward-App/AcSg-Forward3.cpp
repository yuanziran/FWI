/**
 * @file AcSg-Forward3.cpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 三维声波正演模拟
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */


#include"fwi.h"
#include"fwi_cmd.hpp"
#include"fwi_memory.hpp"
#include"fwi_io.hpp"
#include"fwi_tool.hpp"
#include"fwi_load_source_func.hpp"
#include"fwi_coeff.hpp"
#include"acoustic_forward_3d.h"
#include"fwi_load_point_source.hpp"
int main(int argc, char**argv)
{
	fwi_cmd_t cmd;
	ParseCommandLine(cmd,argc,argv);
	DebugCommandLine(cmd);
	int NX = 0;
	int NY = 0;
	int NZ = 0;
	int NT = 0;
	int CDO_ORDER = 0;
	int ABC_LAYERS = 0;
	float DX = 0.0f;
	float DY = 0.0f;
	float DZ = 0.0f;
	float DT = 0.0f;
	//*********************所有的维度参数***************
	GetCmdValueEx(cmd,"NX",0,NX);
	GetCmdValueEx(cmd, "NY", 0, NY);
	GetCmdValueEx(cmd,"NZ",0,NZ);
	GetCmdValueEx(cmd, "NT", 0, NT);
	GetCmdValueEx(cmd, "DX", 0, DX);
	GetCmdValueEx(cmd, "DY", 0, DY);
	GetCmdValueEx(cmd, "DZ", 0, DZ);
	GetCmdValueEx(cmd, "DT", 0, DT);
	GetCmdValueEx(cmd, "CDO_ORDER", 0, CDO_ORDER);
	GetCmdValueEx(cmd,"ABC_LAYERS",0,ABC_LAYERS);
	int XDIM = NX + 2 * (CDO_ORDER+ABC_LAYERS);
	int YDIM = NY + 2 * (CDO_ORDER + ABC_LAYERS);
	int ZDIM = NZ + 2 * (CDO_ORDER+ABC_LAYERS);
	int XLEN = NX + 2 * ABC_LAYERS;
	int YLEN = NY + 2 * ABC_LAYERS;
	int ZLEN = NZ + 2 * ABC_LAYERS;
	int OFFSET = CDO_ORDER + ABC_LAYERS;
	//*******************速度模型***********************
	float* dense_model = nullptr;
	float* vp_model = nullptr;
	float* raw_dense_model = nullptr;
	float* raw_vp_model = nullptr;
	std::string dense_name;
	std::string vp_name;
	GetCmdValueEx(cmd,"DENSE_NAME",1,dense_name);
	GetCmdValueEx(cmd, "VP_NAME", 0, vp_name);
	dense_model = FWIAllocateMemory(XDIM*YDIM*ZDIM, 1000.0f);
	if (!dense_name.empty())
	{
		raw_dense_model = FWIAllocateMemory(NX*NY*NZ,0.0f);
		Read3DBin(dense_name.c_str(),raw_dense_model,NX,NY,NZ);
		Expand3D(dense_model,raw_dense_model,XDIM,YDIM,ZDIM,NX,NY,NZ,OFFSET);
	}
	vp_model = FWIAllocateMemory(XDIM*YDIM*ZDIM, 3000.0f);
	raw_vp_model = FWIAllocateMemory(NX*NY*NZ, 0.0f);
	Read3DBin(vp_name.c_str(), raw_vp_model, NX,NY,NZ);
	Expand3D(vp_model, raw_vp_model, XDIM,YDIM,ZDIM, NX,NY,NZ, OFFSET);
	//检查模型
	Write3DBin("ac_3d_dense_model", dense_model, XDIM,YDIM,ZDIM);
	Write3DBin("ac_3d_vp_model",vp_model,XDIM,YDIM,ZDIM);

	float vp_info[6] = { 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
	GetCpmlAverageVpModelBox(vp_info,vp_model,XDIM,YDIM,ZDIM,CDO_ORDER,ABC_LAYERS);

	//*************************有限差分系数*********************
	float * sg_cdo_coeff = FWIAllocateMemory(CDO_ORDER);
	SgCdoOp1(sg_cdo_coeff,CDO_ORDER);
	//*********************获取震源参数***************************************
	float MAIN_FREQ = 15.0f;
	float TIME_DELAY = 0.0f;
	float MAX_AMP = 1.0f;
	GetCmdValueEx(cmd,"MAIN_FREQ",1,MAIN_FREQ,15.0f);
	GetCmdValueEx(cmd, "TIME_DELAY", 1, TIME_DELAY, 0.0f);
	GetCmdValueEx(cmd,"MAX_AMP",1,MAX_AMP,1.0f);
	int LOAD_XPOS = 0;
	int LOAD_YPOS = 0;
	int LOAD_ZPOS = 0;
	GetCmdValueEx(cmd,"LOAD_XPOS",1,LOAD_XPOS,NX/2);
	GetCmdValueEx(cmd, "LOAD_XPOS", 1, LOAD_YPOS, NY / 2);
	GetCmdValueEx(cmd, "LOAD_ZPOS", 1, LOAD_ZPOS, NZ / 2);
	int LOAD_GLOBAL_XPOS = LOAD_XPOS + OFFSET;
	int LOAD_GLOBAL_YPOS = LOAD_YPOS + OFFSET;
	int LOAD_GLOBAL_ZPOS = LOAD_ZPOS + OFFSET;
	//准备震源缓冲区
	float* source_buffer = FWIAllocateMemory(NT,0);
	int LOAD_TEND =0;
	GetCmdValueEx(cmd,"LOAD_TEND",1,LOAD_TEND,NT/2);
	//使用雷克子波填充震源
	for (int i = 0; i < LOAD_TEND; ++i)
	{
		source_buffer[i]=RickerSource(MAIN_FREQ,DT*i,TIME_DELAY,MAX_AMP);
	}
	//****************************CPML吸收边界系数***********************
	float* x1_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* x2_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* y1_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* y2_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* z1_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* z2_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* x1_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* x2_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* y1_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* y2_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* z1_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* z2_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	//吸收边界系数
	const float DAMPING_RATE = 1E-6f;
	Cpmlfun(x1_cpml_alpha, x1_cpml_beta, ABC_LAYERS, DX, DT, MAIN_FREQ, vp_info[0], DAMPING_RATE);
	Cpmlfun(x2_cpml_alpha, x2_cpml_beta, ABC_LAYERS, DX, DT, MAIN_FREQ, vp_info[1], DAMPING_RATE);
	Cpmlfun(y1_cpml_alpha, y1_cpml_beta, ABC_LAYERS, DY, DT, MAIN_FREQ, vp_info[2], DAMPING_RATE);
	Cpmlfun(y2_cpml_alpha, y2_cpml_beta, ABC_LAYERS, DY, DT, MAIN_FREQ, vp_info[3], DAMPING_RATE);
	Cpmlfun(z1_cpml_alpha, z1_cpml_beta, ABC_LAYERS, DZ, DT, MAIN_FREQ, vp_info[4], DAMPING_RATE);
	Cpmlfun(z2_cpml_alpha, z2_cpml_beta, ABC_LAYERS, DZ, DT, MAIN_FREQ, vp_info[5], DAMPING_RATE);
	//吸收边界辅助变量
	float* vx_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN,0.0f);
	float* vx_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
	float* stress_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
	float* stress_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);

	float* vy_dy_phi1= FWIAllocateMemory(ABC_LAYERS*ZLEN*XLEN, 0.0f);
	float* vy_dy_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*XLEN, 0.0f);
	float* stress_dy_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*XLEN, 0.0f);
	float* stress_dy_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*XLEN, 0.0f);

	float* vz_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* vz_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* stress_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* stress_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	//*********************检波器测线*********************
	//和地表平行的二维测面
	int ZPLANE_ZPOS = 0;
	int ZPLANE_XBEGIN = 0;
	int ZPLANE_XINTERVAL = 1;
	int ZPLANE_XEND = NX;
	int ZPLANE_YBEGIN = 0;
	int ZPLANE_YINTERVAL = 1;
	int ZPLANE_YEND = NY;
	std::string ZPLANE_NAME;
	GetCmdValueEx(cmd,"ZPLANE_ZPOS",1,ZPLANE_ZPOS,0);
	GetCmdValueEx(cmd, "ZPLANE_XBEGIN", 1, ZPLANE_XBEGIN, 0);
	GetCmdValueEx(cmd, "ZPLANE_XINTERVAL", 1, ZPLANE_XINTERVAL, 1);
	GetCmdValueEx(cmd, "ZPLANE_XEND", 1, ZPLANE_XEND, NX);
	GetCmdValueEx(cmd, "ZPLANE_YBEGIN", 1, ZPLANE_YBEGIN, 0);
	GetCmdValueEx(cmd, "ZPLANE_YINTERVAL", 1, ZPLANE_YINTERVAL, 1);
	GetCmdValueEx(cmd, "ZPLANE_YEND", 1, ZPLANE_YEND, NY);
	GetCmdValueEx(cmd,"ZPLANE_NAME",1,ZPLANE_NAME,"zplane");

	//检波器缓冲区
	float* stress_zplane_geo = FWIAllocateMemory(NT*NX*NY,0.0f);
	float* vx_zplane_geo = FWIAllocateMemory(NT*NX*NY, 0.0f);
	float* vy_zplane_geo = FWIAllocateMemory(NT*NX*NY, 0.0f);
	float* vz_zplane_geo = FWIAllocateMemory(NT*NX*NY, 0.0f);

	//申请波场
	float* stress = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
	float* vx = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	float* vy = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	float* vz = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	TIME_BEGIN("Ac3D")
	for (int t = 0; t < NT; ++t)
	{
		std::cout << t << std::endl;
		LoadPointSourceOp::LoadAcStressSource3D(stress,XDIM,YDIM,ZDIM,LOAD_GLOBAL_XPOS, LOAD_GLOBAL_YPOS,LOAD_GLOBAL_ZPOS,source_buffer[t]);
		AcSgStressToVelocityCPML3D(vx,vy,vz,stress,dense_model,XDIM,YDIM,ZDIM,XLEN,YLEN,ZLEN,CDO_ORDER,sg_cdo_coeff,
		ABC_LAYERS, 
			x1_cpml_alpha, x1_cpml_beta, x2_cpml_alpha, x2_cpml_beta,
			y1_cpml_alpha, y1_cpml_beta, y2_cpml_alpha, y2_cpml_beta,
			z1_cpml_alpha, z1_cpml_beta, z2_cpml_alpha, z2_cpml_beta,
			stress_dx_phi1,stress_dx_phi2,
			stress_dy_phi1, stress_dy_phi2,
			stress_dz_phi1,stress_dz_phi2,
			DX,DY,DZ,DT);
		AcSgVelocityToStressCPML3D(stress,vx,vy,vz,dense_model,vp_model,XDIM,YDIM,ZDIM,XLEN,YLEN,ZLEN,CDO_ORDER,
			sg_cdo_coeff,ABC_LAYERS,x1_cpml_alpha,x1_cpml_beta,x2_cpml_alpha,x2_cpml_beta,
			y1_cpml_alpha, y1_cpml_beta, y2_cpml_alpha, y2_cpml_beta,
			z1_cpml_alpha,z1_cpml_beta,z2_cpml_alpha,z2_cpml_beta,
			vx_dx_phi1,vx_dx_phi2,
			vy_dy_phi1, vy_dy_phi2,
			vz_dz_phi1,vz_dz_phi2,
			DX,DY,DZ,DT);
		//添加检波器记录的
		for (int i = 0; i < NX; ++i)
		for(int j=0;j<NY;++j)
		{
			const size_t record_to_pos = (t * NX + i)*NY+j;
			const size_t record_from_pos = ((i + OFFSET)*YDIM + j+OFFSET)*ZDIM+ZPLANE_ZPOS+OFFSET;
			stress_zplane_geo[record_to_pos] = stress[record_from_pos];
			vx_zplane_geo[record_to_pos] = vx[record_from_pos];
			vy_zplane_geo[record_to_pos] = vy[record_from_pos];
			vz_zplane_geo[record_to_pos] = vz[record_from_pos];
		}
	}
	TIME_END("Ac3D")
	Write3DBinEx(ZPLANE_NAME.c_str(),"-stress-",stress_zplane_geo,NT,NX,NY);
	Write3DBinEx(ZPLANE_NAME.c_str(), "-vx-", vx_zplane_geo, NT, NX, NY);
	Write3DBinEx(ZPLANE_NAME.c_str(), "-vy-", vy_zplane_geo, NT, NX, NY);
	Write3DBinEx(ZPLANE_NAME.c_str(), "-vz-", vz_zplane_geo, NT, NX, NY);

	//清理内存
	DebugFreeMemory(&sg_cdo_coeff);
	DebugFreeMemory(&x1_cpml_alpha);
	DebugFreeMemory(&x1_cpml_beta);
	DebugFreeMemory(&x2_cpml_alpha);
	DebugFreeMemory(&x2_cpml_beta);
	DebugFreeMemory(&y1_cpml_alpha);
	DebugFreeMemory(&y1_cpml_beta);
	DebugFreeMemory(&y2_cpml_alpha);
	DebugFreeMemory(&y2_cpml_beta);
	DebugFreeMemory(&z1_cpml_alpha);
	DebugFreeMemory(&z1_cpml_beta);
	DebugFreeMemory(&z2_cpml_alpha);
	DebugFreeMemory(&z2_cpml_beta);
	DebugFreeMemory(&vx_dx_phi1);
	DebugFreeMemory(&vx_dx_phi2);
	DebugFreeMemory(&vy_dy_phi1);
	DebugFreeMemory(&vy_dy_phi2);
	DebugFreeMemory(&vz_dz_phi1);
	DebugFreeMemory(&vz_dz_phi2);
	DebugFreeMemory(&stress_dx_phi1);
	DebugFreeMemory(&stress_dx_phi2);
	DebugFreeMemory(&stress_dy_phi1);
	DebugFreeMemory(&stress_dy_phi2);
	DebugFreeMemory(&stress_dz_phi1);
	DebugFreeMemory(&stress_dz_phi2);

	DebugFreeMemory(&raw_dense_model);
	DebugFreeMemory(&raw_vp_model);
	DebugFreeMemory(&dense_model);
	DebugFreeMemory(&vp_model);
	DebugFreeMemory(&source_buffer);
	DebugFreeMemory(&stress);
	DebugFreeMemory(&vx);
	DebugFreeMemory(&vy);
	DebugFreeMemory(&vz);
	DebugFreeMemory(&stress_zplane_geo);
	DebugFreeMemory(&vx_zplane_geo);
	DebugFreeMemory(&vy_zplane_geo);
	DebugFreeMemory(&vz_zplane_geo);
	return 0;
}

