/**
 * @file Acoustic-Forward.cpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 一阶速度应力方程交错网格形式
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
#include"acoustic_forward_2d.h"
#include"fwi_load_point_source.hpp"
int main(int argc, char**argv)
{
	fwi_cmd_t cmd;
	ParseCommandLine(cmd,argc,argv);
	DebugCommandLine(cmd);
	int NX = 0;
	int NZ = 0;
	int NT = 0;
	int CDO_ORDER = 0;
	int ABC_LAYERS = 0;
	float DX = 0.0f;
	float DZ = 0.0f;
	float DT = 0.0f;
	//*********************所有的维度参数***************
	GetCmdValueEx(cmd,"NX",0,NX);
	GetCmdValueEx(cmd,"NZ",0,NZ);
	GetCmdValueEx(cmd, "NT", 0, NT);
	GetCmdValueEx(cmd, "DX", 0, DX);
	GetCmdValueEx(cmd, "DZ", 0, DZ);
	GetCmdValueEx(cmd, "DT", 0, DT);
	GetCmdValueEx(cmd, "CDO_ORDER", 0, CDO_ORDER);
	GetCmdValueEx(cmd,"ABC_LAYERS",0,ABC_LAYERS);
	int XDIM = NX + 2 * (CDO_ORDER+ABC_LAYERS);
	int ZDIM = NZ + 2 * (CDO_ORDER+ABC_LAYERS);
	int XLEN = NX + 2 * ABC_LAYERS;
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
	dense_model = FWIAllocateMemory(XDIM*ZDIM, 1000.0f);
	if (!dense_name.empty())
	{
		raw_dense_model = FWIAllocateMemory(NX*NZ,0.0f);
		Read2DBin(dense_name.c_str(),raw_dense_model,NX,NZ);
		Expand2D(dense_model,raw_dense_model,XDIM,ZDIM,NX,NZ,OFFSET);
	}
	vp_model = FWIAllocateMemory(XDIM*ZDIM, 3000.0f);
	raw_vp_model = FWIAllocateMemory(NX*NZ, 0.0f);
	Read2DBin(vp_name.c_str(), raw_vp_model, NX, NZ);
	Expand2D(vp_model, raw_vp_model, XDIM, ZDIM, NX, NZ, OFFSET);
	//检查模型
	//Write2DBin("ac_2d_dense_model", dense_model, XDIM, ZDIM);
	//Write2DBin("ac_2d_vp_model",vp_model,XDIM,ZDIM);

	float vp_info[4] = { 0.0f,0.0f,0.0f,0.0f };
	GetCpmlAverageVpModelRect(vp_info,vp_model,XDIM,ZDIM,CDO_ORDER,ABC_LAYERS);

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
	int LOAD_ZPOS = 0;
	GetCmdValueEx(cmd,"LOAD_XPOS",1,LOAD_XPOS,NX/2);
	GetCmdValueEx(cmd, "LOAD_ZPOS", 1, LOAD_ZPOS, NZ / 2);
	int LOAD_GLOBAL_XPOS = LOAD_XPOS + OFFSET;
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
	float* z1_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* z2_cpml_alpha = FWIAllocateMemory(ABC_LAYERS);
	float* x1_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* x2_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* z1_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	float* z2_cpml_beta = FWIAllocateMemory(ABC_LAYERS);
	//吸收边界系数
	const float DAMPING_RATE = 1E-6f;
	Cpmlfun(x1_cpml_alpha, x1_cpml_beta, ABC_LAYERS, DX, DT, MAIN_FREQ, vp_info[0], DAMPING_RATE);
	Cpmlfun(x2_cpml_alpha, x2_cpml_beta, ABC_LAYERS, DX, DT, MAIN_FREQ, vp_info[1], DAMPING_RATE);
	Cpmlfun(z1_cpml_alpha, z1_cpml_beta, ABC_LAYERS, DZ, DT, MAIN_FREQ, vp_info[2], DAMPING_RATE);
	Cpmlfun(z2_cpml_alpha, z2_cpml_beta, ABC_LAYERS, DZ, DT, MAIN_FREQ, vp_info[3], DAMPING_RATE);
	//吸收边界辅助变量
	float* vx_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN,0.0f);
	float* vx_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN, 0.0f);
	float* stress_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN, 0.0f);
	float* stress_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN, 0.0f);

	float* vz_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN, 0.0f);
	float* vz_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN, 0.0f);
	float* stress_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN, 0.0f);
	float* stress_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN, 0.0f);
	//*********************检波器测线*********************
	//和地表平行的测线
	int XLINE_ZPOS = 0;
	int XLINE_XBEGIN = 0;
	int XLINE_XINTERVAL = 1;
	int XLINE_XEND = NX;
	std::string XLINE_NAME;
	GetCmdValueEx(cmd,"XLINE_ZPOS",1,XLINE_ZPOS,0);
	GetCmdValueEx(cmd, "XLINE_XBEGIN", 1, XLINE_XBEGIN, 0);
	GetCmdValueEx(cmd, "XLINE_XINTERVAL", 1, XLINE_XINTERVAL, 1);
	GetCmdValueEx(cmd, "XLINE_XEND", 1, XLINE_XEND, NX);
	GetCmdValueEx(cmd,"XLINE_NAME",1,XLINE_NAME,"xline");
	//垂直的测线
	int ZLINE_XPOS = LOAD_XPOS;
	int ZLINE_ZBEGIN = 0;
	int ZLINE_ZINTERVAL = 1;
	int ZLINE_ZEND = NZ;
	std::string ZLINE_NAME;
	GetCmdValueEx(cmd, "ZLINE_XPOS", 1, ZLINE_XPOS, 0);
	GetCmdValueEx(cmd, "ZLINE_ZBEGIN", 1, ZLINE_ZBEGIN, 0);
	GetCmdValueEx(cmd, "ZLINE_ZINTERVAL", 1, ZLINE_ZINTERVAL, 1);
	GetCmdValueEx(cmd, "ZLINE_ZEND", 1, ZLINE_ZEND, NZ);
	GetCmdValueEx(cmd, "ZLINE_NAME", 1, ZLINE_NAME, "zline");

	//检波器缓冲区
	float* stress_xline_geo = FWIAllocateMemory(NT*NX,0.0f);
	float* vx_xline_geo = FWIAllocateMemory(NT*NX, 0.0f);
	float* vz_xline_geo = FWIAllocateMemory(NT*NX, 0.0f);

	float* stress_zline_geo = FWIAllocateMemory(NT*NZ, 0.0f);
	float* vx_zline_geo = FWIAllocateMemory(NT*NZ, 0.0f);
	float* vz_zline_geo = FWIAllocateMemory(NT*NZ, 0.0f);
	//申请波场
	float* stress = FWIAllocateMemory(XDIM*ZDIM,0.0f);
	float* vx = FWIAllocateMemory(XDIM*ZDIM, 0.0f);
	float* vz = FWIAllocateMemory(XDIM*ZDIM, 0.0f);
	TIME_BEGIN("Ac2D")
	for (int t = 0; t < NT; ++t)
	{
		LoadPointSourceOp::LoadAcStressSource(stress,XDIM,ZDIM,LOAD_GLOBAL_XPOS,LOAD_GLOBAL_ZPOS,source_buffer[t]);
		AcSgStressToVelocityCPML2D(vx,vz,stress,dense_model,XDIM,ZDIM,XLEN,ZLEN,CDO_ORDER,sg_cdo_coeff,
		ABC_LAYERS, x1_cpml_alpha, x1_cpml_beta, x2_cpml_alpha, x2_cpml_beta,
			z1_cpml_alpha, z1_cpml_beta, z2_cpml_alpha, z2_cpml_beta,stress_dx_phi1,stress_dx_phi2,
			stress_dz_phi1,stress_dz_phi2,DX,DZ,DT);
		AcSgVelocityToStressCPML2D(stress,vx,vz,dense_model,vp_model,XDIM,ZDIM,XLEN,ZLEN,CDO_ORDER,
			sg_cdo_coeff,ABC_LAYERS,x1_cpml_alpha,x1_cpml_beta,x2_cpml_alpha,x2_cpml_beta,
			z1_cpml_alpha,z1_cpml_beta,z2_cpml_alpha,z2_cpml_beta,vx_dx_phi1,vx_dx_phi2,
			vz_dz_phi1,vz_dz_phi2,DX,DZ,DT);
		//添加检波器记录的
		for (int i = 0; i < NX; ++i)
		{
			const size_t record_to_pos = t * NX + i;
			const size_t record_from_pos = (i + OFFSET)*ZDIM + XLINE_ZPOS+OFFSET;
			stress_xline_geo[record_to_pos] = stress[record_from_pos];
			vx_xline_geo[record_to_pos] = vx[record_from_pos];
			vz_xline_geo[record_to_pos] = vz[record_from_pos];
		}
		for (int j = 0; j < NZ; ++j)
		{
			const size_t record_to_pos = t * NZ + j;
			const size_t record_from_pos = (ZLINE_XPOS + OFFSET)*ZDIM + j + OFFSET;
			stress_zline_geo[record_to_pos] = stress[record_from_pos];
			vx_zline_geo[record_to_pos] = vx[record_from_pos];
			vz_zline_geo[record_to_pos] = vz[record_from_pos];
		}
	}
	TIME_END("Ac2D")
	Write2DBinEx(XLINE_NAME.c_str(),"-stress-",stress_xline_geo,NT,NX);
	Write2DBinEx(XLINE_NAME.c_str(), "-vx-", vx_xline_geo, NT, NX);
	Write2DBinEx(XLINE_NAME.c_str(), "-vz-", vz_xline_geo, NT, NX);

	Write2DBinEx(ZLINE_NAME.c_str(), "-stress-", stress_zline_geo, NT, NZ);
	Write2DBinEx(ZLINE_NAME.c_str(), "-vx-", vx_zline_geo, NT, NZ);
	Write2DBinEx(ZLINE_NAME.c_str(), "-vz-", vz_zline_geo, NT, NZ);

	//清理内存
	DebugFreeMemory(&sg_cdo_coeff);
	DebugFreeMemory(&x1_cpml_alpha);
	DebugFreeMemory(&x1_cpml_beta);
	DebugFreeMemory(&x2_cpml_alpha);
	DebugFreeMemory(&x2_cpml_beta);
	DebugFreeMemory(&z1_cpml_alpha);
	DebugFreeMemory(&z1_cpml_beta);
	DebugFreeMemory(&z2_cpml_alpha);
	DebugFreeMemory(&z2_cpml_beta);
	DebugFreeMemory(&vx_dx_phi1);
	DebugFreeMemory(&vx_dx_phi2);
	DebugFreeMemory(&vz_dz_phi1);
	DebugFreeMemory(&vz_dz_phi2);
	DebugFreeMemory(&stress_dx_phi1);
	DebugFreeMemory(&stress_dx_phi2);
	DebugFreeMemory(&stress_dz_phi1);
	DebugFreeMemory(&stress_dz_phi2);
	DebugFreeMemory(&raw_dense_model);
	DebugFreeMemory(&raw_vp_model);
	DebugFreeMemory(&dense_model);
	DebugFreeMemory(&vp_model);
	DebugFreeMemory(&source_buffer);
	DebugFreeMemory(&stress);
	DebugFreeMemory(&vx);
	DebugFreeMemory(&vz);
	DebugFreeMemory(&stress_xline_geo);
	DebugFreeMemory(&vx_xline_geo);
	DebugFreeMemory(&vz_xline_geo);
	DebugFreeMemory(&stress_zline_geo);
	DebugFreeMemory(&vx_zline_geo);
	DebugFreeMemory(&vz_zline_geo);
	return 0;
}