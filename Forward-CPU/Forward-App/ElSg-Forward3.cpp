/**
 * @file ElSg3D.cpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 弹性波一阶速度应力方程正演模拟(三维)
 * @version 0.1
 * @date 2021-06-29
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
#include"elastic_forward_3d.h"
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
    float DY=0.0f;
	float DZ = 0.0f;
	float DT = 0.0f;
	//*********************所有的维度参数***************
	GetCmdValueEx(cmd,"NX",0,NX);
    GetCmdValueEx(cmd,"NY",0,NY);
	GetCmdValueEx(cmd,"NZ",0,NZ);
	GetCmdValueEx(cmd, "NT", 0, NT);
	GetCmdValueEx(cmd, "DX", 0, DX);
    GetCmdValueEx(cmd, "DY", 0, DY);
	GetCmdValueEx(cmd, "DZ", 0, DZ);
	GetCmdValueEx(cmd, "DT", 0, DT);
	GetCmdValueEx(cmd, "CDO_ORDER", 0, CDO_ORDER);
	GetCmdValueEx(cmd,"ABC_LAYERS",0,ABC_LAYERS);
	int XDIM = NX + 2 * (CDO_ORDER+ABC_LAYERS);
    int YDIM=NY+2*(CDO_ORDER+ABC_LAYERS);
	int ZDIM = NZ + 2 * (CDO_ORDER+ABC_LAYERS);
	int XLEN = NX + 2 * ABC_LAYERS;
    int YLEN=NY+2*ABC_LAYERS;
	int ZLEN = NZ + 2 * ABC_LAYERS;
	int OFFSET = CDO_ORDER + ABC_LAYERS;
	//*******************速度模型***********************
	float* dense_model = nullptr;
	float* vp_model = nullptr;
	float* temp_model = nullptr;
	std::string dense_name;
	std::string vp_name;
    std::string vs_name;
	GetCmdValueEx(cmd,"DENSE_NAME",1,dense_name);
	GetCmdValueEx(cmd, "VP_NAME", 0, vp_name);
    GetCmdValueEx(cmd, "VS_NAME", 1, vs_name);
	dense_model = FWIAllocateMemory(XDIM*YDIM*ZDIM, 1000.0f);
    temp_model=FWIAllocateMemory(NX*NY*NZ,0.0f);
	if (!dense_name.empty())
	{
		Read3DBin(dense_name.c_str(),temp_model,NX,NY,NZ);
		Expand3D(dense_model,temp_model,XDIM,YDIM,ZDIM,NX,NY,NZ,OFFSET);
	}
	vp_model = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	Read3DBin(vp_name.c_str(), temp_model, NX,NY,NZ);
	Expand3D(vp_model, temp_model, XDIM,YDIM,ZDIM, NX,NY,NZ, OFFSET);
    float* vs_model = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	if(!vs_name.empty())
    {
        Read3DBin(vs_name.c_str(),temp_model,NX,NY,NZ);
        Expand3D(vs_model,vs_model,XDIM,YDIM,ZDIM,NX,NY,NZ,OFFSET);
    }
    else
    {
        float scalar=1.0f/1.73f;
        for(int i=0;i<(XDIM*YDIM*ZDIM);++i)
        {
            vs_model[i]=scalar*vp_model[i];
        }
    }
    //检查模型
	//Write2DBin("el_2d_dense_model", dense_model, XDIM, ZDIM);
	//Write2DBin("el_2d_vp_model",vp_model,XDIM,ZDIM);
    //Write2DBin("el_2d_vs_model",vs_model,XDIM,ZDIM);
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
    int LOAD_YPOS=0;
	int LOAD_ZPOS = 0;
	GetCmdValueEx(cmd,"LOAD_XPOS",1,LOAD_XPOS,NX/2);
    GetCmdValueEx(cmd, "LOAD_ZPOS", 1, LOAD_YPOS, NY / 2);
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
    float* vy_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN,0.0f);
	float* vy_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
    float* vz_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN,0.0f);
	float* vz_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
	float* sigma_x_dx_phi1 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
	float* sigma_x_dx_phi2 = FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
    float* tau_y_dx_phi1=FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
    float* tau_y_dx_phi2=FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
    float* tau_z_dx_phi1=FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);
    float* tau_z_dx_phi2=FWIAllocateMemory(ABC_LAYERS*ZLEN*YLEN, 0.0f);

    float* vy_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* vy_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* vx_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* vx_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* vz_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* vz_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* sigma_y_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* sigma_y_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* tau_x_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* tau_x_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* tau_z_dy_phi1=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);
    float* tau_z_dy_phi2=FWIAllocateMemory(ABC_LAYERS*XLEN*ZLEN,0.0f);

    float* vx_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* vx_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
    float* vy_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* vy_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* vz_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* vz_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* sigma_z_dz_phi1 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* sigma_z_dz_phi2 = FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
    float* tau_x_dz_phi1= FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
    float* tau_x_dz_phi2= FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
	float* tau_y_dz_phi1= FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
    float* tau_y_dz_phi2= FWIAllocateMemory(ABC_LAYERS*XLEN*YLEN, 0.0f);
    //*********************检波器测线**********************************
	//和地表平行的测线
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
    int ZPLANE_RECORD_XLEN=GetLen(ZPLANE_XBEGIN,ZPLANE_XINTERVAL,ZPLANE_XEND);
    int ZPLANE_RECORD_YLEN=GetLen(ZPLANE_YBEGIN,ZPLANE_YINTERVAL,ZPLANE_YEND);
	//检波器缓冲区
	//float* sigma_xx_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN,0.0f);
    //float* sigma_zz_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN,0.0f);
    //float* tau_y_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN,0.0f);
	float* vx_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN, 0.0f);
    float* vy_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN, 0.0f);
	float* vz_xline_geo = FWIAllocateMemory(NT*ZPLANE_RECORD_XLEN*ZPLANE_RECORD_YLEN, 0.0f);

	
	//申请波场
	float* sigma_xx = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
    float* sigma_yy = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
    float* sigma_zz = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
    float* tau_x = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
    float* tau_y = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
    float* tau_z = FWIAllocateMemory(XDIM*YDIM*ZDIM,0.0f);
	float* vx = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
    float* vy = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	float* vz = FWIAllocateMemory(XDIM*YDIM*ZDIM, 0.0f);
	TIME_BEGIN("El2D")
	for (int t = 0; t < NT; ++t)
	{
		//LoadPointSourceOp::LoadAcStressSource(stress,XDIM,ZDIM,LOAD_GLOBAL_XPOS,LOAD_GLOBAL_ZPOS,source_buffer[t]);
		LoadPointSourceOp::LoadPSource(vx,vz,XDIM,ZDIM,LOAD_GLOBAL_XPOS,LOAD_GLOBAL_ZPOS,source_buffer[t]);
        ElSgStressToVelocityCPML3D(vx,vy,vz,sigma_xx,sigma_yy,sigma_zz,tau_x,tau_y,tau_z,
        dense_model,XDIM,YDIM,ZDIM,XLEN,YLEN,ZLEN,CDO_ORDER,sg_cdo_coeff,
		ABC_LAYERS, x1_cpml_alpha, x1_cpml_beta, x2_cpml_alpha, x2_cpml_beta,
        y1_cpml_alpha, y1_cpml_beta, y2_cpml_alpha, y2_cpml_beta,
		z1_cpml_alpha, z1_cpml_beta, z2_cpml_alpha, z2_cpml_beta,
        sigma_x_dx_phi1,sigma_x_dx_phi2,tau_y_dx_phi1,tau_y_dx_phi2,tau_z_dx_phi1,tau_z_dx_phi2,
		sigma_y_dy_phi1,sigma_y_dy_phi2,tau_x_dy_phi1,tau_x_dy_phi2,tau_z_dy_phi1,tau_z_dy_phi2,
        sigma_z_dz_phi1,sigma_z_dz_phi2,tau_x_dz_phi1,tau_x_dz_phi2,tau_y_dx_phi1,tau_y_dx_phi2,
        DX,DY,DZ,DT);
		ElSgVelocityToStressCPML3D(sigma_xx,sigma_yy,sigma_zz,tau_x,tau_y,tau_z,vx,vy,vz,
        dense_model,vp_model,vs_model,XDIM,YDIM,ZDIM,XLEN,YLEN,ZLEN,CDO_ORDER,
			sg_cdo_coeff,ABC_LAYERS,
            x1_cpml_alpha,x1_cpml_beta,x2_cpml_alpha,x2_cpml_beta,
            y1_cpml_alpha,y1_cpml_beta,y2_cpml_alpha,y2_cpml_beta,
			z1_cpml_alpha,z1_cpml_beta,z2_cpml_alpha,z2_cpml_beta,
            vx_dx_phi1,vx_dx_phi2,vy_dx_phi1,vy_dx_phi2,vz_dx_phi1,vz_dx_phi2,
			vy_dy_phi1,vy_dy_phi2,vx_dy_phi1,vx_dy_phi2,vz_dy_phi1,vz_dy_phi2,
            vz_dz_phi1,vz_dz_phi2,vx_dz_phi1,vx_dz_phi2,vy_dz_phi1,vy_dz_phi2,
            DX,DY,DZ,DT);
		//添加检波器记录的
		for (int i = 0; i < ZPLANE_RECORD_XLEN; ++i)
        for(int j=0;j<ZPLANE_RECORD_YLEN;++j)
		{
			const size_t record_to_pos = (t * ZPLANE_RECORD_XLEN + i)*ZPLANE_RECORD_YLEN+j;
			const size_t record_from_pos = ((i*ZPLANE_XINTERVAL + OFFSET+ZPLANE_XBEGIN)*YDIM + ZPLANE_YBEGIN+j*ZPLANE_YINTERVAL+OFFSET)*ZDIM+
            OFFSET+ZPLANE_ZPOS;
			//sigma_xx_xline_geo[record_to_pos] = sigma_xx[record_from_pos];
            //sigma_zz_xline_geo[record_to_pos] = sigma_zz[record_from_pos];
            //tau_y_xline_geo[record_to_pos] = sigma_xx[record_from_pos];
            vx_xline_geo[record_to_pos]=vx[record_from_pos];
            vy_xline_geo[record_to_pos]=vy[record_from_pos];
            vz_xline_geo[record_to_pos]=vz[record_from_pos];
		}
	
	}
	TIME_END("El3D")
	//Write2DBinEx(XLINE_NAME.c_str(),"-sigma-xx-",sigma_xx_xline_geo,NT,RECORD_XLEN);
	//Write2DBinEx(XLINE_NAME.c_str(), "-sigma-zz-", sigma_zz_xline_geo, NT, RECORD_XLEN);
	//Write2DBinEx(XLINE_NAME.c_str(), "-tau-y-", tau_y_xline_geo, NT, RECORD_XLEN);
    Write3DBinEx(ZPLANE_NAME.c_str(), "-vx-", vx_xline_geo, NT, ZPLANE_RECORD_XLEN,ZPLANE_RECORD_YLEN);
    Write3DBinEx(ZPLANE_NAME.c_str(), "-vy-", vy_xline_geo, NT, ZPLANE_RECORD_XLEN,ZPLANE_RECORD_YLEN);
    Write3DBinEx(ZPLANE_NAME.c_str(), "-vz-", vz_xline_geo, NT, ZPLANE_RECORD_XLEN,ZPLANE_RECORD_YLEN);


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
    DebugFreeMemory(&vy_dx_phi1);
	DebugFreeMemory(&vy_dx_phi2);
    DebugFreeMemory(&vz_dx_phi1);
    DebugFreeMemory(&vz_dx_phi2);
    DebugFreeMemory(&vy_dy_phi1);
    DebugFreeMemory(&vy_dy_phi2);
    DebugFreeMemory(&vx_dy_phi1);
    DebugFreeMemory(&vx_dy_phi2);
    DebugFreeMemory(&vz_dy_phi1);
    DebugFreeMemory(&vz_dy_phi2);
    DebugFreeMemory(&vy_dz_phi1);
	DebugFreeMemory(&vy_dz_phi2);
    DebugFreeMemory(&vx_dz_phi1);
	DebugFreeMemory(&vx_dz_phi2);
	DebugFreeMemory(&vz_dz_phi1);
	DebugFreeMemory(&vz_dz_phi2);
	DebugFreeMemory(&sigma_x_dx_phi1);
	DebugFreeMemory(&sigma_x_dx_phi2);
    DebugFreeMemory(&tau_y_dx_phi1);
    DebugFreeMemory(&tau_y_dx_phi2);
    DebugFreeMemory(&tau_z_dx_phi1);
    DebugFreeMemory(&tau_z_dx_phi2);
    DebugFreeMemory(&sigma_y_dy_phi1);
	DebugFreeMemory(&sigma_y_dy_phi2);
    DebugFreeMemory(&tau_x_dy_phi1);
    DebugFreeMemory(&tau_x_dy_phi2);
    DebugFreeMemory(&tau_z_dy_phi1);
    DebugFreeMemory(&tau_z_dy_phi2);
	DebugFreeMemory(&sigma_z_dz_phi1);
	DebugFreeMemory(&sigma_z_dz_phi2);
    DebugFreeMemory(&tau_x_dz_phi1);
    DebugFreeMemory(&tau_x_dz_phi2);
    DebugFreeMemory(&tau_y_dz_phi1);
    DebugFreeMemory(&tau_y_dz_phi2);
	DebugFreeMemory(&dense_model);
	DebugFreeMemory(&vp_model);
    DebugFreeMemory(&vs_model);
	DebugFreeMemory(&source_buffer);
	DebugFreeMemory(&sigma_xx);
    DebugFreeMemory(&sigma_yy);
    DebugFreeMemory(&sigma_zz);
    DebugFreeMemory(&tau_x);
    DebugFreeMemory(&tau_y);
    DebugFreeMemory(&tau_z);
	DebugFreeMemory(&vx);
    DebugFreeMemory(&vy);
	DebugFreeMemory(&vz);
	DebugFreeMemory(&vx_xline_geo);
    DebugFreeMemory(&vy_xline_geo);
	DebugFreeMemory(&vz_xline_geo);
	return 0;
}

