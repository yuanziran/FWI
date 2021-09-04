/**
 * @file fwi_load_point_source.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 常用的点震源加载方式
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once 
/**
 * @brief 点震源加载
 * 
 */
class LoadPointSourceOp
{
    public:
    //应力震源加载
    static inline void LoadAcStressSource(float* const stress,const int xdim,const int zdim,const int xpos,const int zpos,const float value)
    {
        stress[(xpos*zdim+zpos)]+=value;
    }

    static inline void LoadElStressSource(float* const sigma_xx,float* const sigma_zz,float* const tau_y,const int xdim,const int zdim,const int xpos,const int zpos,const float value)
    {
        const size_t global_pos=(xpos*zdim+zpos);
        sigma_xx[global_pos]+=value;
        sigma_zz[global_pos]+=value;
        tau_y[global_pos]+=value;
    }

    //加载P波震源质点震动速度震源
    static inline void LoadPSource(float* const vx,float* const vz,const int xdim,const int zdim,const int xpos,const int zpos,const float value)
    {
        const size_t global_pos=(xpos*zdim+zpos); 
        vx[global_pos+zdim]+=value;
        vx[global_pos+zdim+1]+=value;
        vx[global_pos-zdim]-=value;
        vx[global_pos-zdim+1]-=value;
        vz[global_pos+1]+=value;
        vz[global_pos+1-zdim]+=value;
        vz[global_pos-1]-=value;
        vz[global_pos-1-zdim]-=value;
    }

    //加载S波震源
    static inline void LoadSSource(float* const vx,float* const vz,const int xdim,const int zdim,const int xpos,const int zpos,const float value)
    {
        const size_t global_pos=(xpos*zdim+zpos);
        vx[global_pos+1]+=value;
        vx[global_pos+1-zdim]+=value;
        vx[global_pos-1]-=value;
        vx[global_pos-1-zdim]-=value;
        vz[global_pos+zdim]+=value;
        vz[global_pos+zdim-1]+=value;
        vz[global_pos-zdim]-=value;
        vz[global_pos-zdim-1]-=value;
    }

     //加载爆炸震源
    static inline void LoadEnergySource(float* const vx,float* const vz,const int xdim,const int zdim,const int xpos,const int zpos,const float value)
    {
        const size_t global_pos=(xpos*zdim+zpos);
        vx[global_pos]+=value;
        vx[global_pos+1]+=value;
        vx[global_pos+zdim]+=value;
        vx[global_pos+zdim+1]+=value;
        vz[global_pos]+=value;
        vz[global_pos+1]+=value;
        vz[global_pos+zdim]+=value;
        vz[global_pos+zdim+1]+=value;
    }

	//加载三维声波应力震源
	static inline void LoadAcStressSource3D(float* const stress, const int xdim,const int ydim, const int zdim, const int xpos,const int ypos, const int zpos, const float value)
	{
		stress[(xpos*ydim + ypos)*zdim+zpos] += value;
	}

};