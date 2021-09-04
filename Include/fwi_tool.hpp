/**
 * @file fwi_tool.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 一些有用的工具
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once 
#include<assert.h>
#include<iostream>
#include<math.h>
//扩充模型的边界
static inline void Expand2D(float* const array,const float* const kernel,const int full_xdim,const int full_zdim,const int naked_xdim,const int naked_zdim,
const int layers,const int methods=0)
{
    assert((full_xdim==(naked_xdim+2*layers)));
    assert((full_zdim==(naked_zdim+2*layers)));
    #pragma omp parallel for collapse(2)
    for(int i=0;i<naked_xdim;++i)
    for(int j=0;j<naked_zdim;++j)
    {
        const size_t full_pos=(i+layers)*full_zdim+j+layers;
        const size_t nake_pos=i*naked_zdim+j;
        array[full_pos]=kernel[nake_pos];
    }
    #pragma omp parallel for collapse(2) 
    for(int i=0;i<layers;++i)
    for(int j=0;j<naked_zdim;++j)
    {
        array[(i*full_zdim)+j+layers]=array[(layers*full_zdim)+j+layers];
        array[((i+layers+naked_xdim)*full_zdim+j+layers)]=array[(full_xdim-layers-1)*full_zdim+j+layers];
    }
    #pragma omp parallel for collapse(2)
    for(int i=0;i<layers;++i)
    for(int j=0;j<full_xdim;++j)
    {
        array[j*full_zdim+i]=array[j*full_zdim+layers];
        array[j*full_zdim+i+naked_zdim+layers]=array[j*full_zdim+full_zdim-layers-1];
    }
}

//拉普拉斯滤波算子
static inline void LapacianFilter2D(const float* const  in,float* const out,const float* const coeff,const int full_xdim,const int full_zdim,
const int naked_xdim,const int naked_zdim,const int order,const float dx,const float dz)
{
    const float dx_2_inv=1/(dx*dx);
    const float dz_2_inv=1/(dz*dz);
    #pragma omp parallel for collapse(2)
    for(int i=0;i<naked_xdim;++i)
    for(int j=0;j<naked_zdim;++j)
    {
        const size_t full_pos=(i+order)*full_zdim+j+order;
        const size_t nake_pos=i*naked_zdim+j;
        out[nake_pos]=coeff[0]*in[full_pos]*(dx_2_inv+dz_2_inv);
        for(int k=1;k<=order;++k)
        {
            out[nake_pos]+=coeff[k]*((in[full_pos+k*full_zdim]+in[full_pos-k*full_zdim])*dx_2_inv+(in[full_pos+k]+in[full_pos-k])*dz_2_inv);
        }
    }
}


static inline float GetAverageModelRect(const float* const data, const int xdim, const int zdim, const int xbegin, const int zbegin, const int xlen, const int zlen)
{
	const float scalar = 1.0f / (xlen*zlen);
	float average = 0.0f;
	for(int i=0;i<xlen;++i)
		for (int j = 0; j < zlen; ++j)
		{
			average += data[(i+xbegin)*zdim+j+zbegin];
		}
	return average * scalar;
}


static inline void GetCpmlAverageVpModelRect(float vp_info[], const float* const data, const int xdim, const int zdim, const int cdo_order, const int abc_layers)
{
	constexpr int x1_vp = 0;
	constexpr int x2_vp = 1;
	constexpr int z1_vp = 2;
	constexpr int z2_vp = 3;
	const int xlen = xdim - 2 * cdo_order;
	const int zlen = zdim - 2 * cdo_order;
	vp_info[z1_vp] = GetAverageModelRect(data,xdim,zdim,cdo_order,cdo_order,xlen,abc_layers);
	vp_info[x1_vp]= GetAverageModelRect(data, xdim, zdim, cdo_order, cdo_order,abc_layers,zlen);
	vp_info[z2_vp] = GetAverageModelRect(data, xdim, zdim,cdo_order, zlen - abc_layers, xlen, abc_layers);
	vp_info[x2_vp] = GetAverageModelRect(data, xdim, zdim, xlen - abc_layers,cdo_order, abc_layers, zlen);
	std::cout << "cpml average vp info : " << std::endl
		<< "vp  x1  : " << vp_info[x1_vp] << std::endl
		<< "vp  x2  : " << vp_info[x2_vp] << std::endl
		<< "vp  z1  : " << vp_info[z1_vp] << std::endl
		<< "vp  z2  : " << vp_info[z2_vp] << std::endl;
}


//**********************************************3维操作*****************************************
static inline void Expand3D(float* const array, const float* const kernel, const int full_xdim, const int full_ydim,const int full_zdim, const int naked_xdim,const int naked_ydim, const int naked_zdim,
	const int layers, const int methods = 0)
{
	assert((full_xdim == (naked_xdim + 2 * layers)));
	assert((full_ydim == (naked_ydim + 2 * layers)));
	assert((full_zdim == (naked_zdim + 2 * layers)));

	#pragma omp parallel for collapse(3)
	for (int i = 0; i < naked_xdim; ++i)
		for (int j = 0; j < naked_ydim; ++j)
			for(int k=0; k<naked_zdim;++k)
			{
				const size_t full_pos = ((i + layers)*full_ydim + j + layers)*full_zdim+k+layers;
				const size_t nake_pos = (i * naked_ydim + j)*naked_zdim+k;
				array[full_pos] = kernel[nake_pos];
			}
	//填充X方向的边界
	#pragma omp parallel for collapse(3) 
	for (int i = 0; i < layers; ++i)
		for (int j = 0; j < naked_ydim; ++j)
			for(int k=0;k<naked_zdim;++k)
			{
				array[((i*full_ydim) + j + layers)*full_zdim+k+layers] = array[((layers*full_ydim) + j + layers)*full_zdim+k+layers];
				array[(((i + layers + naked_xdim)*full_ydim + j + layers)*full_zdim+k+layers)] = array[((full_xdim - layers - 1)*full_ydim + j + layers)*full_zdim+k+layers];
			}
	//填充Y方向的边界
	#pragma omp parallel for collapse(3)
	for (int i = 0; i < layers; ++i)
		for (int j = 0; j < full_xdim; ++j)
			for(int k=0;k<naked_zdim;++k)
			{
				array[j*full_ydim*full_zdim + i*full_zdim+k+layers] = array[j*full_ydim*full_zdim + layers*full_zdim+ k + layers];
				array[j*full_ydim*full_zdim + (i + naked_ydim + layers)*full_zdim+k+layers] = array[j*full_ydim*full_zdim + (full_ydim - layers - 1)*full_zdim+k+layers];
			}
	//填充Z方向的边界
	#pragma omp parallel for collapse(3)
	for(int i=0;i<layers;++i)
		for(int j=0;j<full_xdim;++j)
			for (int k = 0; k < full_ydim; ++k)
			{
				array[(j*full_ydim + k)*full_zdim + i] = array[(j*full_ydim + k)*full_zdim + layers];
				array[(j*full_ydim + k)*full_zdim + i + layers + naked_zdim] = array[(j*full_ydim + k)*full_zdim + full_zdim - layers - 1];
			}
}

//3维拉普拉斯滤波
static inline void LapacianFilter3D(const float* const  in, float* const out, const float* const coeff, const int full_xdim,const int full_ydim, const int full_zdim,
	const int naked_xdim,const int naked_ydim, const int naked_zdim, const int order, const float dx,const float dy, const float dz)
{
	const float dx_2_inv = 1 / (dx*dx);
	const float dy_2_inv = 1 / (dy*dy);
	const float dz_2_inv = 1 / (dz*dz);
	#pragma omp parallel for collapse(3)
	for (int i = 0; i < naked_xdim; ++i)
		for (int j = 0; j < naked_ydim; ++j)
			for(int k=0;k<naked_zdim;++k)
			{
				const size_t full_pos = ((i + order)*full_ydim + j + order)*full_zdim+k+order;
				const size_t nake_pos = (i * naked_ydim + j)*naked_zdim+k;
				out[nake_pos] = coeff[0] * in[full_pos] * (dx_2_inv + dy_2_inv+dz_2_inv);
				for (int ll = 1; ll <= order; ++k)
				{
					out[nake_pos] += coeff[ll] * (
					(in[full_pos + ll *full_ydim*full_zdim] + in[full_pos - k *full_ydim* full_zdim])*dx_2_inv 
					+ (in[full_pos+ll*full_zdim]+in[full_pos-ll*full_zdim])*dy_2_inv+
					(in[full_pos + k] + in[full_pos - k])*dz_2_inv
					);
				}
		}
}


static inline float GetAverageModelBox(const float* const data, const int xdim,const int ydim, const int zdim, const int xbegin,const int ybegin ,const int zbegin, const int xlen,const int ylen,const int zlen)
{
	const float scalar = 1.0f / (xlen*ylen*zlen);
	float average = 0.0f;
	for (int i = 0; i < xlen; ++i)
		for (int j = 0; j < ylen; ++j)
			for(int k=0;k<zlen;++k)
			{
				average += data[((i + xbegin)*ydim + j + ybegin)*zdim+k+zbegin];
			}
	return average * scalar;
}


static inline void GetCpmlAverageVpModelBox(float vp_info[], const float* const data, const int xdim,const int ydim ,const int zdim, const int cdo_order, const int abc_layers)
{
	constexpr int x1_vp = 0;
	constexpr int x2_vp = 1;
	constexpr int y1_vp = 2;
	constexpr int y2_vp = 3;
	constexpr int z1_vp = 4;
	constexpr int z2_vp = 5;
	const int xlen = xdim - 2 * cdo_order;
	const int ylen = ydim - 2 * cdo_order;
	const int zlen = zdim - 2 * cdo_order;
	vp_info[x1_vp] = GetAverageModelBox(data,xdim,ydim,zdim,cdo_order,cdo_order,cdo_order,abc_layers,ylen,zlen);
	vp_info[y1_vp] = GetAverageModelBox(data, xdim, ydim, zdim, cdo_order, cdo_order, cdo_order, xlen, abc_layers, zlen);
	vp_info[z1_vp] = GetAverageModelBox(data, xdim, ydim, zdim, cdo_order, cdo_order, cdo_order, xlen, ylen, abc_layers);
	vp_info[x2_vp] = GetAverageModelBox(data, xdim, ydim, zdim, xlen-abc_layers, cdo_order, cdo_order, abc_layers, ylen, zlen);
	vp_info[y2_vp] = GetAverageModelBox(data, xdim, ydim, zdim, cdo_order, ylen-abc_layers, cdo_order, xlen, abc_layers, zlen);
	vp_info[z2_vp] = GetAverageModelBox(data, xdim, ydim, zdim, cdo_order, cdo_order,zlen-abc_layers, xlen, ylen, abc_layers);
	//vp_info[z1_vp] = GetAverageModelRect(data, xdim, zdim, cdo_order, cdo_order, xlen, abc_layers);
	//vp_info[x1_vp] = GetAverageModelRect(data, xdim, zdim, cdo_order, cdo_order, abc_layers, zlen);
	//vp_info[z2_vp] = GetAverageModelRect(data, xdim, zdim, cdo_order, zlen - abc_layers, xlen, abc_layers);
	//vp_info[x2_vp] = GetAverageModelRect(data, xdim, zdim, xlen - abc_layers, cdo_order, abc_layers, zlen);
	std::cout << "cpml average vp info : " << std::endl
		<< "vp  x1  : " << vp_info[x1_vp] << std::endl
		<< "vp  x2  : " << vp_info[x2_vp] << std::endl
		<< "vp  y1  : " << vp_info[y1_vp] << std::endl
		<< "vp  y2  : " << vp_info[y2_vp] << std::endl
		<< "vp  z1  : " << vp_info[z1_vp] << std::endl
		<< "vp  z2  : " << vp_info[z2_vp] << std::endl;
}


static inline int GetLen(const int begin,const int interval,const int end)
{
	return (end-1-begin)/interval+1;
}



