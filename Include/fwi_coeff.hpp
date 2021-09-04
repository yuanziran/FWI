/**
 * @file fwi_coeff.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 有限差分系数和卷积一致吸收边界条件
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once 
#include<math.h>

//交错网格一阶有限差分系数
static inline void SgCdoOp1(float* coeffs, const unsigned int size)
{
    //std::cout << "**********SgCdoCoeff************************" << std::endl;
    float sign_coeff = 1.0f;
    for (unsigned int m = 1; m <= size; m++) {
        float c = sign_coeff / (m + m - 1.0f);
        float m2 = (m + m - 1.0f) * (m + m - 1.0f);
        for (unsigned int i = 1; i <= size; i++) {
            if (i == m) continue;
            float i2 = (i + i - 1.0f) * (i + i - 1.0f);
            c *= i2 / fabs(i2 - m2);
        }
        coeffs[m - 1] = c;
        sign_coeff *= -1.0;
        //std::cout << "coeff [" << m << "] =" << coeffs[m - 1] << std::endl;
    }
}

//常规网格一阶有限差分
static inline void RgCdoOp1(float* coeffs,const unsigned int size)
{
    float sign_coeff = 1.0;
    for (unsigned int m=1; m<=size; m++) {
        float c = sign_coeff / (m + m);
        float m2 = (float)(m) * m;
        for (unsigned int i=1; i<=size; i++) {
            if (i == m) continue;
            float i2 = (float)(i) * i;
            c *= i2 / fabs(i2 - m2);
        }
        coeffs[m-1] = c;
        sign_coeff *= -1.0;
    }
}

//常规网格二阶有限差分系数
static inline void RgCdoOp2(float* const coeffs, const unsigned int size)
{
    unsigned int radius = size - 1;
    coeffs[0] = 0.0f;
    for (unsigned int i = 1; i <= radius; ++i)
    {
        float i2 = float(i * i);
        float a = 1.0f / i2;
        if (i % 2 == 0) a = -a;
        for (unsigned int j = 1; j <= radius; ++j)
        {
            if (i != j)
            {
                a *= fabsf(1.0f / (1.0f - i2 / float(j * j)));
            }
        }
        coeffs[i] = a;
        coeffs[0] -= 2.0f * coeffs[i];
    }
}


//卷积一致吸收边界条件
static inline void Cpmlfun(float* const cpmlAlpha, float* const cpmlBeta, const int abcLayers, const float dh, const float dt, const float maxFreq, const float maxVelocity, const float dampingRate, const int m = 3)
{
    const float d_max = -(m + 1) / (2 * abcLayers * dh) * maxVelocity * logf(dampingRate);
    const float alpha_max = 3.1415926f * maxFreq;
    for (int i = 0; i < abcLayers; ++i)
    {
        float frac = float(i) / abcLayers;
        float alpha = alpha_max * (1.0f - frac);
        if (alpha < 0.0f) alpha = 0.0f;
        float d = d_max * powf(frac, float(m));
        cpmlAlpha[abcLayers - i - 1] = expf(-(d + alpha) * dt);
        cpmlBeta[abcLayers - i - 1] = d / ((d + alpha)) * (cpmlAlpha[abcLayers - i - 1] - 1.0f);
    }
}



//卷积一致吸收边界条件改进版本
//该代码没有测试通过
//吸收边界效果不好
static inline void CpmlfunEx(float* const cpmlAlpha, float* const  cpmlBeta,float* const kappa, const int abcLayers, const float dh, const float dt, const float maxFreq, const float maxVelocity, const float dampingRate, const int m = 3)
{
    const float d_max = -(m + 1) / (2 * abcLayers * dh) * maxVelocity * logf(dampingRate);
    const float alpha_max = 3.1415926f * maxFreq;
    for (int i = 0; i < abcLayers; ++i)
    {
        float frac = float(i) / abcLayers;
        float alpha = alpha_max * (1.0f - frac);
        kappa[abcLayers-i-1]=1+frac;
        if (alpha < 0.0f) alpha = 0.0f;
        float d = d_max * powf(frac, float(m));
        cpmlAlpha[abcLayers - i - 1] = expf(-(d/kappa[abcLayers-i-1] + alpha) * dt);
        cpmlBeta[abcLayers - i - 1] = d / (kappa[abcLayers-i-1]*(d + kappa[abcLayers-i-1]*alpha)) * (cpmlAlpha[abcLayers - i - 1] - 1.0f);
		kappa[abcLayers - i - 1] = 1.0f / kappa[abcLayers - i - 1];
	}
}


//海绵吸收边界条件
static inline void SponealFun(float* const sponeal,const int layers)
{
    //const float scalar=1.0f/(expf(1.0f)-1);
    for(int i=layers;i>0;--i)
    {
        const float dh=(layers-i)/(float)layers;
        //sponeal[layers-i]=(expf(-dh*dh)-1.0f)*scalar;
        sponeal[layers-i]=expf(-0.7f*dh*dh);
        std::cout<<sponeal[layers-i]<<std::endl;
    }
}


static void inline ApplySponealAbc2D(float* const wave,const int xdim,const int zdim,
const int order,const int layers,const float* const abc_coeff)
{
    const int xlen=xdim-order-order;
    const int zlen=zdim-order-order;
    #pragma omp parallel for collapse(2)
    for(int i=0;i<layers;++i)
    for(int j=0;j<zlen;++j)
    {
        const size_t x_offset_1=(order+i)*zdim+j+order;
        const size_t x_offset_2=(xdim-order-layers+i)*zdim+j+order;
        wave[x_offset_1]*=abc_coeff[layers-i-1];
        wave[x_offset_2]*=abc_coeff[i];
    }
    #pragma omp parallel for collapse(2)
    for(int i=0;i<layers;++i)
    for(int j=0;j<xlen;++j)
    {
        const size_t z_offset_1=(j+order)*zdim+i+order;
        const size_t z_offset_2=(j+order)*zdim+zdim-layers+i-order;
        wave[z_offset_1]*=abc_coeff[layers-i-1];
        wave[z_offset_2]*=abc_coeff[i];
    }
}


static void inline ApplySponealAbc2DEx(
    float* const wave,
    float* const wave_dx,
    float* const wave_dz,
    const int xdim,const int zdim,
    const int order,const int layers,const float* const abc_coeff)
{
    const int xlen=xdim-order-order;
    const int zlen=zdim-order-order;
    #pragma omp parallel for collapse(2)
    for(int i=0;i<layers;++i)
    for(int j=0;j<zlen;++j)
    {
        const size_t x_offset_1=(order+i)*zdim+j+order;
        const size_t x_offset_2=(xdim-order-layers+i)*zdim+j+order;
        wave[x_offset_1]*=abc_coeff[i];
        wave[x_offset_2]*=abc_coeff[i];
        wave_dx[x_offset_1]*=abc_coeff[i];
        wave_dx[x_offset_2]*=abc_coeff[i];
        wave_dz[x_offset_1]*=abc_coeff[i];
        wave_dz[x_offset_2]*=abc_coeff[i];
    }
    #pragma omp parallel for collapse(2)
    for(int i=0;i<layers;++i)
    for(int j=0;j<xlen;++j)
    {
        const size_t z_offset_1=(j+order)*zdim+i+order;
        const size_t z_offset_2=(j+order)*zdim+zdim-layers+i-order;
        wave[z_offset_1]*=abc_coeff[i];
        wave[z_offset_2]*=abc_coeff[i];
        wave_dx[z_offset_1]*=abc_coeff[i];
        wave_dx[z_offset_2]*=abc_coeff[i];
        wave_dz[z_offset_1]*=abc_coeff[i];
        wave_dz[z_offset_2]*=abc_coeff[i];
    }
}
