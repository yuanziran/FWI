/**
 * @file acoustic_forward_2d.h
 * @author Zhenghong Guo (you@domain.com)
 * @brief 二维各向同性声波一阶速度应力方程交错网格正演模拟程序
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once  
#ifdef __cplusplus
extern "C"{
#endif
#include"fwi.h"


/**
 * @brief 二维声波一阶速度应力方程交错网格正演模拟质点震动速度计算应力
 * 
 * @param stress        压力波场
 * @param vx            质点震动速度的水平分量
 * @param vz            质点震动速度的垂直分量
 * @param dense         密度参数
 * @param vp            纵波速度
 * @param xdim          x维度
 * @param zdim          z维度
 * @param xbegin        进行有限差分的x起始维度
 * @param zbegin        进行有限差分的z起始维度
 * @param xlen          有限差分计算长度
 * @param zlen          有限差分计算长度
 * @param order         有限差分的阶数
 * @param coeff         有限差分系数
 * @param dx            空间采样间隔
 * @param dz            空间采样间隔
 * @param dt            时间采样间隔
 * @return 
 */
DLL_API void  AcSgVelocityToStress2D(
    float*  const           stress,
    const float*  const     vx,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const int               xdim,
    const int               zdim,
    const int               xbegin,
    const int               zbegin,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const float             dx,
    const float             dz,
    const float             dt
);

/**
 * @brief 从应力计算质点震动速度交错网格版本
 * 参数的物理意义同上
 * @param vx 
 * @param vz 
 * @param stress 
 * @param dense 
 * @param xdim 
 * @param zdim 
 * @param xbegin 
 * @param zbegin 
 * @param xlen 
 * @param zlen 
 * @param order 
 * @param coeff 
 * @param dx 
 * @param dz 
 * @param dt 
 * @return DLL_API 
 */
DLL_API void  AcSgStressToVelocity2D(
    float*  const           vx,
    float*  const           vz,
    const float* const      stress,
    const float* const      dense,
    const int               xdim,
    const int               zdim,
    const int               xbegin,
    const int               zbegin,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const float             dx,
    const float             dz,
    const float             dt
);


//********************************************卷积一致吸收边界条件****************************
DLL_API void  AcSgVelocityToStressCPML2D(
    float*  const           stress,
    const float*  const     vx,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const int               xdim,
    const int               zdim,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    float*       const      vx_dx_phi1,
    float*       const      vx_dx_phi2,
    float*       const      vz_dz_phi1,
    float*       const      vz_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);


DLL_API void  AcSgStressToVelocityCPML2D(
    float*  const           vx,
    float*  const           vz,
    const float* const      stress,
    const float* const      dense,
    const int               xdim,
    const int               zdim,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    float*       const      stress_dx_phi1,
    float*       const      stress_dx_phi2,
    float*       const      stress_dz_phi1,
    float*       const      stress_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);

//**************************************************吸收边界效果不好没有测试成功************************
DLL_API void  AcSgVelocityToStressCPML2DEx(
    float*  const           stress,
    const float*  const     vx,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const int               xdim,
    const int               zdim,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x1_cpml_kappa,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      x2_cpml_kappa,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z1_cpml_kappa,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    const float* const      z2_cpml_kappa,
    float*       const      vx_dx_phi1,
    float*       const      vx_dx_phi2,
    float*       const      vz_dz_phi1,
    float*       const      vz_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);


DLL_API void  AcSgStressToVelocityCPML2DEx(
    float*  const           vx,
    float*  const           vz,
    const float* const      stress,
    const float* const      dense,
    const int               xdim,
    const int               zdim,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x1_cpml_kappa,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      x2_cpml_kappa,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z1_cpml_kappa,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    const float* const      z2_cpml_kappa,
    float*       const      stress_dx_phi1,
    float*       const      stress_dx_phi2,
    float*       const      stress_dz_phi1,
    float*       const      stress_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);



#ifdef __cplusplus
}
#endif