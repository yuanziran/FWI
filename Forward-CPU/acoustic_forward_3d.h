/**
 * @file acoustic_forward_3d.h
 * @author Zhenghong Guo (you@domain.com)
 * @brief 三维各向同性声波一阶速度应力方程交错网格正演模拟程序
 * @version 0.1
 * @date 2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once  
#ifdef __cplusplus
extern "C"{
#endif
#include"fwi.h"
#include<stddef.h>

DLL_API void  AcSgVelocityToStress3D(
    float*  const           stress,
    const float*  const     vx,
    const float*  const     vy,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const int               xdim,
    const int               ydim,
    const int               zdim,
    const int               xbegin,
    const int               ybegin,
    const int               zbegin,
    const int               xlen,
    const int               ylen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const float             dx,
    const float             dy,
    const float             dz,
    const float             dt
);


DLL_API void  AcSgStressToVelocity3D(
    float*  const           vx,
    float*  const           vy,
    float*  const           vz,
    const float* const      stress,
    const float* const      dense,
    const int               xdim,
    const int               ydim,
    const int               zdim,
    const int               xbegin,
    const int               ybegin,
    const int               zbegin,
    const int               xlen,
    const int               ylen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const float             dx,
    const float             dy,
    const float             dz,
    const float             dt
);


//********************************************************************3维度卷积一致吸收边界边界CPML*****************************
 DLL_API void  AcSgVelocityToStressCPML3D(
    float*  const           stress,
    const float*  const     vx,
    const float*  const     vy,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const int               xdim,
    const int               ydim,
    const int               zdim,
    const int               xlen,
    const int               ylen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      y1_cpml_alpha,
    const float* const      y1_cpml_beta,
    const float* const      y2_cpml_alpha,
    const float* const      y2_cpml_beta,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    float*       const      vx_dx_phi1,
    float*       const      vx_dx_phi2,
    float*       const      vy_dy_phi1,
    float*       const      vy_dy_phi2,
    float*       const      vz_dz_phi1,
    float*       const      vz_dz_phi2,
    const float             dx,
    const float             dy,
    const float             dz,
    const float             dt
);


DLL_API void  AcSgStressToVelocityCPML3D(
    float*  const           vx,
    float*  const           vy,
    float*  const           vz,
    const float* const      stress,
    const float* const      dense,
    const int               xdim,
    const int               ydim,
    const int               zdim,
    const int               xlen,
    const int               ylen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
    const int               layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      y1_cpml_alpha,
    const float* const      y1_cpml_beta,
    const float* const      y2_cpml_alpha,
    const float* const      y2_cpml_beta,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    float*       const      stress_dx_phi1,
    float*       const      stress_dx_phi2,
    float*       const      stress_dy_phi1,
    float*       const      stress_dy_phi2,
    float*       const      stress_dz_phi1,
    float*       const      stress_dz_phi2,
    const float             dx,
    const float             dy,
    const float             dz,
    const float             dt
);



#ifdef __cplusplus
}
#endif
