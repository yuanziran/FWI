/**
 * @file elastic_forward_2d.h
 * @author Zhenghong Guo (you@domain.com)
 * @brief  弹性波一阶速度应力方程交错网格形式
 * @version 0.1
 * @date 2021-06-28
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

DLL_API void  ElSgVelocityToStress2D(
    float*  const           sigma_xx,
    float*  const           sigma_zz,
    float*  const           tau_y,
    const float*  const     vx,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const float* const      vs,
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

DLL_API void  ElSgStressToVelocity2D(
    float*  const           vx,
    float*  const           vz,
    float*  const           sigma_xx,
    float*  const           sigma_zz,
    float*  const           tau_y,
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


//********************************************CPML吸收边界条件*************************
DLL_API void  ElSgVelocityToStressCPML2D(
    float*  const           sigma_xx,
    float*  const           sigma_zz,
    float*  const           tau_y,
    const float*  const     vx,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const float* const      vs,
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
    float*       const      vz_dx_phi1,
    float*       const      vz_dx_phi2,
    float*       const      vx_dz_phi1,
    float*       const      vx_dz_phi2,
    float*       const      vz_dz_phi1,
    float*       const      vz_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);

DLL_API void  ElSgStressToVelocityCPML2D(
    float*  const           vx,
    float*  const           vz,
    float*  const           sigma_xx,
    float*  const           sigma_zz,
    float*  const           tau_y,
    const float* const      dense,
    const int               xdim,
    const int               zdim,
    const int               xlen,
    const int               zlen,
    const int               order,
    const float* const      coeff,
     const int              layers,
    const float* const      x1_cpml_alpha,
    const float* const      x1_cpml_beta,
    const float* const      x2_cpml_alpha,
    const float* const      x2_cpml_beta,
    const float* const      z1_cpml_alpha,
    const float* const      z1_cpml_beta,
    const float* const      z2_cpml_alpha,
    const float* const      z2_cpml_beta,
    float*       const      sigma_xx_dx_phi1,
    float*       const      sigma_xx_dx_phi2,
    float*       const      tau_y_dx_phi1,
    float*       const      tau_y_dx_phi2,
    float*       const      sigma_zz_dz_phi1,
    float*       const      sigma_zz_dz_phi2,
    float*       const      tau_y_dz_phi1,
    float*       const      tau_y_dz_phi2,
    const float             dx,
    const float             dz,
    const float             dt
);




#ifdef __cplusplus
}
#endif