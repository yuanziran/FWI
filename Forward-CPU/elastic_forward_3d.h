/**
 * @file elastic_forward_3d.h
 * @author Zhenghong Guo (you@domain.com)
 * @brief 三维弹性波一阶速度应力交错网格形式
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

DLL_API void  ElSgVelocityToStress3D(
    float*  const           sigma_xx,
    float*  const           sigma_yy,
    float*  const           sigma_zz,
    float*  const           tau_x,
    float*  const           tau_y,
    float*  const           tau_z,
    const float*  const     vx,
    const float*  const     vy,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const float* const      vs,
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

DLL_API void  ElSgStressToVelocity3D(
    float*  const                 vx,
    float*  const                 vy,
    float*  const                 vz,
    const float*  const           sigma_xx,
    const float*  const           sigma_yy,
    const float*  const           sigma_zz,
    const float*  const           tau_x,
    const float*  const           tau_y,
    const float*  const           tau_z,
    const float*  const           dense,
    const int                     xdim,
    const int                     ydim,
    const int                     zdim,
    const int                     xbegin,
    const int                     ybegin,
    const int                     zbegin,
    const int                     xlen,
    const int                     ylen,
    const int                     zlen,
    const int                     order,
    const float* const            coeff,
    const float                   dx,
    const float                   dy,
    const float                   dz,
    const float                   dt
);


//**************************************CPML吸收边界条件**********************
DLL_API void  ElSgVelocityToStressCPML3D(
    float*  const           sigma_xx,
    float*  const           sigma_yy,
    float*  const           sigma_zz,
    float*  const           tau_x,
    float*  const           tau_y,
    float*  const           tau_z,
    const float*  const     vx,
    const float*  const     vy,
    const float*  const     vz,
    const float* const      dense,
    const float* const      vp,
    const float* const      vs,
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
    float*       const      vy_dx_phi1,
    float*       const      vy_dx_phi2,
    float*       const      vz_dx_phi1,
    float*       const      vz_dx_phi2,
    float*       const      vy_dy_phi1,
    float*       const      vy_dy_phi2,
    float*       const      vx_dy_phi1,
    float*       const      vx_dy_phi2,
    float*       const      vz_dy_phi1,
    float*       const      vz_dy_phi2,
    float*       const      vz_dz_phi1,
    float*       const      vz_dz_phi2,
    float*       const      vx_dz_phi1,
    float*       const      vx_dz_phi2,
    float*       const      vy_dz_phi1,
    float*       const      vy_dz_phi2,
    const float             dx,
    const float             dy,
    const float             dz,
    const float             dt
);


DLL_API void  ElSgStressToVelocityCPML3D(
    float*  const                 vx,
    float*  const                 vy,
    float*  const                 vz,
    const float*  const           sigma_xx,
    const float*  const           sigma_yy,
    const float*  const           sigma_zz,
    const float*  const           tau_x,
    const float*  const           tau_y,
    const float*  const           tau_z,
    const float*  const           dense,
    const int                     xdim,
    const int                     ydim,
    const int                     zdim,
    const int                     xlen,
    const int                     ylen,
    const int                     zlen,
    const int                     order,
    const float* const            coeff,
    const int                     layers,
    const float* const            x1_cpml_alpha,
    const float* const            x1_cpml_beta,
    const float* const            x2_cpml_alpha,
    const float* const            x2_cpml_beta,
    const float* const            y1_cpml_alpha,
    const float* const            y1_cpml_beta,
    const float* const            y2_cpml_alpha,
    const float* const            y2_cpml_beta,
    const float* const            z1_cpml_alpha,
    const float* const            z1_cpml_beta,
    const float* const            z2_cpml_alpha,
    const float* const            z2_cpml_beta, 
    float*       const            sigma_xx_dx_phi1,
    float*       const            sigma_xx_dx_phi2,
    float*       const            tau_y_dx_phi1,
    float*       const            tau_y_dx_phi2,
    float*       const            tau_z_dx_phi1,
    float*       const            tau_z_dx_phi2,
    float*       const            sigma_yy_dy_phi1,
    float*       const            sigma_yy_dy_phi2,
    float*       const            tau_x_dy_phi1,
    float*       const            tau_x_dy_phi2,
    float*       const            tau_z_dy_phi1,
    float*       const            tau_z_dy_phi2,
    float*       const            sigma_zz_dz_phi1,
    float*       const            sigma_zz_dz_phi2,
    float*       const            tau_x_dz_phi1,
    float*       const            tau_x_dz_phi2,
    float*       const            tau_y_dz_phi1,
    float*       const            tau_y_dz_phi2,
    const float                   dx,
    const float                   dy,
    const float                   dz,
    const float                   dt
);

#ifdef __cplusplus
}
#endif