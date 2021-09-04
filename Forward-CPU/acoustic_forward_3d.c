/**
 * @file acoustic_forward_3d.c
 * @author Zhenghong Guo (you@domain.com)
 * @brief 三维声波各向同性波动方程交错网格
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include"acoustic_forward_3d.h"

void  AcSgVelocityToStress3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+xbegin)*ydim+jj+ybegin)*zdim+kk*zbegin;
        float _vxx=0.0f;
        float _vyy=0.0f;
        float _vzz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _vxx+=coeff[ll]*(vx[pos+ll*yzdim]-vx[pos-(ll+1)*yzdim]);
            _vyy+=coeff[ll]*(vy[pos+ll*zdim]-vy[pos-(ll+1)*zdim]);
            _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
        }
        _vxx/=dx;_vyy/=dy;_vzz/=dz;
        stress[pos]+=dt*dense[pos]*vp[pos]*vp[pos]*(_vxx+_vyy+_vzz);
    }
}

void  AcSgStressToVelocity3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+xbegin)*ydim+jj+ybegin)*zdim+kk*zbegin;
        float _px=0.0f;
        float _py=0.0f;
        float _pz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _px+=coeff[ll]*(stress[pos+(ll+1)*yzdim]-stress[pos-ll*yzdim]);
            _py+=coeff[ll]*(stress[pos+(ll+1)*zdim]-stress[pos-ll*zdim]);
            _pz+=coeff[ll]*(stress[pos+(ll+1)]-stress[pos-ll]);
        }
        _px/=dx;_py/=dy;_pz/=dz;
        vx[pos]+=dt*2.0f/(dense[pos]+dense[pos+yzdim])*_px;
        vy[pos]+=dt*2.0f/(dense[pos]+dense[pos+zdim])*_py;
        vz[pos]+=dt*2.0f/(dense[pos]+dense[pos+1])*_pz;
    }
}


void  AcSgVelocityToStressCPML3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+order)*ydim+jj+order)*zdim+kk*order;
        float _vxx=0.0f;
        float _vyy=0.0f;
        float _vzz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _vxx+=coeff[ll]*(vx[pos+ll*yzdim]-vx[pos-(ll+1)*yzdim]);
            _vyy+=coeff[ll]*(vy[pos+ll*zdim]-vy[pos-(ll+1)*zdim]);
            _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
        }
        _vxx/=dx;_vyy/=dy;_vzz/=dz;
        const float scalar=dt*dense[pos]*vp[pos]*vp[pos];
        stress[pos]+=scalar*(_vxx+_vyy+_vzz);
        //进入吸收边界
        if(ii<layers)//进入x方向的吸收边界
        {
            //[layers][ylen][zlen]
            const size_t cpml_pos=(ii*ylen+jj)*zlen+kk;
            stress[pos]+=scalar*vx_dx_phi1[cpml_pos];
            vx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vxx;
        }
        else if(ii>(xlen-layers-1))
        {
            const size_t cpml_pos=((xlen-ii-1)*ylen+jj)*zlen+kk;
            stress[pos]+=scalar*vx_dx_phi2[cpml_pos];
            vx_dx_phi2[cpml_pos]=vx_dx_phi2[cpml_pos]*x2_cpml_alpha[xlen-ii-1]+x2_cpml_beta[xlen-ii-1]*_vxx;
        }

        //进入y方向的吸收边界
        if(jj<layers)
        {
            //[layers][xlen][zlen]
            const size_t cpml_pos=(jj*xlen+ii)*zlen+kk;
            stress[pos]+=scalar*vy_dy_phi1[cpml_pos];
            vy_dy_phi1[cpml_pos]=vy_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_vyy;
            
        }
        else if(jj>(ylen-layers-1))
        {
            const size_t cpml_pos=((ylen-1-jj)*xlen+ii)*zlen+kk;
            stress[pos]+=scalar*vy_dy_phi2[cpml_pos];
            vy_dy_phi2[cpml_pos]=vy_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-1-jj]+y2_cpml_beta[ylen-1-jj]*_vyy;
        }
        //进入Z方向的吸收边界
        if(kk<layers)
        {
            //[layers][xlen][ylen]
            const size_t cpml_pos=(kk*xlen+ii)*ylen+jj;
            stress[pos]+=scalar*vz_dz_phi1[cpml_pos];
            vz_dz_phi1[cpml_pos]=vz_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_vzz;
        }
        else if(kk>(zlen-layers-1))
        {
            const size_t cpml_pos=((zlen-kk-1)*xlen+ii)*ylen+jj;
            stress[pos]+=scalar*vz_dz_phi2[cpml_pos];
            vz_dz_phi2[cpml_pos]=vz_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_vzz;
        }
    }
}

void  AcSgStressToVelocityCPML3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+order)*ydim+jj+order)*zdim+kk*order;
        float _px=0.0f;
        float _py=0.0f;
        float _pz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _px+=coeff[ll]*(stress[pos+(ll+1)*yzdim]-stress[pos-ll*yzdim]);
            _py+=coeff[ll]*(stress[pos+(ll+1)*zdim]-stress[pos-ll*zdim]);
            _pz+=coeff[ll]*(stress[pos+(ll+1)]-stress[pos-ll]);
        }
        _px/=dx;_py/=dy;_pz/=dz;
        float x_scalar=dt*2.0f/(dense[pos]+dense[pos+ydim*zdim]);
        float y_scalar=dt*2.0f/(dense[pos]+dense[pos+zdim]);
        float z_scalar=dt*2.0f/(dense[pos]+dense[pos+1]);
        vx[pos]+=x_scalar*_px;
        vy[pos]+=y_scalar*_py;
        vz[pos]+=z_scalar*_pz;
        if(ii<layers)
        {
            const size_t cpml_pos=(ii*ylen+jj)*zlen+kk;
            vx[pos]+=x_scalar*stress_dx_phi1[cpml_pos];
            stress_dx_phi1[cpml_pos]=stress_dx_phi1[cpml_pos]*x1_cpml_alpha[ii]+z1_cpml_beta[ii]*_px;
        }
        else if(ii>(xlen-layers-1))
        {
            const size_t cpml_pos=((xlen-ii-1)*ylen+jj)*zlen+kk;
            vx[pos]+=x_scalar*stress_dx_phi2[cpml_pos];
            stress_dx_phi2[cpml_pos]=stress_dx_phi2[cpml_pos]*x2_cpml_alpha[xlen-ii-1]+x2_cpml_beta[xlen-ii-1]*_px;
        }

        if(jj<layers)
        {
             const size_t cpml_pos=(jj*xlen+ii)*zlen+kk;
             vy[pos]+=y_scalar*stress_dy_phi1[cpml_pos];
             stress_dy_phi1[cpml_pos]=stress_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y2_cpml_beta[jj]*_py;
        }
        else if(jj>(ylen-layers-1))
        {
            const size_t cpml_pos=((ylen-1-jj)*xlen+ii)*zlen+kk;
            vy[pos]+=y_scalar*stress_dy_phi2[cpml_pos];
            stress_dy_phi2[cpml_pos]=stress_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-1-jj]+y2_cpml_beta[ylen-1-jj]*_py;
        }

        if(kk<layers)
        {
             const size_t cpml_pos=(kk*xlen+ii)*ylen+jj;
             vz[pos]+=z_scalar*stress_dz_phi1[cpml_pos];
             stress_dz_phi1[cpml_pos]=stress_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z2_cpml_beta[kk]*_pz;
        }
        else if(kk>(zlen-layers-1))
        {
            const size_t cpml_pos=((zlen-kk-1)*xlen+ii)*ylen+jj;
            vz[pos]+=z_scalar*stress_dz_phi2[cpml_pos];
            stress_dz_phi2[cpml_pos]=stress_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_pz;
        }
    }
}


