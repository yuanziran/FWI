#include"elastic_forward_3d.h"

 void  ElSgVelocityToStress3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+xbegin)*ydim+jj+ybegin)*zdim+kk+zbegin;
        float _vxx=0.0f;
        float _vxy=0.0f;
        float _vxz=0.0f;
        float _vyx=0.0f;
        float _vyy=0.0f;
        float _vyz=0.0f;
        float _vzx=0.0f;
        float _vzy=0.0f;
        float _vzz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _vxx+=coeff[ll]*(vx[pos+ll*yzdim]-vx[pos-(ll+1)*yzdim]);
            _vyx+=coeff[ll]*(vy[pos+(ll+1)*yzdim]-vy[pos-ll*yzdim]);
            _vzx+=coeff[ll]*(vz[pos+(ll+1)*yzdim]-vz[pos-ll*yzdim]);

            _vxy+=coeff[ll]*(vx[pos+(ll+1)*zdim]-vx[pos-ll*zdim]);
            _vyy+=coeff[ll]*(vy[pos+ll*zdim]-vy[pos-(ll+1)*zdim]);
            _vzy+=coeff[ll]*(vz[pos+(ll+1)*zdim]-vz[pos-ll*zdim]);

            _vxz+=coeff[ll]*(vx[pos+ll+1]-vx[pos-ll]);
            _vyz+=coeff[ll]*(vy[pos+ll+1]-vy[pos-ll]);
            _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-1-ll]);
        }
        _vxx/=dx;_vyx/=dx;_vzx/=dx;
        _vxy/=dy;_vyy/=dy;_vzy/=dy;
        _vxz/=dz;_vyz/=dz;_vzz/=dz;
        const float lambda_plus_2mu=dense[pos]*vp[pos]*vs[pos];
        const float mu=dense[pos]*vs[pos]*vs[pos];
        const float lambda=lambda_plus_2mu-2.0f*mu;
        sigma_xx[pos]+=dt*(lambda_plus_2mu*_vxx+lambda*(_vyy+_vzz));
        sigma_yy[pos]+=dt*(lambda_plus_2mu*_vyy+lambda*(_vxx+_vzz));
        sigma_zz[pos]+=dt*(lambda_plus_2mu*_vzz+lambda*(_vxx+_vyy));
        tau_x[pos]+=dt*mu*(_vyz+_vzy);
        tau_y[pos]+=dt*mu*(_vxz+_vzx);
        tau_z[pos]+=dt*mu*(_vxy+_vyx);
    }
}


void  ElSgStressToVelocity3D(
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
)
{
     const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+xbegin)*ydim+jj+ybegin)*zdim+kk+zbegin;
        float _sigma_x_x=0.0f;
        float _sigma_y_y=0.0f;
        float _sigma_z_z=0.0f;
        float _tau_x_y=0.0f;
        float _tau_x_z=0.0f;
        float _tau_y_x=0.0f;
        float _tau_y_z=0.0f;
        float _tau_z_x=0.0f;
        float _tau_z_y=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _sigma_x_x+=coeff[ll]*(sigma_xx[pos+(ll+1)*yzdim]-sigma_xx[pos-ll*yzdim]);
            _tau_y_x+=coeff[ll]*(tau_y[pos+ll*yzdim]-tau_y[pos-(ll+1)*yzdim]);
            _tau_z_x+=coeff[ll]*(tau_z[pos+ll*yzdim]-tau_z[pos-(ll+1)*yzdim]);

            _sigma_y_y+=coeff[ll]*(sigma_yy[pos+(ll+1)*zdim]-sigma_yy[pos-ll*zdim]);
            _tau_x_y+=coeff[ll]*(tau_x[pos+ll*zdim]-tau_x[pos-(ll+1)*zdim]);
            _tau_z_y+=coeff[ll]*(tau_z[pos+ll*zdim]-tau_z[pos-(ll+1)*zdim]);

            _sigma_z_z+=coeff[ll]*(sigma_zz[pos+ll+1]-sigma_zz[pos-ll]);
            _tau_x_z+=coeff[ll]*(tau_x[pos+ll]-tau_x[pos-1-ll]);
            _tau_y_z+=coeff[ll]*(tau_y[pos+ll]-tau_y[pos-1-ll]);
        }
        _sigma_x_x/=dx; _sigma_y_y/=dy; _sigma_z_z/=dz;
        _tau_z_x/=dx;   _tau_x_y/=dy;   _tau_x_z/=dz;
        _tau_y_x/=dx;   _tau_z_y/=dy;   _tau_y_z/=dz;
        vx[pos]+=dt*2.0f/(dense[pos]+dense[pos+yzdim])*(_sigma_x_x+_tau_y_z+_tau_z_y);
        vy[pos]+=dt*2.0f/(dense[pos]+dense[pos+zdim])*(_sigma_y_y+_tau_x_z+_tau_z_x);
        vz[pos]+=dt*2.0f/(dense[pos]+dense[pos+1])*(_sigma_z_z+_tau_x_y+_tau_y_x);
    }
}



//***************************************************************CPML 吸收边界条件****************************************
void  ElSgVelocityToStressCPML3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+order)*ydim+jj+order)*zdim+kk+order;
        float _vxx=0.0f;
        float _vxy=0.0f;
        float _vxz=0.0f;
        float _vyx=0.0f;
        float _vyy=0.0f;
        float _vyz=0.0f;
        float _vzx=0.0f;
        float _vzy=0.0f;
        float _vzz=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _vxx+=coeff[ll]*(vx[pos+ll*yzdim]-vx[pos-(ll+1)*yzdim]);
            _vyx+=coeff[ll]*(vy[pos+(ll+1)*yzdim]-vy[pos-ll*yzdim]);
            _vzx+=coeff[ll]*(vz[pos+(ll+1)*yzdim]-vz[pos-ll*yzdim]);

            _vxy+=coeff[ll]*(vx[pos+(ll+1)*zdim]-vx[pos-ll*zdim]);
            _vyy+=coeff[ll]*(vy[pos+ll*zdim]-vy[pos-(ll+1)*zdim]);
            _vzy+=coeff[ll]*(vz[pos+(ll+1)*zdim]-vz[pos-ll*zdim]);

            _vxz+=coeff[ll]*(vx[pos+ll+1]-vx[pos-ll]);
            _vyz+=coeff[ll]*(vy[pos+ll+1]-vy[pos-ll]);
            _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-1-ll]);
        }
        _vxx/=dx;_vyx/=dx;_vzx/=dx;
        _vxy/=dy;_vyy/=dy;_vzy/=dy;
        _vxz/=dz;_vyz/=dz;_vzz/=dz;
        const float lambda_plus_2mu=dense[pos]*vp[pos]*vp[pos];
        const float mu=dense[pos]*vs[pos]*vs[pos];
        const float lambda=lambda_plus_2mu-mu-mu;
        sigma_xx[pos]+=dt*(lambda_plus_2mu*_vxx+lambda*(_vyy+_vzz));
        sigma_yy[pos]+=dt*(lambda_plus_2mu*_vyy+lambda*(_vxx+_vzz));
        sigma_zz[pos]+=dt*(lambda_plus_2mu*_vzz+lambda*(_vxx+_vyy));
        tau_x[pos]+=dt*mu*(_vyz+_vzy);
        tau_y[pos]+=dt*mu*(_vxz+_vzx);
        tau_z[pos]+=dt*mu*(_vxy+_vyx);

        //进入吸收边界条件
         if(ii<layers)//进入x方向的吸收边界
        {
            //[layers][ylen][zlen]
            const size_t cpml_pos=(ii*ylen+jj)*zlen+kk;
            sigma_xx[pos]+=dt*lambda_plus_2mu*vx_dx_phi1[cpml_pos];
            sigma_yy[pos]+=dt*lambda*vx_dx_phi1[cpml_pos];
            sigma_zz[pos]+=dt*lambda*vx_dx_phi1[cpml_pos];
            tau_y[pos]+=dt*mu*vz_dx_phi1[cpml_pos];
            tau_z[pos]+=dt*mu*vy_dx_phi1[cpml_pos];
            vx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vxx;
            vy_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vy_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vyx;
            vz_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vz_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vzx; 
        }
        else if(ii>(xlen-layers-1))
        {
            const size_t cpml_pos=((xlen-ii-1)*ylen+jj)*zlen+kk;
            sigma_xx[pos]+=dt*lambda_plus_2mu*vx_dx_phi2[cpml_pos];
            sigma_yy[pos]+=dt*lambda*vx_dx_phi2[cpml_pos];
            sigma_zz[pos]+=dt*lambda*vx_dx_phi2[cpml_pos];
            tau_y[pos]+=dt*mu*vz_dx_phi2[cpml_pos];
            tau_z[pos]+=dt*mu*vy_dx_phi2[cpml_pos];
            vx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen - ii - 1]*_vxx;
            vy_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vy_dx_phi2[cpml_pos]+x2_cpml_beta[xlen - ii - 1]*_vyx;
            vz_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vz_dx_phi2[cpml_pos]+x2_cpml_beta[xlen - ii - 1]*_vzx;
        }

        if(jj<layers)
        {
            const size_t cpml_pos=(jj*xlen+ii)*zlen+kk;
            sigma_yy[pos]+=dt*lambda_plus_2mu*vy_dy_phi1[cpml_pos];
            sigma_xx[pos]+=dt*lambda*vy_dy_phi1[cpml_pos];
            sigma_zz[pos]+=dt*lambda*vz_dy_phi1[cpml_pos];
            tau_x[pos]+=dt*mu*vz_dy_phi1[cpml_pos];
            tau_z[pos]+=dt*mu*vx_dy_phi1[cpml_pos];
            vy_dy_phi1[cpml_pos]=vy_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_vyy;
            vx_dy_phi1[cpml_pos]=vx_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_vxy;
            vz_dy_phi1[cpml_pos]=vz_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_vzy;
           
        }
        else if(jj>(ylen-layers-1))
        {
            const size_t cpml_pos=((ylen-1-jj)*xlen+ii)*zlen+kk;
            sigma_yy[pos]+=dt*lambda_plus_2mu*vy_dy_phi2[cpml_pos];
            sigma_xx[pos]+=dt*lambda*vy_dy_phi2[cpml_pos];
            sigma_zz[pos]+=dt*lambda*vz_dy_phi2[cpml_pos];
            tau_x[pos]+=dt*mu*vz_dy_phi2[cpml_pos];
            tau_z[pos]+=dt*mu*vx_dy_phi2[cpml_pos];
            vy_dy_phi2[cpml_pos]=vy_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_vyy;
            vx_dy_phi2[cpml_pos]=vx_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_vxy;
            vz_dy_phi2[cpml_pos]=vz_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_vzy;
        }

        if(kk<layers)
        {
            const size_t cpml_pos=(kk*xlen+ii)*ylen+jj;
            sigma_zz[pos]+=dt*lambda_plus_2mu*vz_dz_phi1[cpml_pos];
            sigma_xx[pos]+=dt*lambda*vz_dz_phi1[cpml_pos];
            sigma_yy[pos]+=dt*lambda*vz_dz_phi1[cpml_pos];
            tau_x[pos]+=dt*mu*vy_dz_phi1[cpml_pos];
            tau_y[pos]+=dt*mu*vx_dz_phi1[cpml_pos];
            vz_dz_phi1[cpml_pos]=vz_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_vzz;
            vy_dz_phi1[cpml_pos]=vy_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_vyz;
            vx_dz_phi1[cpml_pos]=vx_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_vxz;
            
        }
        else if(kk>(zlen-layers-1))
        {
            const size_t cpml_pos=((zlen-kk-1)*xlen+ii)*ylen+jj;
            sigma_zz[pos]+=dt*lambda_plus_2mu*vz_dz_phi2[cpml_pos];
            sigma_xx[pos]+=dt*lambda*vz_dz_phi2[cpml_pos];
            sigma_yy[pos]+=dt*lambda*vz_dz_phi2[cpml_pos];
            tau_x[pos]+=dt*mu*vy_dz_phi2[cpml_pos];
            tau_y[pos]+=dt*mu*vx_dz_phi2[cpml_pos];
            vz_dz_phi2[cpml_pos]=vz_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_vzz;
            vy_dz_phi2[cpml_pos]=vy_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_vyz;
            vx_dz_phi2[cpml_pos]=vx_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_vxz;
        }
    }
}



void  ElSgStressToVelocityCPML3D(
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
)
{
    const size_t yzdim=ydim*zdim;
    #pragma omp parallel for collapse(3)
    for(int ii=0;ii<xlen;++ii)
    for(int jj=0;jj<ylen;++jj)
    for(int kk=0;kk<zlen;++kk)
    {
        const size_t pos=((ii+order)*ydim+jj+order)*zdim+kk+order;
        float _sigma_x_x=0.0f;
        float _sigma_y_y=0.0f;
        float _sigma_z_z=0.0f;
        float _tau_x_y=0.0f;
        float _tau_x_z=0.0f;
        float _tau_y_x=0.0f;
        float _tau_y_z=0.0f;
        float _tau_z_x=0.0f;
        float _tau_z_y=0.0f;
        for(int ll=0;ll<order;++ll)
        {
            _sigma_x_x+=coeff[ll]*(sigma_xx[pos+(ll+1)*yzdim]-sigma_xx[pos-ll*yzdim]);
            _tau_y_x+=coeff[ll]*(tau_y[pos+ll*yzdim]-tau_y[pos-(ll+1)*yzdim]);
            _tau_z_x+=coeff[ll]*(tau_z[pos+ll*yzdim]-tau_z[pos-(ll+1)*yzdim]);

            _sigma_y_y+=coeff[ll]*(sigma_yy[pos+(ll+1)*zdim]-sigma_yy[pos-ll*zdim]);
            _tau_x_y+=coeff[ll]*(tau_x[pos+ll*zdim]-tau_x[pos-(ll+1)*zdim]);
            _tau_z_y+=coeff[ll]*(tau_z[pos+ll*zdim]-tau_z[pos-(ll+1)*zdim]);

            _sigma_z_z+=coeff[ll]*(sigma_zz[pos+ll+1]-sigma_zz[pos-ll]);
            _tau_x_z+=coeff[ll]*(tau_x[pos+ll]-tau_x[pos-1-ll]);
            _tau_y_z+=coeff[ll]*(tau_y[pos+ll]-tau_y[pos-1-ll]);
        }
        _sigma_x_x/=dx; _sigma_y_y/=dy; _sigma_z_z/=dz;
        _tau_z_x/=dx;   _tau_x_y/=dy;   _tau_x_z/=dz;
        _tau_y_x/=dx;   _tau_z_y/=dy;   _tau_y_z/=dz;
        const float x_dense=dt*2.0f/(dense[pos]+dense[pos+yzdim]);
        const float y_dense=dt*2.0f/(dense[pos]+dense[pos+zdim]);
        const float z_dense=dt*2.0f/(dense[pos]+dense[pos+1]);
        vx[pos]+=x_dense*(_sigma_x_x+_tau_y_z+_tau_z_y);
        vy[pos]+=y_dense*(_sigma_y_y+_tau_x_z+_tau_z_x);
        vz[pos]+=z_dense*(_sigma_z_z+_tau_x_y+_tau_y_x);
        if(ii<layers)
        {
            const size_t cpml_pos=(ii*ylen+jj)*zlen+kk;
            vx[pos]+=x_dense*sigma_xx_dx_phi1[cpml_pos];
            vy[pos]+=y_dense*tau_z_dx_phi1[cpml_pos];
            vz[pos]+=z_dense*tau_y_dx_phi1[cpml_pos];
            sigma_xx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*sigma_xx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_sigma_x_x;
            tau_y_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*tau_y_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_tau_y_x;
            tau_z_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*tau_z_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_tau_z_x;
           
        }
        else if(ii>(xlen-layers-1))
        {
            const size_t cpml_pos=((xlen-ii-1)*ylen+jj)*zlen+kk;
            vx[pos]+=x_dense*sigma_xx_dx_phi2[cpml_pos];
            vy[pos]+=y_dense*tau_z_dx_phi2[cpml_pos];
            vz[pos]+=z_dense*tau_y_dx_phi2[cpml_pos];
            sigma_xx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*sigma_xx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_sigma_x_x;
            tau_y_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*tau_y_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_tau_y_x;
            tau_z_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*tau_z_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_tau_z_x;
           
        }

        if(jj<layers)
        {
            const size_t cpml_pos=(jj*xlen+ii)*zlen+kk;
            vx[pos]+=x_dense*tau_z_dy_phi1[cpml_pos];
            vy[pos]+=y_dense*sigma_yy_dy_phi1[cpml_pos];
            vz[pos]+=z_dense*tau_x_dy_phi1[cpml_pos];
            sigma_yy_dy_phi1[cpml_pos]=sigma_yy_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_sigma_y_y;
            tau_x_dy_phi1[cpml_pos]=tau_x_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_tau_x_y;
            tau_z_dy_phi1[cpml_pos]=tau_z_dy_phi1[cpml_pos]*y1_cpml_alpha[jj]+y1_cpml_beta[jj]*_tau_z_y;
           
        }
        else if(jj>(ylen-layers-1))
        {
            const size_t cpml_pos=((ylen-1-jj)*xlen+ii)*zlen+kk;
            vx[pos]+=x_dense*tau_z_dy_phi2[cpml_pos];
            vy[pos]+=y_dense*sigma_yy_dy_phi2[cpml_pos];
            vz[pos]+=z_dense*tau_x_dy_phi2[cpml_pos];
            sigma_yy_dy_phi2[cpml_pos]=sigma_yy_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_sigma_y_y;
            tau_x_dy_phi2[cpml_pos]=tau_x_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_tau_x_y;
            tau_z_dy_phi2[cpml_pos]=tau_z_dy_phi2[cpml_pos]*y2_cpml_alpha[ylen-jj-1]+y2_cpml_beta[ylen-jj-1]*_tau_z_y;
           
        }

        if(kk<layers)
        {
             const size_t cpml_pos=(kk*xlen+ii)*ylen+jj;
             vx[pos]+=x_dense*tau_y_dz_phi1[cpml_pos];
             vy[pos]+=y_dense*tau_x_dz_phi1[cpml_pos];
             vz[pos]+=z_dense*sigma_zz_dz_phi1[cpml_pos];
             sigma_zz_dz_phi1[cpml_pos]=sigma_zz_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_sigma_z_z;
             tau_x_dz_phi1[cpml_pos]=tau_x_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_tau_x_z;
             tau_y_dz_phi1[cpml_pos]=tau_y_dz_phi1[cpml_pos]*z1_cpml_alpha[kk]+z1_cpml_beta[kk]*_tau_y_z;
            
        }
        else if(kk>(zlen-layers-1))
        {
            const size_t cpml_pos=((zlen-kk-1)*xlen+ii)*ylen+jj;
            vx[pos]+=x_dense*tau_y_dz_phi2[cpml_pos];
            vy[pos]+=y_dense*tau_x_dz_phi2[cpml_pos];
            vz[pos]+=z_dense*sigma_zz_dz_phi2[cpml_pos];
            sigma_zz_dz_phi2[cpml_pos]=sigma_zz_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_sigma_z_z;
            tau_x_dz_phi2[cpml_pos]=tau_x_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_tau_x_z;
            tau_y_dz_phi2[cpml_pos]=tau_y_dz_phi2[cpml_pos]*z2_cpml_alpha[zlen-kk-1]+z2_cpml_beta[zlen-kk-1]*_tau_y_z;
            
        }
    }
}