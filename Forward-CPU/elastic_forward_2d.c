#include"elastic_forward_2d.h"

void  ElSgVelocityToStress2D(
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
)
{
    #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+xbegin)*zdim+(jj+zbegin);
          float _vxx=0.0f;
          float _vzz=0.0f;
          float _vxz=0.0f;
          float _vzx=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _vxx+=coeff[ll]*(vx[pos+ll*zdim]-vx[pos-(ll+1)*zdim]);
              _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
              _vxz+=coeff[ll]*(vx[pos+ll+1]-vx[pos-ll]);
              _vzx+=coeff[ll]*(vz[pos+(ll+1)*zdim]-vz[pos-ll*zdim]);
          }
          _vxx/=dx;_vzz/=dz;_vzx/=dx;_vxz/=dz;
          const float lambda_plus_2mu=dense[pos]*vp[pos]*vs[pos];
          const float mu=dense[pos]*vs[pos]*vs[pos];
          const float lambda=lambda_plus_2mu-mu-mu;
          sigma_xx[pos]+=dt*(lambda_plus_2mu*_vxx+lambda*_vzz);
          sigma_zz[pos]+=dt*(lambda_plus_2mu*_vzz+lambda*_vxx);
          tau_y[pos]+=dt*mu*(_vxz+_vzx);
      }
}


void  ElSgStressToVelocity2D(
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
)
{
     #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+xbegin)*zdim+(jj+zbegin);
          float _px=0.0f;
          float _pz=0.0f;
          float _tx=0.0f;
          float _tz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _px+=coeff[ll]*(sigma_xx[pos+(ll+1)*zdim]-sigma_xx[pos-ll*zdim]);
              _pz+=coeff[ll]*(sigma_zz[pos+ll+1]-sigma_zz[pos-ll]);
              _tx+=coeff[ll]*(tau_y[pos+ll*zdim]-tau_y[pos-(ll+1)*zdim]);
              _tz+=coeff[ll]*(tau_y[pos+ll]-tau_y[pos-ll-1]);
          }
          _px/=dx;_pz/=dz;_tx/=dx;_tz/=dz;
          vx[pos]+=2.0f*dt/(dense[pos+zdim]+dense[pos])*(_px+_tz);
          vz[pos]+=2.0f*dt/(dense[pos+1]+dense[pos])*(_pz+_tx);
      }
}



void  ElSgVelocityToStressCPML2D(
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
)
{
     #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _vxx=0.0f;
          float _vzz=0.0f;
          float _vxz=0.0f;
          float _vzx=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _vxx+=coeff[ll]*(vx[pos+ll*zdim]-vx[pos-(ll+1)*zdim]);
              _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
              _vxz+=coeff[ll]*(vx[pos+ll+1]-vx[pos-ll]);
              _vzx+=coeff[ll]*(vz[pos+(ll+1)*zdim]-vz[pos-ll*zdim]);
          }
          _vxx/=dx;_vzz/=dz;_vzx/=dx;_vxz/=dz;
          const float lambda_plus_2mu=dense[pos]*vp[pos]*vs[pos];
          const float mu=dense[pos]*vs[pos]*vs[pos];
          const float lambda=lambda_plus_2mu-mu-mu;
          sigma_xx[pos]+=dt*(lambda_plus_2mu*_vxx+lambda*_vzz);
          sigma_zz[pos]+=dt*(lambda_plus_2mu*_vzz+lambda*_vxx);
          tau_y[pos]+=dt*mu*(_vxz+_vzx);
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              sigma_xx[pos]+=dt*lambda_plus_2mu*vx_dx_phi1[cpml_pos];
              sigma_zz[pos]+=dt*lambda*vx_dx_phi1[cpml_pos];
              tau_y[pos]+=dt*mu*vz_dx_phi1[cpml_pos];
              vx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vxx;
              vz_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vz_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vzx;
          }
           else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              sigma_xx[pos]+=dt*lambda_plus_2mu*vx_dx_phi2[cpml_pos];
              sigma_zz[pos]+=dt*lambda*vx_dx_phi2[cpml_pos];
              tau_y[pos]+=dt*mu*vz_dx_phi2[cpml_pos];
              vx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_vxx;
              vz_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vz_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_vzx;
          }

          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              sigma_xx[pos]+=dt*lambda*vz_dz_phi1[cpml_pos];
              sigma_zz[pos]+=dt*lambda_plus_2mu*vz_dz_phi1[cpml_pos];
              tau_y[pos]+=dt*mu*vx_dz_phi1[cpml_pos];
              vz_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*vz_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_vzz;
              vx_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*vx_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_vxz;
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              sigma_xx[pos]+=dt*lambda*vz_dz_phi2[cpml_pos];
              sigma_zz[pos]+=dt*lambda_plus_2mu*vz_dz_phi2[cpml_pos];
              tau_y[pos]+=dt*mu*vx_dz_phi2[cpml_pos];
              vz_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*vz_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_vzz;
              vx_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*vx_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_vxz;
          }
      }
}

void  ElSgStressToVelocityCPML2D(
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
)
{
     #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _px=0.0f;
          float _pz=0.0f;
          float _tx=0.0f;
          float _tz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _px+=coeff[ll]*(sigma_xx[pos+(ll+1)*zdim]-sigma_xx[pos-ll*zdim]);
              _pz+=coeff[ll]*(sigma_zz[pos+ll+1]-sigma_zz[pos-ll]);
              _tx+=coeff[ll]*(tau_y[pos+ll*zdim]-tau_y[pos-(ll+1)*zdim]);
              _tz+=coeff[ll]*(tau_y[pos+ll]-tau_y[pos-ll-1]);
          }
          _px/=dx;_pz/=dz;_tx/=dx;_tz/=dz;
          const float x_scalar=2.0f*dt/(dense[pos+zdim]+dense[pos]);
          const float z_scalar=2.0f*dt/(dense[pos+1]+dense[pos]);
          vx[pos]+=x_scalar*(_px+_tz);
          vz[pos]+=z_scalar*(_pz+_tx);
          //进入吸收边界
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              vx[pos]+=x_scalar*(sigma_xx_dx_phi1[cpml_pos]);
              vz[pos]+=z_scalar*(tau_y_dx_phi1[cpml_pos]);
              sigma_xx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*sigma_xx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_px;
              tau_y_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*tau_y_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_tx;
          }
          else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              vx[pos]+=x_scalar*(sigma_xx_dx_phi2[cpml_pos]);
              vz[pos]+=z_scalar*(tau_y_dx_phi2[cpml_pos]);
              sigma_xx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*sigma_xx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_px;
              tau_y_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*tau_y_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_tx;
          }
          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              vx[pos]+=x_scalar*(tau_y_dz_phi1[cpml_pos]);
              vz[pos]+=z_scalar*(sigma_zz_dz_phi1[cpml_pos]);
              tau_y_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*tau_y_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_tz;
              sigma_zz_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*sigma_zz_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_pz;
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              vx[pos]+=x_scalar*(tau_y_dz_phi2[cpml_pos]);
              vz[pos]+=z_scalar*(sigma_zz_dz_phi2[cpml_pos]);
              tau_y_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*tau_y_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_tz;
              sigma_zz_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*sigma_zz_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_pz;
          }
      }
}
