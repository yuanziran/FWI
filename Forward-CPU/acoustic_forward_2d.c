#include"acoustic_forward_2d.h"

void  AcSgVelocityToStress2D(
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
)
{
      #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+xbegin)*zdim+(jj+zbegin);
          float _vxx=0.0f;
          float _vzz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _vxx+=coeff[ll]*(vx[pos+ll*zdim]-vx[pos-(ll+1)*zdim]);
              _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
          }
          _vxx/=dx;_vzz/=dz;
          stress[pos]+=dt*dense[pos]*vp[pos]*vp[pos]*(_vxx+_vzz);
      }
}


void  AcSgStressToVelocity2D(
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
)
{
      #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+xbegin)*zdim+(jj+zbegin);
          float _px=0.0f;
          float _pz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _px+=coeff[ll]*(stress[pos+(ll+1)*zdim]-stress[pos-ll*zdim]);
              _pz+=coeff[ll]*(stress[pos+ll+1]-stress[pos-ll]);
          }
          _px/=dx;_pz/=dz;
          vx[pos]+=2.0f*dt/(dense[pos+zdim]+dense[pos])*_px;
          vz[pos]+=2.0f*dt/(dense[pos+1]+dense[pos])*_pz;
      }
}

//**************************************************卷积一致吸收边界条件********************************
void  AcSgVelocityToStressCPML2D(
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
)
{
      #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _vxx=0.0f;
          float _vzz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _vxx+=coeff[ll]*(vx[pos+ll*zdim]-vx[pos-(ll+1)*zdim]);
              _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
          }
          _vxx/=dx;_vzz/=dz;
          const float k=dt*dense[pos]*vp[pos]*vp[pos];
          stress[pos]+=k*(_vxx+_vzz);
          //进入吸收边界
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              stress[pos]+=k*vx_dx_phi1[cpml_pos];
              vx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vxx;
          }
          else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              stress[pos]+=k*vx_dx_phi2[cpml_pos];
              vx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_vxx;
          }
          //进入吸收边界
          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              stress[pos]+=k*vz_dz_phi1[cpml_pos];
              vz_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*vz_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_vzz;
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              stress[pos]+=k*vz_dz_phi2[cpml_pos];
              vz_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*vz_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_vzz;
          }
      }
}

void  AcSgStressToVelocityCPML2D(
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
)
{
      #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _px=0.0f;
          float _pz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _px+=coeff[ll]*(stress[pos+(ll+1)*zdim]-stress[pos-ll*zdim]);
              _pz+=coeff[ll]*(stress[pos+ll+1]-stress[pos-ll]);
          }
          const float x_scalar=2.0f*dt/(dense[pos+zdim]+dense[pos]);
          const float z_scalar=2.0f*dt/(dense[pos+1]+dense[pos]);
          vx[pos]+=x_scalar*_px;
          vz[pos]+=z_scalar*_pz;
          //进入吸收边界条件
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              vx[pos]+=x_scalar*stress_dx_phi1[cpml_pos];
              stress_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*stress_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_px;
          }
          else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              vx[pos]+=x_scalar*stress_dx_phi2[cpml_pos];
              stress_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*stress_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_px;
          }
          //进入吸收边界
          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              vz[pos]+=z_scalar*stress_dz_phi1[cpml_pos];
              stress_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*stress_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_pz;
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              vz[pos]+=z_scalar*stress_dz_phi2[cpml_pos];
              stress_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*stress_dz_phi2[cpml_pos]+z1_cpml_beta[zlen-jj-1]*_pz;
          }
      }
}



void  AcSgVelocityToStressCPML2DEx(
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
)
{
     #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _vxx=0.0f;
          float _vzz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _vxx+=coeff[ll]*(vx[pos+ll*zdim]-vx[pos-(ll+1)*zdim]);
              _vzz+=coeff[ll]*(vz[pos+ll]-vz[pos-ll-1]);
          }
          _vxx/=dx;_vzz/=dz;
          const float k=dt*dense[pos]*vp[pos]*vp[pos];
          stress[pos]+=k*(_vxx+_vzz);
          //进入吸收边界
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              stress[pos]+=k*(vx_dx_phi1[cpml_pos]+(x1_cpml_kappa[ii]-1.0f)*_vxx);
              //vx_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*vx_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_vxx*x1_cpml_kappa[ii];
			  vx_dx_phi1[cpml_pos] = x1_cpml_alpha[ii] * vx_dx_phi1[cpml_pos] + x1_cpml_beta[ii] * _vxx;
          }
          else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              //stress[pos]+=k*(vx_dx_phi2[cpml_pos]+(x2_cpml_kappa[xlen-1-ii]-1.0f)*_vxx);
			  stress[pos] += k * (vx_dx_phi2[cpml_pos]);
              vx_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*vx_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_vxx*x2_cpml_kappa[xlen-1-ii];
			  //vx_dx_phi2[cpml_pos] = x2_cpml_alpha[xlen - ii - 1] * vx_dx_phi2[cpml_pos] + x2_cpml_beta[xlen - ii - 1] * _vxx;
          }
          //进入吸收边界
          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              stress[pos]+=k*(vz_dz_phi1[cpml_pos]+(z1_cpml_kappa[jj]-1.0f)*_vzz);
              vz_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*vz_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_vzz*z1_cpml_kappa[jj];
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              stress[pos]+=k*(vz_dz_phi2[cpml_pos]+(z2_cpml_kappa[zlen-1-jj]-1.0f)*_vzz);
              vz_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*vz_dz_phi2[cpml_pos]+z2_cpml_beta[zlen-jj-1]*_vzz*z2_cpml_kappa[zlen-1-jj];
          }
      }
}


void  AcSgStressToVelocityCPML2DEx(
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
)
{
      #pragma omp parallel for collapse(2)
      for(int ii=0;ii<xlen;++ii)
      for(int jj=0;jj<zlen;++jj)
      {
          const size_t pos=(ii+order)*zdim+(jj+order);
          float _px=0.0f;
          float _pz=0.0f;
          for(int ll=0;ll<order;++ll)
          {
              _px+=coeff[ll]*(stress[pos+(ll+1)*zdim]-stress[pos-ll*zdim]);
              _pz+=coeff[ll]*(stress[pos+ll+1]-stress[pos-ll]);
          }
          const float x_scalar=2.0f*dt/(dense[pos+zdim]+dense[pos]);
          const float z_scalar=2.0f*dt/(dense[pos+1]+dense[pos]);
          vx[pos]+=x_scalar*_px;
          vz[pos]+=z_scalar*_pz;
          //进入吸收边界条件
          if(ii<layers)
          {
              const size_t cpml_pos=ii*zlen+jj;
              vx[pos]+=x_scalar*(stress_dx_phi1[cpml_pos]+(x1_cpml_kappa[ii]-1.0f)*_px);
              //stress_dx_phi1[cpml_pos]=x1_cpml_alpha[ii]*stress_dx_phi1[cpml_pos]+x1_cpml_beta[ii]*_px*x1_cpml_kappa[ii];
			  stress_dx_phi1[cpml_pos] = x1_cpml_alpha[ii] * stress_dx_phi1[cpml_pos] + x1_cpml_beta[ii] * _px;
          }
          else if((ii>(xlen-layers-1)))
          {
              const size_t cpml_pos=(xlen-1-ii)*zlen+jj;
              vx[pos]+=x_scalar*(stress_dx_phi2[cpml_pos]);
              stress_dx_phi2[cpml_pos]=x2_cpml_alpha[xlen-ii-1]*stress_dx_phi2[cpml_pos]+x2_cpml_beta[xlen-ii-1]*_px*x2_cpml_kappa[xlen-ii-1];
			  //stress_dx_phi2[cpml_pos] = x2_cpml_alpha[xlen - ii - 1] * stress_dx_phi2[cpml_pos] + x2_cpml_beta[xlen - ii - 1] * _px;
          }
          //进入吸收边界
          if(jj<layers)
          {
              const size_t cpml_pos=jj*xlen+ii;
              vz[pos]+=z_scalar*(stress_dz_phi1[cpml_pos]+(z1_cpml_kappa[jj]-1.0f)*_pz);
              stress_dz_phi1[cpml_pos]=z1_cpml_alpha[jj]*stress_dz_phi1[cpml_pos]+z1_cpml_beta[jj]*_pz*z1_cpml_kappa[jj];
          }
          else if((jj>(zlen-layers-1)))
          {
              const size_t cpml_pos=(zlen-1-jj)*xlen+ii;
              vz[pos]+=z_scalar*(stress_dz_phi2[cpml_pos]+(z2_cpml_kappa[zlen-jj-1]-1.0f)*_pz);
              stress_dz_phi2[cpml_pos]=z2_cpml_alpha[zlen-jj-1]*stress_dz_phi2[cpml_pos]+z1_cpml_beta[zlen-jj-1]*_pz*z2_cpml_kappa[zlen-jj-1];
          }
      }
}