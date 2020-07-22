#include"model.hpp"
#include<fstream>
#include<iostream>
using namespace FWI;

void model:: ReadModelData2D(const std::string & path,float* data,const size_t& dim1,const size_t& dim2)
{
    std::ifstream in(path.c_str(),std::ios_base::binary);
    assert(in.is_open());
    size_t xdim=0;
    size_t ydim=0;
    in.read((char*)&xdim,sizeof(xdim));
    in.read((char*)&ydim,sizeof(ydim));
    assert((xdim==dim1)&&(ydim==dim2));
    in.read((char*)data,sizeof(float)*xdim*ydim);
}

float* model:: ReadModelData2DEx(const std::string& path, size_t& dim1, size_t& dim2)
{
	std::ifstream in(path.c_str(), std::ios_base::binary);
	assert(in.is_open());
	in.read((char*)&dim1,sizeof(size_t));
	in.read((char*)&dim2,sizeof(size_t));
	float* data = new float[dim1*dim2];
	in.read((char*)data,sizeof(float)*dim1*dim2);
	return data;
}

void model:: WriteModelData2D(const std::string& path,const float* data,const size_t& dim1,const size_t& dim2)
{
    std::ofstream out(path.c_str(),std::ios_base::binary|std::ios_base::trunc);
    assert(out.is_open());
    out.write((const char*)&dim1,sizeof(size_t));
    out.write((const char*)&dim2,sizeof(size_t));
    out.write((const char*)data,sizeof(float)*dim1*dim2);
}

void model::CopyToArray(const float value,float* data,const size_t& len)
{
    #pragma omp parallel for
    for(size_t i=0;i<len;++i)
    {
        data[i]=value;
    }
}


float model::AcVp2K(const float dense,const float vp)
{
     return dense*vp*vp;
}

float model::AcK2Vp(const float dense,const float k)
{
    return sqrtf(k/dense);
}

void model::VpVs2LamudaMu(const float dense,const float vp,const float vs,float &lamuda,float &mu)
{
    mu=vs*vs*dense;
    lamuda=vp*vp*dense+2*mu;
}

void model:: LamudaMu2VpVs(const float dense,const float lambda,const float mu,float &vp,float &vs)
{
    vp=sqrtf((lambda+2*mu)/dense);
    vs=sqrtf(mu/dense);
}

 model::model(int argc,char**argv)
 {
     ParseModelCommandLine(argc,argv,m_cmd);
 }

void model::debug()const
{
    DebugModelCommandLine(m_cmd);
}

model::~model()
{
    if(m_cmd.type==MODEL_TYPE::AC_MODEL)
    {
        if(m_cmd.ac_model_param.dense) DestoryCString(m_cmd.ac_model_param.dense);
        if(m_cmd.ac_model_param.k_or_vp) DestoryCString(m_cmd.ac_model_param.k_or_vp);
    }
    else if(m_cmd.type==MODEL_TYPE::EL_MODEL)
    {
        if(m_cmd.el_model_param.dense) DestoryCString(m_cmd.el_model_param.dense);
        if(m_cmd.el_model_param.vp_or_lambda) DestoryCString(m_cmd.el_model_param.vp_or_lambda);
        if(m_cmd.el_model_param.vs_or_mu) DestoryCString(m_cmd.el_model_param.vs_or_mu);
    }
}


void model::loadAcVelocityModel(const float dense,const float vp)
{
    m_cmd.type=MODEL_TYPE::AC_FIXED_VALUE_MODEL;
    m_cmd.flag.open_velocity_model_flag=1;
    m_cmd.flag.open_valid_dense_flag=1;
    m_cmd.ac_model_fixed_value_param.dense=dense;
    m_cmd.ac_model_fixed_value_param.k_or_vp=vp;
}

void model::loadAcStressVelocityModel(const float dense,const float k)
{
    m_cmd.type=MODEL_TYPE::AC_FIXED_VALUE_MODEL;
    m_cmd.flag.open_valid_dense_flag=1;
    m_cmd.flag.open_velocity_model_flag=0;
    m_cmd.ac_model_fixed_value_param.dense=dense;
    m_cmd.ac_model_fixed_value_param.k_or_vp=k;
}

void model::loadAcVelocityModel(const char* dense,const char* vp)
{
    m_cmd.type=MODEL_TYPE::AC_MODEL;
    m_cmd.flag.open_velocity_model_flag=1;
    if(dense)
    {
        m_cmd.flag.open_valid_dense_flag=1;
        if(m_cmd.ac_model_param.dense) DestoryCString(m_cmd.ac_model_param.dense);
        m_cmd.ac_model_param.dense=AllocateAndInitCString(dense);
    }
    if(m_cmd.ac_model_param.k_or_vp) DestoryCString(m_cmd.ac_model_param.k_or_vp);
    m_cmd.ac_model_param.k_or_vp=AllocateAndInitCString(vp);
}
        
void model::loadAcStressVelocityModel(const char* dense,const char* k)
{
    m_cmd.type=MODEL_TYPE::AC_MODEL;
    m_cmd.flag.open_velocity_model_flag=0;
    if(dense)
    {
        m_cmd.flag.open_valid_dense_flag=1;
        if(m_cmd.ac_model_param.dense) DestoryCString(m_cmd.ac_model_param.dense);
        m_cmd.ac_model_param.dense=AllocateAndInitCString(dense);
    }
    if(m_cmd.ac_model_param.k_or_vp) DestoryCString(m_cmd.ac_model_param.k_or_vp);
    m_cmd.ac_model_param.k_or_vp=AllocateAndInitCString(k);
}

void model::loadElVelocityModel(const float dense,const float vp,const float vs)
{
    m_cmd.type=MODEL_TYPE::EL_FIXED_VALUE_MODEL;
    m_cmd.flag.open_velocity_model_flag=1;
    m_cmd.flag.open_valid_dense_flag=1;
    m_cmd.el_model_fixed_value_param.dense=dense;
    m_cmd.el_model_fixed_value_param.vp_or_lambda=vp;
    m_cmd.el_model_fixed_value_param.vs_or_mu=vs;
}

void model::loadElStressVelocityModel(const float dense,const float lambda,const float mu)
{
    m_cmd.type=MODEL_TYPE::EL_FIXED_VALUE_MODEL;
    m_cmd.flag.open_velocity_model_flag=0;
    m_cmd.flag.open_valid_dense_flag=1;
    m_cmd.el_model_fixed_value_param.dense=dense;
    m_cmd.el_model_fixed_value_param.vp_or_lambda=lambda;
    m_cmd.el_model_fixed_value_param.vs_or_mu=mu;
}

 void model::loadElVelocityModel(const char* dense,const char* vp,const char* vs)
 {
     m_cmd.type=MODEL_TYPE::EL_MODEL;
     m_cmd.flag.open_velocity_model_flag=1;
     if(dense)
     {
         m_cmd.flag.open_valid_dense_flag=1;
         if(m_cmd.el_model_param.dense) DestoryCString(m_cmd.el_model_param.dense);
         m_cmd.el_model_param.dense=AllocateAndInitCString(dense);
     }
     if(m_cmd.el_model_param.vp_or_lambda) DestoryCString(m_cmd.el_model_param.vp_or_lambda);
     m_cmd.el_model_param.vp_or_lambda=AllocateAndInitCString(vp);
     if(m_cmd.el_model_param.vs_or_mu) DestoryCString(m_cmd.el_model_param.vs_or_mu);
     m_cmd.el_model_param.vs_or_mu=AllocateAndInitCString(vs);
 }

void model:: loadElStressVelocityModel(const char* dense,const char* lambda,const char* mu)
{
    m_cmd.type=MODEL_TYPE::EL_MODEL;
    m_cmd.flag.open_velocity_model_flag=0;
    if(dense)
    {
         m_cmd.flag.open_valid_dense_flag=1;
         if(m_cmd.el_model_param.dense) DestoryCString(m_cmd.el_model_param.dense);
         m_cmd.el_model_param.dense=AllocateAndInitCString(dense);
    }
    if(m_cmd.el_model_param.vp_or_lambda) DestoryCString(m_cmd.el_model_param.vp_or_lambda);
    m_cmd.el_model_param.vp_or_lambda=AllocateAndInitCString(lambda);
    if(m_cmd.el_model_param.vs_or_mu) DestoryCString(m_cmd.el_model_param.vs_or_mu);
    m_cmd.el_model_param.vs_or_mu=AllocateAndInitCString(mu);
}

void model:: setAcVelocityModel2D(float* const dense,float* const vp,const size_t& dim1,const size_t& dim2)
{
    const size_t buffer_len=dim1*dim2;
    if(m_cmd.type==MODEL_TYPE::AC_FIXED_VALUE_MODEL)
    {
       CopyToArray(m_cmd.ac_model_fixed_value_param.dense,dense,buffer_len);
       if(m_cmd.flag.open_velocity_model_flag)
       {
           CopyToArray(m_cmd.ac_model_fixed_value_param.k_or_vp,vp,buffer_len);
       }
       else
       {
           float _vp=AcK2Vp(m_cmd.ac_model_fixed_value_param.dense,m_cmd.ac_model_fixed_value_param.k_or_vp);
           CopyToArray(_vp,vp,buffer_len);
       }
    }
    else if(m_cmd.type==MODEL_TYPE::AC_MODEL && m_cmd.flag.open_velocity_model_flag==1)
    {
        if(m_cmd.flag.open_velocity_model_flag)
        {
            ReadModelData2D(m_cmd.ac_model_param.dense,dense,dim1,dim2);
        }
        ReadModelData2D(m_cmd.ac_model_param.k_or_vp,vp,dim1,dim2);
    }
    else
    {
        std::cerr<<"bad type or flag"<<std::endl;
        exit(EXIT_FAILURE);
    }
}  

void model:: setAcStressVelocityModel2D(float* const dense,float* const k,const size_t& dim1,const size_t& dim2)
{
    const size_t buffer_len=dim1*dim2;
    if(m_cmd.type==MODEL_TYPE::AC_FIXED_VALUE_MODEL)
    {
       CopyToArray(m_cmd.ac_model_fixed_value_param.dense,dense,buffer_len);
       if(m_cmd.flag.open_velocity_model_flag)
       {
           float _k=AcVp2K(m_cmd.ac_model_fixed_value_param.dense,m_cmd.ac_model_fixed_value_param.k_or_vp);
           CopyToArray(_k,k,buffer_len);
       }
       else
       {
           CopyToArray(m_cmd.ac_model_fixed_value_param.k_or_vp,k,buffer_len);
       }
    }
    else if(m_cmd.type==MODEL_TYPE::AC_MODEL && m_cmd.flag.open_velocity_model_flag==0)
    {
        if(m_cmd.flag.open_velocity_model_flag)
        {
            ReadModelData2D(m_cmd.ac_model_param.dense,dense,dim1,dim2);
        }
        else
        {
            CopyToArray(1.0f,dense,buffer_len);
        }
        
        ReadModelData2D(m_cmd.ac_model_param.k_or_vp,k,dim1,dim2);
    }
    else
    {
        std::cerr<<"bad type or flag"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

void model::setElVelocityModel2D(float* const dense,float* const vp,float* const vs,const size_t& dim1,const size_t& dim2)
{
    const size_t buffer_len=dim1*dim2;
    if(m_cmd.type==MODEL_TYPE::EL_FIXED_VALUE_MODEL)
    {
        CopyToArray(m_cmd.el_model_fixed_value_param.dense,dense,buffer_len);
        if(m_cmd.flag.open_velocity_model_flag==1)
        {
            CopyToArray(m_cmd.el_model_fixed_value_param.vp_or_lambda,vp,buffer_len);
            CopyToArray(m_cmd.el_model_fixed_value_param.vs_or_mu,vs,buffer_len);
        }
        else
        {
            float _vp,_vs;
            LamudaMu2VpVs(m_cmd.el_model_fixed_value_param.dense,m_cmd.el_model_fixed_value_param.vp_or_lambda,
            m_cmd.el_model_fixed_value_param.vs_or_mu,_vp,_vs);
            CopyToArray(_vp,vp,buffer_len);
            CopyToArray(_vs,vs,buffer_len);
        }
    }
    else if(m_cmd.type==MODEL_TYPE::EL_MODEL && m_cmd.flag.open_velocity_model_flag==1)
    {
        if(m_cmd.flag.open_valid_dense_flag==1)
        {
            ReadModelData2D(m_cmd.el_model_param.dense,dense,dim1,dim2);
        }
        else
        {
            CopyToArray(1.0f,dense,buffer_len);
        }
        ReadModelData2D(m_cmd.el_model_param.vp_or_lambda,vp,dim1,dim2);
        ReadModelData2D(m_cmd.el_model_param.vs_or_mu,vs,dim1,dim2);
    }
    else
    {
        std::cerr<<"bad type or flag"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

 void model::setElStressVelocityModel2D(float* const dense,float* const lambda,float* const mu,const size_t& dim1,const size_t& dim2)
 {
     const size_t buffer_len=dim1*dim2;
    if(m_cmd.type==MODEL_TYPE::EL_FIXED_VALUE_MODEL)
    {
        CopyToArray(m_cmd.el_model_fixed_value_param.dense,dense,buffer_len);
        if(m_cmd.flag.open_velocity_model_flag==1)
        {
            float _lambda,_mu;
            VpVs2LamudaMu(m_cmd.el_model_fixed_value_param.dense,m_cmd.el_model_fixed_value_param.vp_or_lambda,
            m_cmd.el_model_fixed_value_param.vs_or_mu,_lambda,_mu);
            CopyToArray(_lambda,lambda,buffer_len);
            CopyToArray(_mu,mu,buffer_len);
        }
        else
        {
            CopyToArray(m_cmd.el_model_fixed_value_param.vp_or_lambda,lambda,buffer_len);
            CopyToArray(m_cmd.el_model_fixed_value_param.vs_or_mu,mu,buffer_len);
        }
    }
    else if(m_cmd.type==MODEL_TYPE::EL_MODEL && m_cmd.flag.open_velocity_model_flag==0)
    {
        if(m_cmd.flag.open_valid_dense_flag==1)
        {
            ReadModelData2D(m_cmd.el_model_param.dense,dense,dim1,dim2);
        }
        else
        {
            CopyToArray(1.0f,dense,buffer_len);
        }
        ReadModelData2D(m_cmd.el_model_param.vp_or_lambda,lambda,dim1,dim2);
        ReadModelData2D(m_cmd.el_model_param.vs_or_mu,mu,dim1,dim2);
    }
    else
    {
        std::cerr<<"bad type or flag"<<std::endl;
        exit(EXIT_FAILURE);
    }
 }