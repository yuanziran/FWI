#pragma once
#include"fwi_cpu_sys_info.h"
#include<string>
#include"model_cmd.hpp"
namespace FWI
{
    /**
     * 模型文件的格式定义
     * size_t   dim1
     * size_t   dim2
     * float*   data
    */
    class DLL_API model
    {
        public:
        static void ReadModelData2D(const std::string& path,float* data,const size_t& dim1,const size_t& dim2);
		/**
		*使用New 申请内存
		* */
		static float* ReadModelData2DEx(const std::string& path,size_t& dim1,size_t& dim2);
        static void WriteModelData2D(const std::string& path,const float* data,const size_t& dim1,const size_t& dim2);
        static void CopyToArray(const float value,float* data,const size_t& len);
        static float AcVp2K(const float dense,const float vp);
        static float AcK2Vp(const float dense,const float k);
        static void VpVs2LamudaMu(const float dense,const float vp,const float vs,float &lamuda,float &mu);
        static void LamudaMu2VpVs(const float dense,const float lamuda,const float mu,float &vp,float &vs);
        
        protected:
        model_cmd   m_cmd;
        public:
        model(int argc,char**argv);
        void debug()const;
        ~model();
        public:
        void loadAcVelocityModel(const float dense,const float vp);                                         ///<加载声波速度模型
        void loadAcStressVelocityModel(const float dense,const float k);                                    ///<加载声波应力速度模型
        void loadAcVelocityModel(const char* dense,const char* vp);
        void loadAcStressVelocityModel(const char* dense,const char* k);
        void loadElVelocityModel(const float dense,const float vp,const float vs);                          ///<加载弹性波速度模型
        void loadElStressVelocityModel(const float dense,const float lambda,const float mu);                ///<加载弹性波应力速度模型
        void loadElVelocityModel(const char* dense,const char* vp,const char* vs);
        void loadElStressVelocityModel(const char* dense,const char* lambda,const char* mu);
        void setAcVelocityModel2D(float* const dense,float* const vp,const size_t& dim1,const size_t& dim2);                            
        void setAcStressVelocityModel2D(float* const dense,float* const k,const size_t& dim1,const size_t& dim2);
        void setElVelocityModel2D(float* const dense,float* const vp,float* const vs,const size_t& dim1,const size_t& dim2);
        void setElStressVelocityModel2D(float* const dense,float* const lambda,float* const mu,const size_t& dim1,const size_t& dim2);

    };
}