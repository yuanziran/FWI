#pragma once
#include"fwi_cpu_sys_info.h"
#include"mdoel_shape.hpp"
#include<map>
#include<string>
#include<fstream>
namespace FWI
{
    class  model_script
    {
        float*          m_data;                         ///<模型数据
        size_t          m_xdim;                         ///<维度大小
        size_t          m_ydim;                         ///<维度大小
        std::string     out_file_path;                  ///<输出模型文件
        float           m_background_value;             ///<背景值
        public:
        model_script(const char* path);
        ~model_script();
        unsigned int readScriptHeader(std::ifstream &in);
        void debug()const;
        void run(std::ifstream& in,unsigned int &lineNo);
        protected:
        std::map<std::string,model_shape_2d*> *m_interpreter;
        void init();
        void uninit();
    };

	DLL_API void InterpretorModelScript(const char* script_name);
}