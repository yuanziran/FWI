#pragma once
#include"fwi_cpu_sys_info.h"
namespace FWI
{ 
    class DLL_API source_func_op
    {
        public:
        virtual float operator()(const float t)const=0;///<震源函数
        virtual void  debug()const=0;///<调试接口
        public:
        inline virtual ~source_func_op(){}
    };
}