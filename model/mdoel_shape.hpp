#pragma once
#include"fwi_cpu_sys_info.h"
#include<string>
#include<sstream>
 namespace FWI
 {
     class DLL_API model_shape_2d
     {
       public:
       virtual bool isValid()const=0;                                               ///<验证参数是否合法
       virtual void parseCommand(std::istringstream& in)=0;                         ///<解析命令行
       virtual bool draw(float* data,const size_t& xdim,const size_t& ydim)const=0; ///<绘制图像
       virtual std::string getCommandExample()const=0;                              ///<给一个命令实例
       virtual void debug()const=0;                                                 ///<打印参数
       virtual ~model_shape_2d(){}    
     };
 }