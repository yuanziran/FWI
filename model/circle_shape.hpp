#pragma once
#include"mdoel_shape.hpp"
namespace FWI
{
    class DLL_API circle_shape final : public model_shape_2d
    {
        float  m_x_pos;
        float  m_y_pos;
        float  m_r;
        float  m_value;
        public:
        bool isValid()const;                                               ///<验证参数是否合法
        void parseCommand(std::istringstream& in);                   ///<解析命令行
        bool draw(float* data,const size_t& xdim,const size_t& ydim)const; ///<绘制图像
        std::string getCommandExample()const
        {
            return "circle xpos ypos r value";
        }                              
        void debug()const; 
    }; 
}