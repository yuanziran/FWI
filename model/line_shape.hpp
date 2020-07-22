#pragma once
#include"mdoel_shape.hpp"
namespace FWI
{
    class DLL_API line_shape_2d  : public model_shape_2d
    {
        float   m_x_begin;
        float   m_y_begin;
        float   m_x_end;
        float   m_y_end;
        float   m_value;
        public:
        bool isValid()const;                                             
        void parseCommand(std::istringstream& in);                         
        bool draw(float* data,const size_t& xdim,const size_t& ydim)const;
        std::string getCommandExample()const
        {
            return "line xbegin ybegin xend yend value";
        }                              
		void debug()const;
    }; 
}