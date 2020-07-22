#include"line_shape.hpp"
#include<iostream>
using namespace FWI;

bool line_shape_2d :: isValid()const
{
    auto bound_in=[](float a,float b,float c){
        return (c>=a && c<=b && a<=b);
    };
    if((bound_in(0.0f,m_x_end,m_x_begin) && bound_in(m_x_begin,1.0f,m_x_end) &&
    bound_in(0.0f,m_y_end,m_y_begin) && bound_in(m_y_begin,1.0f,m_y_end)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

 void line_shape_2d::parseCommand(std::istringstream& in)
 {
     in>>m_x_begin>>m_y_begin>>m_x_end>>m_y_end>>m_value;
 }

 bool line_shape_2d::draw(float* data,const size_t& xdim,const size_t& ydim)const
 {
    if(!isValid()) return false;
    size_t x_begin=size_t(m_x_begin*xdim);
    size_t y_begin=size_t(m_y_begin*ydim);
    size_t x_end=size_t(m_x_end*xdim);
    size_t y_end=size_t(m_y_end*ydim);
    if(x_begin==x_end)
    {
        for(size_t i=y_begin;i<y_end;++i)
        {
            data[x_begin*ydim+i]=m_value;
        }
    }
    else if(y_begin==y_end)
    {
        for(size_t i=x_begin;i<y_end;++i)
        {
            data[i*ydim+y_begin]=m_value;
        }
    }
    else
    {
        float ration=(m_y_end-m_y_begin)/(m_x_end-m_x_begin);
        for(size_t i=x_begin;i<x_end;++i)
        {
            data[i*ydim+size_t((i-x_begin)*ration)+y_begin]=m_value;
        }
    }
    return true;
 }

void line_shape_2d::debug()const
{
    std::cout<<"line : "<<std::endl
             <<"xbegin="<<m_x_begin<<std::endl
             <<"ybegin="<<m_y_begin<<std::endl
             <<"xend  ="<<m_x_end<<std::endl
             <<"yend  ="<<m_y_end<<std::endl
             <<"value ="<<m_value<<std::endl;
}