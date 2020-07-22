#include"rect_shape.hpp"
#include<iostream>
using namespace FWI;
bool rect_shape :: isValid()const
{
    auto bound_in=[](float a,float b,float c){
        return (c>=a && c<=b && a<=b);
    };
    if((bound_in(0.0f,m_x_end,m_x_begin) && bound_in(m_x_begin,1.0f,m_x_end) &&
    bound_in(0.0f,m_y_end,m_y_begin) && bound_in(m_y_begin,1.0f,m_y_end))) return true;
    else 
    {
        return false;
    }

}

void rect_shape:: parseCommand(std::istringstream& in)
{
    in>>m_x_begin>>m_y_begin>>m_x_end>>m_y_end>>m_value;
}

bool rect_shape :: draw(float* data,const size_t& xdim,const size_t& ydim)const
{
    if(!isValid()) return false;
    size_t xmin=size_t(m_x_begin*xdim);
    size_t ymin=size_t(m_y_begin*ydim);
    size_t xmax=size_t(m_x_end*xdim);
    size_t ymax=size_t(m_y_end*ydim);
    for(size_t i=xmin;i<xmax;++i)
    for(size_t j=ymin;j<ymax;++j)
    {
        data[i*ydim+j]=m_value;
    }
    return true;
}

void rect_shape:: debug()const
{
    DEBUG_TAG
    DEBUG_TAG
    std::cout<<"rect : "<<std::endl
    <<"xbegin ="<<m_x_begin<<std::endl
    <<"ybegin ="<<m_y_begin<<std::endl
    <<"xend   ="<<m_x_end<<std::endl
    <<"yend   ="<<m_y_end<<std::endl
    <<"value  ="<<m_value<<std::endl;
    DEBUG_TAG
    DEBUG_TAG
}