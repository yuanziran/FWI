#include"circle_shape.hpp"
#include<iostream>
using namespace FWI;

bool circle_shape :: isValid()const
 {
    if(m_r>=0.0f) return true;
    else return false;
 }

void circle_shape :: parseCommand(std::istringstream& in)
{
    in>>m_x_pos>>m_y_pos>>m_r>>m_value;
}

bool circle_shape :: draw(float* data,const size_t& xdim,const size_t& ydim)const
{
    if(isValid())
    {
        float xmin=m_x_pos-m_r;
        float ymin=m_y_pos-m_r;
        float xmax=m_x_pos+m_r;
        float ymax=m_y_pos+m_r;
        if(xmin<0.0f) xmin=0.0f;
        if(ymin<0.0f) ymin=0.0f;
        if(xmax>1.0f) xmax=1.0f;
        if(ymax>1.0f) ymax=1.0f;
        size_t xlow=size_t(xmin*xdim);
        size_t xheight=size_t(xmax*xdim);
        size_t ylow=size_t(ymin*ydim);
        size_t yheight=size_t(ymax*ydim);
        float r_2=m_r*m_r-m_x_pos*m_x_pos-m_y_pos*m_y_pos;
        auto fun=[this,r_2](float x,float y)
        {
            return ((x*x+y*y-2.0f*(x*m_x_pos+y*m_y_pos))<=r_2);
        };
        for(size_t i=xlow;i<xheight;++i)
        for(size_t j=ylow;j<yheight;++j)
        {
            float x=float(i)/xdim;
            float y=float(j)/ydim;
            if(fun(x,y))
            {
                data[i*ydim+j]=m_value;
            }
        }
        return true;
    }
    return false;
}

void circle_shape :: debug()const
{
    DEBUG_TAG
    DEBUG_TAG
    std::cout<<"circle shape : "<<std::endl
    <<"xpos ="<<m_x_pos<<std::endl
    <<"ypos ="<<m_y_pos<<std::endl
    <<"r    ="<<m_r<<std::endl
    <<"value="<<m_value<<std::endl;
    DEBUG_TAG
    DEBUG_TAG
}
