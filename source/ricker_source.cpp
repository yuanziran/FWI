#include"ricker_source.hpp"
#include"fwi_cpu_sys_info.h"
#include<iostream>
using namespace FWI;

ricker_source::ricker_source(const float freq,const float timeDelay,const float peak,const bool isImprovedRicker)
{
    this->peak=peak;
    this->freq=freq;
    this->time_delay=timeDelay;
    this->is_improved_ricker_source=isImprovedRicker;
}

float ricker_source::operator()(const float t)const
{
    if(is_improved_ricker_source)
    {
        return improvedRicker(t);
    }
    else
    {
        return ricker(t);  
    }
}

void ricker_source::debug()const
{
    DEBUG_TAG
    DEBUG_TAG
    std::cout<<"ricker source debug"<<std::endl
    <<"ricker freq               : "<<freq<<std::endl
    <<"time delay                : "<<time_delay<<std::endl
    <<"peak value                : "<<peak<<std::endl
    <<"is improved ricker source : "<<std::boolalpha<<is_improved_ricker_source<<std::endl;
    DEBUG_TAG
    DEBUG_TAG
}

float ricker_source::ricker(const float time)const
{
    float t=3.1415926f*freq*(time-time_delay);
    float t_2=t*t;
    return peak*(1.0f-2.0f*(t_2)*expf(-t_2));
}

float ricker_source::improvedRicker(const float time)const
{
     const float a=3.0f*sqrtf(6.0f);
     const float b=13.0f*expf(-13.5f);
     constexpr float pi_2=2*3.1415926f;
     float phase=freq/(pi_2)*time-a;
     float t=0.25f*phase*phase;
     return ((t-0.5f)*expf(-t)-b)/(0.5f+b);
}
