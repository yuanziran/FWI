#include"decay_source.hpp"
#include"fwi_cpu_sys_info.h"
#include<iostream>
using namespace FWI;

decay_sin_source::decay_sin_source(const float freq,const float timeDelay,const float peak,const bool isDecaySinSource)
{
    this->peak=peak;
    this->freq=freq;
    this->time_delay=timeDelay;
    this->is_decay_sin_source=isDecaySinSource;
}

float decay_sin_source::operator()(const float t)const
{
    if(is_decay_sin_source)
    {
        return decaySinSource(t);
    }
    else
    {
        return decayCosSource(t);
    }
}

void decay_sin_source::debug()const
{
    DEBUG_TAG
    DEBUG_TAG
    std::cout<<"decay sin source debug"<<std::endl
    <<"ricker freq               : "<<freq<<std::endl
    <<"time delay                : "<<time_delay<<std::endl
    <<"peak value                : "<<peak<<std::endl
    <<"is decay sin source       : "<<std::boolalpha<<is_decay_sin_source<<std::endl;
    DEBUG_TAG
    DEBUG_TAG
}

float decay_sin_source::decaySinSource(const float time)const
{
    const float t1=2.0f*3.1415926f*freq*(time-time_delay);
    return peak*sinf(t1)*expf(-2.0f*(freq*freq*(time-time_delay)*(time-time_delay)));
}

float decay_sin_source::decayCosSource(const float time)const
{
    const float t1=2.0f*3.1415926f*freq*(time-time_delay);
    return peak*cosf(t1)*expf(-2.0f*(freq*freq*(time-time_delay)*(time-time_delay)));
}