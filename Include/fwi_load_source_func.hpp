/**
 * @file fwi_load_source_func.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 常用的震源函数
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once  
#include<math.h>
//雷克子波震源
static inline float RickerSource(const float main_freq,const float time,const float time_delay=0.0f,const float max_amp=1.0f)
{
    constexpr float pi=3.1415926f;
    float t=pi*main_freq*(time-time_delay);
    float t_2=t*t;
    return max_amp*(1.0f-2.0f*t_2)*expf(-t_2);
}


