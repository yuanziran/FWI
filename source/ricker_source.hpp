#pragma once
#include"source_func_op.hpp"
namespace FWI
{
    class DLL_API ricker_source final : public source_func_op
    {
        public:
        ricker_source(const float freq=5.0f,const float timeDelay=0.0f,const float peak=1.0f,const bool isImprovedRicker=false);
        virtual float operator()(const float t)const override;
        virtual void  debug()const override;
        public:
        inline virtual ~ricker_source(){} 
        public:
        float ricker(const float time)const;
        float improvedRicker(const float time)const;
        protected:
        float   peak;
        float   freq;
        float   time_delay;
        bool    is_improved_ricker_source;
    };
}