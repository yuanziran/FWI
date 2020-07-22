#pragma once
#include"source_func_op.hpp"
namespace FWI
{
     class DLL_API decay_sin_source final : public source_func_op
    {
        public:
        decay_sin_source(const float freq=10.0f,const float timeDelay=0.0f,const float peak=1.0f,const bool isDecaySinSource=true);
        virtual float operator()(const float t)const override;
        virtual void  debug()const override;
        public:
        inline virtual ~decay_sin_source(){}
        public:
        float  decaySinSource(const float time)const;
        float  decayCosSource(const float time)const;
        protected:
        float   peak;
        float   freq;
        float   time_delay;
        bool    is_decay_sin_source;
    };
}