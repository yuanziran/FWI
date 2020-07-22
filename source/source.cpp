#include"source.hpp"
#include"common_cmd.hpp"
#include"source_func_op.hpp"
#include"ricker_source.hpp"
#include"decay_source.hpp"
#include<iostream>
using namespace FWI;

source_buffer  FWI::MakeSourceBuffer(const common_cmd_core& info,const source_cmd& cmd)
{
    source_buffer buffer;
    unsigned int common_dim=(info.cdo_order+info.abc_layers);
    size_t zdim=info.input_dim.z_dim+(common_dim<<1);
    if(info.input_dim.y_dim==0)
    {
        buffer.load_pos=(cmd.load_x_pos+common_dim)*zdim+(cmd.load_z_pos+common_dim);
        #ifdef MPI_PARALLEL
        buffer.load_offset=(cmd.shot_x_interval)*zdim+cmd.shot_z_interval;
        #endif
    }
    else
    {
        size_t ydim=info.input_dim.y_dim+(common_dim<<1);
        buffer.load_pos=((cmd.load_x_pos+common_dim)*ydim+cmd.load_y_pos+common_dim)*zdim+cmd.load_z_pos+common_dim;
        #ifdef MPI_PARALLEL
        buffer.load_offset=((cmd.shot_x_interval)*ydim+cmd.shot_y_interval)*zdim+cmd.shot_z_interval;
        #endif
    }
    buffer.load_t_pos=cmd.load_t_begin;
    buffer.source_buffer_len=cmd.load_t_end-cmd.load_t_begin;
    buffer.source_buffer=new float[buffer.source_buffer_len];
    source_func_op *op=nullptr;
    switch(cmd.source_func_type)
    {
        case SOURCE_FUNC_TYPE::RICKER_TYPE:
            op=new ricker_source(cmd.main_freq,cmd.time_delay,cmd.max_amp,false);
            break;
        case SOURCE_FUNC_TYPE::IMPROVED_RICKER_TYPE:
            op=new  ricker_source(cmd.main_freq,cmd.time_delay,cmd.max_amp,true);
            break;
        case SOURCE_FUNC_TYPE::DECAY_SIN_TYPE:
            op=new decay_sin_source(cmd.main_freq,cmd.time_delay,cmd.max_amp,true);
            break;
        case SOURCE_FUNC_TYPE::DECAY_COS_TYPE:
            op=new decay_sin_source(cmd.main_freq,cmd.time_delay,cmd.max_amp,false);
            break;
        default:
            std::cerr<<"bad type "<<(int&)cmd.source_func_type<<std::endl;
            exit(EXIT_FAILURE);
    }
    for(unsigned int i=0;i<buffer.source_buffer_len;++i)
    {
        buffer.source_buffer[i]=(*op)(i*info.sample_dim.dt);
    }
    delete op;
    buffer.load_type=cmd.source_load_type;
    return buffer;
}

void FWI::DebugSourceBuffer(const source_buffer& buffer)
{
    DEBUG_TAG
    DEBUG_TAG
    std::cout<<"Source Buffer Debug!"<<std::endl
    <<"load pos     : "<<buffer.load_pos<<std::endl
    #ifdef MPI_PARALLEL
    <<"load offset  : "<<buffer.load_offset<<std::endl
    #endif
    <<"load t begin : "<<buffer.load_t_pos<<std::endl
    <<"buffer len   : "<<buffer.source_buffer_len<<std::endl;
    for(unsigned int i=0;i<buffer.source_buffer_len;++i)
    {
        std::cout<<"[ "<<i<<" ]="<<buffer.source_buffer[i]<<std::endl;
    }
    DEBUG_TAG
    DEBUG_TAG
}