#pragma once
#include"fwi_cpu_sys_info.h"
#include"fwi_data_type.hpp"
#include"source_cmd.hpp"
namespace FWI
{
    struct common_cmd_core;
    struct source_buffer
    {
        size_t              load_pos;
        #ifdef MPI_PARALLEL
        size_t              load_offset;
        #endif
        size_t              load_t_pos;
        unsigned int        source_buffer_len;
        SOURCE_LOAD_TYPE    load_type; 
        float*              source_buffer;//using new allocate!
    };

    DLL_API source_buffer  MakeSourceBuffer(const common_cmd_core& info,const source_cmd& cmd);//using new allocate memory!
    DLL_API void  DebugSourceBuffer(const source_buffer& buffer);
}
