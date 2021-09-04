/**
 * @file fwi_memory.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 内存分配以及相应的文件处理
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once  
#include"fwi.h"
#ifdef CUDA_PARALLEL 
#include"fwi_gpu.h"
#endif
static inline float* FWIAllocateMemory(const size_t num,float value=0.0f)
{
    if(num==0) return nullptr;
    float* ptr=new float[num];
    assert(ptr!=nullptr);
    for(size_t i=0;i<num;++i)
    {
        ptr[i]=value;
    }
    return ptr;
}


static inline void FWIFreeMemory(float** ptr)
{
    if(*ptr!=nullptr)
    {
        delete[] *ptr;
    }
    *ptr=nullptr; 
}

#define DebugFreeMemory(ptr) \
{\
	if(*ptr==nullptr)\
	{\
		printf("%s is nullptr\n",STR(ptr));\
		printf("where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
	}\
	else\
	{\
		FWIFreeMemory(ptr);\
	}\
}

#define DebugFreeMemoryEx(ptr) \
{\
	if(ptr==nullptr)\
	{\
		printf("%s is nullptr\n",STR(ptr));\
		printf("where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
	}\
	else\
	{\
		FWIFreeMemory(&ptr);\
	}\
}


static inline float* const FWIMemsetMemory(float* const memory,const size_t num,const float value=0.0f)
{
    for(size_t i=0;i<num;++i)
    {
        memory[i]=value;
    }
    return memory;
}


#ifdef CUDA_PARALLEL

static inline float* FWIGPUAllocateMemory(const size_t num,float value=0.0f)
{
    if(num==0) return nullptr;
    float *ptr;
    CUDA_CHECK(cudaMallocManaged((void**)&ptr, sizeof(float)*num));
    for(size_t i=0;i<num;++i)
    {
        ptr[i]=value;
    }
    return ptr;
}


static inline void FWIFreeMemory(float** ptr)
{
    if(*ptr==nullptr)
    {
        printf("%s is nullptr!\n",STR(ptr));
    }
    else
    {
       CUDA_CHECK(cudaFree(*ptr));
    }
    *ptr=nullptr;
}
#endif