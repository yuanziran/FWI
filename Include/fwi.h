/**
 * @file fwi.h
 * @author Zhenghong Guo (you@domain.com)
 * @brief  全局的宏函数以及辅助函数定义
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once
#include<assert.h>
#include<stddef.h>
#include<stdint.h>
#include<string.h>
#include<stdlib.h>
#include<memory.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<malloc.h>
#include<float.h>
#include<time.h>
#include<stdbool.h>
#define STR(call)       #call


#define DEBUG_TAG printf("*************************************************************************\n");

#define PRINT_CSTRING(p) (p?p:"null")

#define RUN_AT(info)    \
{\
    printf("where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
    printf("%s\n",info);\
}


static inline char* AllocateAndInitCString(const char* str)
{
    char* p=(char*)malloc((strlen(str)+1));
    assert(p);
    strncpy(p,str,strlen(str)+1);
    return p;
}

static inline void DestroyCString(char* str)
{
    if(str)
    {
        free(str);
    }
}


static inline bool IsFloatEqual(const float a,const float b)
{
    if(fabsf(a-b)<FLT_EPSILON)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static inline bool IsDoubleEqual(const double a,const float b)
{
    if(fabs(a-b)<DBL_EPSILON)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static inline bool IsLongDoubleEqual(const long double a,const long double b)
{
    if(fabsl(a-b)<LDBL_EPSILON)
    {
        return true;
    }
    else
    {
        return false;
    }
}


//计时工具
#define TIME_BEGIN(name)                                                    \
{                                                                           \
    time_t start,end;                                                       \
    start=time(NULL);                               


#define TIME_END(taskName)                                                  \
    end=time(NULL);                                                         \
    printf("%s time cost : %f\n",taskName,difftime(end,start));             \
}
/**
 *  time_t start,end;
 *  start=time(NULL);
 *   *
 *   *
 *   *
 *   *
 *   *
 *  end=time(NULL)
 *  printf(time=%d\n",difftime(end,start));
 * 
*/


#ifdef _WIN32
# include <io.h>
# include <fcntl.h>
# define SET_BINARY_MODE(handle) setmode(_fileno(handle), O_BINARY)
#else
# define SET_BINARY_MODE(handle) ((void)0)
#endif

#if defined linux
#include<errno.h>
#ifdef OPENMP_PARALLEL
#include<omp.h>
#endif
#define STRCASECMP strcasecmp
#define CHECK(call)         \
{\
    if(!(call))\
    {\
        printf("call name : %s. where in file:%s line :%d func:%s\n Error : %s error:%d\n",STR(call),__FILE__,__LINE__,__FUNCTION__,strerror(errno),errno);\
        exit(EXIT_FAILURE);\
    }\
}

#include<sys/time.h>

inline double CpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec+(double)tp.tv_usec*1.e-6);
}

#define BUG_ON(format,args...) \
{\
     printf("fatal error : where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
    printf(format,##args);\
    exit(EXIT_FAILURE);\
}

#define WARNING_ON(format,args...) \
{\
    printf("fatal error : where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
    printf(format,##args);\
}

#define DLL_API  

#elif defined _WIN32
#include<windows.h>
#include<strsafe.h>
inline void
Win32Perror(LPTSTR lpszFunction)
{
    LPVOID lpMsgBuf;                                                                                                                                                                
    LPVOID lpDisplayBuf;                                                                                                                                                            
    DWORD  dw=GetLastError();                                                                                                                                                       
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER|FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS,NULL,dw,MAKELANGID(LANG_NEUTRAL,SUBLANG_DEFAULT),(LPTSTR)&lpMsgBuf,0,NULL);
    lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT,(lstrlen((LPCTSTR)lpMsgBuf)+lstrlen((LPCTSTR)lpszFunction)+40)*sizeof(TCHAR));                                                  
    StringCchPrintf((LPTSTR)lpDisplayBuf,LocalSize(lpDisplayBuf),TEXT("%s failed with error %d: %s"),STR(call), dw, lpMsgBuf);                                                           
    MessageBox(NULL,(LPCTSTR)lpDisplayBuf,TEXT("Error"), MB_OK);                                                                                                                    
    LocalFree(lpMsgBuf); LocalFree(lpDisplayBuf);                                                                                                                                   
    exit(dw);  
}

#define CHECK(call)\
{\
    if(!(call))\
    {\
        Win32Perror(STR(call));\
    }\
}

#define STRCASECMP _stricmp

#define BUG_ON(format,...) \
{\
    printf("fatal error : where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
    printf(format,__VA_ARGS__);\
    exit(EXIT_FAILURE);\
}

#define WARNING_ON(format,...)\
{\
    printf("warning : where in file :%s in func %s at line : %d\n",__FILE__,__func__,__LINE__);\
    printf(format,__VA_ARGS__);\
}



#ifdef DLL_EXPORT
#define DLL_API  __declspec(dllexport)
#else
#define DLL_API  __declspec(dllimport)
#endif

#endif
