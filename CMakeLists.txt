cmake_minimum_required(VERSION 3.14)
project(FWI VERSION 1.0.0)
set(CMAKE_INSTALL_PREFIX  FWI)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
#设置程序编译后的路径
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "Current Platform : Linux")
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/bin)
    set(LIBRARY_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/lib)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Current Platform : Windows")
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/bin)
    set(LIBRARY_OUTPUT_PATH ${CMAKE_INSTALL_PREFIX}/bin)
endif()


#编译器编译选项控制
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(STATUS "using clang")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "using gcc")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message(STATUS "using icc")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "using visual studio c++")
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    add_compile_options(/FS)
endif()



#查找某些并行库
option(USE_CUDA "USE CUDA PARALLEL" OFF)
if(USE_CUDA)
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        add_definitions(-DCUDA_PARALLEL)
        enable_language(CUDA)
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_50)
    else(CUDA_FOUND)
        message("Can't find Cuda Library")
    endif(CUDA_FOUND)
endif(USE_CUDA)
option(USE_MPI "USE MPI PARALLEL" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        message(STATUS "${MPI_INCLUDE_PATH}")
        message(STATUS "${MPI_C_LIBRARIES}")
        add_definitions(-DMPI_PARALLEL)
        include_directories(${MPI_INCLUDE_PATH})
    endif(MPI_FOUND)
endif(USE_MPI)
option(USE_OPENMP "USE OPENMP PARALLEL" OFF)
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        endif()
        add_definitions(-DOPENMP_PARALLEL)
    else()
        message(STATUS "Can't find OpenMP Library")
    endif()
endif(USE_OPENMP)


include_directories(Include)
add_subdirectory(Forward-CPU)

#程序安装过程
install(TARGETS Forward-CPU AcSg-Forward AcSg-Forward3 ElSg-Forward ElSg-Forward3
LIBRARY DESTINATION lib
RUNTIME DESTINATION bin 
PUBLIC_HEADER DESTINATION include
)