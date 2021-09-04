#编译成动态链接库
SET(target_name Forward-CPU)
FILE(GLOB source_list *.h *.c)
add_library(${target_name} SHARED ${source_list})
target_compile_definitions(${target_name} PRIVATE DLL_EXPORT)
target_include_directories(${target_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#提取动态库的头文件
add_subdirectory(Forward-App)

