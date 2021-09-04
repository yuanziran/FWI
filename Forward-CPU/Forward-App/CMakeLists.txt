file(GLOB source_list *.cpp)
foreach(file_name ${source_list})
    STRING(REGEX REPLACE ".+/(.+)\\..*" "\\1" project_name ${file_name})#提取到文件名
    add_executable(${project_name} ${file_name})
    target_link_libraries(${project_name} PRIVATE Forward-CPU)
endforeach()