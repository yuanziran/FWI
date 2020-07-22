set(target_name tool)
file(GLOB SOURCES *.c *.cpp *.h *.hpp)
add_library(${target_name} STATIC  ${SOURCES})
target_include_directories(${target_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})