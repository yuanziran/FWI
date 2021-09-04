/**
 * @file fwi_io.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 一些有用文件读写工具
 * @version 0.1
 * @date 2021-09-03
 * 
 * @copyright Copyright (c) 2021
 * 
 */


#pragma once
#include<fstream>
#include<iostream>
#include<assert.h>
#include<cstdlib>
#include<string>

static inline void Read2DBin(const char* path, float* const buffer, const int n1, const int n2)
{
    std::ifstream in(path, std::ios_base::binary);
	if (!in.is_open())
	{
		std::cerr << path << " cant open for read!" << std::endl;
		exit(-1);
	}
    in.read((char*)buffer, sizeof(float) * n1 * n2);
    in.close();
}


static inline void Write2DBin(const char* path, const float* const buffer, const int n2, const int n1)
{
    std::string name = path;
    name += std::to_string(n2);
    name += "-";
    name += std::to_string(n1);
    name += ".bin";
    std::ofstream out(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
    assert(out.is_open());
    out.write((const char*)buffer, sizeof(float) * n1 * n2);
    out.close();
}

static inline void Write2DBinEx(const char* path,const char* suffix, const float* const buffer, const int n2, const int n1)
{
    std::string name = path;
    name +=suffix;
    name += std::to_string(n2);
    name += "-";
    name += std::to_string(n1);
    name += ".bin";
    std::ofstream out(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
    assert(out.is_open());
    out.write((const char*)buffer, sizeof(float) * n1 * n2);
    out.close();
}


static inline void Read3DBin(const char* path, float* const buffer, const int n1, const int n2,const int n3)
{
	std::ifstream in(path, std::ios_base::binary);
	if (!in.is_open())
	{
		std::cerr << path << " cant open for read!" << std::endl;
		exit(-1);
	}
	in.read((char*)buffer, sizeof(float) * n1 * n2*n3);
	in.close();
}


static inline void Write3DBin(const char* path, const float* const buffer,const int n3,const int n2, const int n1)
{
	std::string name = path;
	name += std::to_string(n3);
	name += "-";
	name += std::to_string(n2);
	name += "-";
	name += std::to_string(n1);
	name += ".bin";
	std::ofstream out(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
	assert(out.is_open());
	out.write((const char*)buffer, sizeof(float) * n1 * n2*n3);
	out.close();
}

static inline void Write3DBinEx(const char* path, const char* suffix, const float* const buffer,const int n3,const int n2, const int n1)
{
	std::string name = path;
	name += suffix;
	name += std::to_string(n3);
	name += "-";
	name += std::to_string(n2);
	name += "-";
	name += std::to_string(n1);
	name += ".bin";
	std::ofstream out(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
	assert(out.is_open());
	out.write((const char*)buffer, sizeof(float) * n1 * n2);
	out.close();
}


