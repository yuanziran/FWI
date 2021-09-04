/**
 * @file fwi_cmd.hpp
 * @author Zhenghong Guo (you@domain.com)
 * @brief 解析命令行
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once 
#include<string>
#include<map>
#include<iostream>
#include<iterator>

typedef std::map<std::string,std::string> fwi_cmd_t;

static inline void ParseCommandLine(fwi_cmd_t& command, int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string cmd(argv[i]);
        auto pos = cmd.find_first_of('=');
        if (pos == std::string::npos || pos == 0)//����key=value
        {
            std::string s;
            command[cmd] = s;
        }
        else
        {
            command[cmd.substr(0, pos)] = cmd.substr(pos + 1, cmd.size() - pos - 1);
        }
    }
}

static inline void DebugCommandLine(const fwi_cmd_t& command)
{
	for (auto iter = command.begin(); iter != command.end(); ++iter)
	{
		std::cout << "key := " << iter->first << " value := " << iter->second << std::endl;
	}
}

static inline std::string GetCmdValue(const fwi_cmd_t& command, const std::string& key, int flag, const std::string& opt = std::string())
{
    auto iter = command.find(key);
    if (iter == command.end())
    {
        if (flag == 0)//����Ҫ��
        {
            std::cerr << "Key : " << key << " Must Give!" << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << "key : " << key << " Default Value : " << opt << std::endl;
            return opt;
        }
    }
    else
    {
        return iter->second;
    }
}


static inline void GetCmdValueEx(const fwi_cmd_t& command, const std::string& key, int flag, size_t& out, const size_t opt = 0)
{
    auto iter = command.find(key);
    if (iter == command.end())
    {
        if (flag == 0)//����Ҫ��
        {
            std::cerr << "Key : " << key << " Must Give! " << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << "key : " << key << " Default Value : " << opt << std::endl;
            out = opt;
        }
    }
    else
    {
        if (iter->second.empty())
        {
            out = opt;
        }
        else
        {
            sscanf(iter->second.c_str(), "%llu", &out);
        }
    }
}

static inline void GetCmdValueEx(const fwi_cmd_t& command, const std::string& key, int flag, float& out, const float opt = 0.0f)
{
    auto iter = command.find(key);
    if (iter == command.end())
    {
        if (flag == 0)//����Ҫ��
        {
            std::cerr << "Key : " << key << " Must Give!" << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << "key : " << key << " Default Value : " << opt << std::endl;
            out = opt;
        }
    }
    else
    {
        if (iter->second.empty())
        {
            out = opt;
        }
        else
        {
            sscanf(iter->second.c_str(), "%f", &out);
        }
    }
}


static inline void GetCmdValueEx(const fwi_cmd_t& command, const std::string& key, int flag, int& out, const int opt = 0)
{
    auto iter = command.find(key);
    if (iter == command.end())
    {
        if (flag == 0)//����Ҫ��
        {
            std::cerr << "Key : " << key << " Must Give!" << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << "key : " << key << " Default Value : " << opt << std::endl;
            out = opt;
        }
    }
    else
    {
        if (iter->second.empty())
        {
            out = opt;
        }
        else
        {
            sscanf(iter->second.c_str(), "%d", &out);
        }
    }
}

static inline void GetCmdValueEx(const fwi_cmd_t& command, const std::string& key, int flag, std::string& out, const std::string opt = std::string())
{
    auto iter = command.find(key);
    if (iter == command.end())
    {
        if (flag == 0)//����Ҫ��
        {
            std::cerr << "Key : " << key << " Must Give!" << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << "key : " << key << " Default Value : " << opt << std::endl;
            out = opt;
        }
    }
    else
    {
        if (iter->second.empty())
        {
            out = opt;
        }
        else
        {
            out = iter->second;
        }
    }
}

