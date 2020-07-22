#pragma once
#include<string>
namespace FWI
{
	int constexpr  MAX_ARGC = 1024;
    /**
     * @brief 命令行辅助类
    */
	class arg_tool final
	{
		public:
		arg_tool(const char* command,const char* tag=" ");
		arg_tool(const std::string& cmd,const char* tag=" "):arg_tool(cmd.c_str(),tag){} 
		inline int getArgc(void)const{return argc;}
		inline char** getArgv(void)const{return (char**)argv;}
		void   debug()const;
		void   addArgc(const char* opt);
		~arg_tool();	
 		protected:
		char* argv[MAX_ARGC];
		int	  argc;
		int   _argc;
		protected:
		arg_tool();
		void split(const char* optarg,const char* tag=" ");
	};
}