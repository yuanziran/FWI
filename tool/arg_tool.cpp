#include"arg_tool.hpp"
#include<cstring>
#include<cassert>
#include<cstdlib>
#include<cstdio>
#include<cstddef>
#include<memory.h>
using namespace FWI;

arg_tool::arg_tool()
{
    argc=0;
	_argc=0;
	memset(argv,0,sizeof(argv));
}

arg_tool::arg_tool(const char* command,const char* tag):arg_tool()
{
    split(command,tag);
}

arg_tool::~arg_tool()
{
    for(int i=argc;i>_argc;--i)
	{
		delete[] argv[i];	
	}
}

void arg_tool::debug()const
{
    for(int i=0;i<argc;++i)
		printf("argv[%d] : %s\n",i,argv[i]);
}

void arg_tool::addArgc(const char* opt)
{
    size_t len=strlen(opt);
	char* p=new char[len+1];
    strncpy(p,opt,len);
	argv[argc]=p;
	++argc;
    assert(argc<MAX_ARGC);
}

void arg_tool::split(const char* optarg,const char* tag)
{
    char *p=strtok((char*)optarg,tag);
	while(p!=NULL)
	{
 		argv[argc]=p;
		argc++;
		p=strtok(NULL,tag);
	}	
	_argc=argc;
}