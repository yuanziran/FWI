#include"model_script.hpp"
#include"model.hpp"
#include"line_shape.hpp"
#include"rect_shape.hpp"
#include"circle_shape.hpp"
#include<iostream>
#include<regex>
using namespace FWI;

model_script::model_script(const char* path)
{
    m_data=nullptr;
    m_xdim=0;
    m_ydim=0;
    m_interpreter=nullptr;
    init();
    std::ifstream in(path);
    assert(in.is_open());
    unsigned int line_no=readScriptHeader(in);
    debug();
    run(in,line_no);
}

model_script::~model_script()
{
    uninit();
    model::WriteModelData2D(out_file_path,m_data,m_xdim,m_ydim);
    if(m_data) delete[] m_data;
}


unsigned int model_script::readScriptHeader(std::ifstream &in)
{
    std::string first_line;
    unsigned int m_line_no=1;
    std::regex r("\\s*#.*");
    bool flag=true;
    while(std::getline(in,first_line))
    {
        std::cout<<"line : "<<m_line_no<<" : "<<first_line<<std::endl;
        m_line_no++;
        flag=std::regex_match(first_line,r);
        //std::cout<<"is flag : "<<std::boolalpha<<flag<<std::endl;
        if(!flag)
        {
            std::istringstream is(first_line);
            is>>m_xdim>>m_ydim>>m_background_value>>out_file_path;
            debug();
            if(m_data) delete[] m_data;
            m_data=new float[m_xdim*m_ydim];
            if(m_data==nullptr)
            {
                std::cerr<<"new fail"<<std::endl;
            }
            for(size_t i=0;i<m_xdim;++i)
            for(size_t j=0;j<m_ydim;++j)
            {
                m_data[i*m_ydim+j]=m_background_value;
            }
            std::cout<<"in here!"<<std::endl;
            break;
        }
    }
    if(in.eof())
    {
        std::cout<<"file end !"<<std::endl;
        exit(EXIT_FAILURE);
    }
    return m_line_no;
}


void model_script:: init()
{
    if(m_interpreter==nullptr)
    {
        m_interpreter=new std::map<std::string,model_shape_2d*>;
        m_interpreter->insert(std::make_pair(std::string("line"),new line_shape_2d()));
        m_interpreter->insert(std::make_pair(std::string("rect"),new rect_shape()));
        m_interpreter->insert(std::make_pair(std::string("circle"),new circle_shape()));
        //m_interpreter->insert(std::make_pair(std::string("triangle"),new triangle_shape()));
    }
}

void model_script::uninit()
{
    if(m_interpreter)
    {
        for(auto iter=m_interpreter->begin();iter!=m_interpreter->end();++iter)
        {
            if(iter->second) delete iter->second;
        }
        delete m_interpreter;
    }
}



void model_script::debug()const
{
   std::cout<<"xdim    : "<<m_xdim<<std::endl
            <<"ydim    : "<<m_ydim<<std::endl
            <<"value   : "<<m_background_value<<std::endl
            <<"output  : "<<out_file_path<<std::endl;
   if(m_interpreter)
    {
        for(auto iter=m_interpreter->begin();iter!=m_interpreter->end();++iter)
        {
			iter->second->debug();
        }
    }
}

void  model_script::run(std::ifstream& in,unsigned int &lineNo)
{
    std::string line;
    std::string command;
    while(std::getline(in,line))
    {
        std::istringstream is(line);
        is>>command;
		lineNo++;
        std::cout << "line : " << lineNo << " command : " << line << std::endl;
        if(command[0]=='#' ||command.empty())//注释行跳过
        {
            continue;
        }
        auto iter=m_interpreter->find(command);
        if(iter!=m_interpreter->end())
        {
            iter->second->parseCommand(is);
            iter->second->draw(m_data,m_xdim,m_ydim);
        }
        else
        {
            std::cerr<<" error line : "<<lineNo<<"error command : "<<command<<std::endl;
            exit(EXIT_FAILURE);
        }
    }   
}

void FWI::InterpretorModelScript(const char* script_name)
{
    model_script app(script_name);
    return;
}