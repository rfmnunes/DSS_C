#ifndef LOG_INCLUDED
#define LOG_INCLUDED

#include <iostream>
#include <fstream>


class Logger
{
public:
Logger(/*const char* filename = "log.log"*/)
{
//logfile = new std::fstream(filename, std::fstream::out);
}
~Logger()
{
logfile->close();
delete logfile;
}

void set_file(const char* filename = "log.log"){
logfile = new std::fstream(filename, std::fstream::out);
}

template<class T>
Logger& operator<<(const T & msg){

	*logfile << msg;
	std::cout << msg;
	return *this;
}


private:
std::fstream* logfile;
};




#endif