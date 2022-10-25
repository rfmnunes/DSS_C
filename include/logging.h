#ifndef LOGGING_INCLUDED
#define LOGGING_INCLUDED


#include <iostream>		// logger
#include <fstream>		// fl_log
#include <sstream>		// <<

#include "registry.h"
#include "log.h"

class clss_log
{	// int dbg to krige?

public:
	
	int idbg;
	std::string dbgfl;

	std::ofstream fl_log;

	void open_log(registry *reg );
	

	void log_string( const char * logstring );

	void log_int( const int  logint );


}; 




extern clss_log g_log;
extern Logger logger;
//
#endif
