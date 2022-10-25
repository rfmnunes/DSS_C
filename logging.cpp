

#include "logging.h"
#include "boost_include.h"
#include "log.h"

void clss_log::open_log(registry *reg )
{

	reg_key *k;

	//dbgfl  = dbgfl_in;
	//idbg= idbg_in;
	
	k = get_key(reg, (char*)("DEBUG"), (char*)("DBGLEVEL"));
	if (k)
		idbg = get_int(k);
	//else return -1;
		
	if ((k = get_key(reg, (char*)("DEBUG"), (char*)("DBGFILE"))) != NULL)
		dbgfl= get_string(k);
	
	boost::algorithm::trim(dbgfl);

	//Logger logger(dbgfl.c_str());

	//fl_log.open(dbgfl.c_str());
	//	
	//if (!fl_log.is_open()){
	//	logger << dbgfl << " failed to open the log file!";
	//		
	//}


	//log_string( " Debugging file (1/2/3): " );
	//log_string	(dbgfl.c_str());
	//log_string(" ");
	//log_int(idbg);
	//log_string( "\n");


	//////////////////////////////////////////////////////////////////
//
//
//    logging::add_console_log(std::clog, keywords::format = /*"%TimeStamp%:*/"%Message%");
//
//    // One can also use lambda expressions to setup filters and formatters
//    logging::add_file_log
//    (
//		dbgfl,
//		//"sample.log",
//       // keywords::filter = expr::attr< severity_level >("Severity") >= normal,
//        keywords::format = expr::stream
//         //   << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d, %H:%M:%S.%f")
//         //   << " [" << expr::format_date_time< attrs::timer::value_type >("Uptime", "%O:%M:%S")
//         //   << "] [" << expr::format_named_scope("Scope", keywords::format = "%n (%f:%l)")
//            << "] <" << expr::attr< severity_level >("Severity")
//            << "> " << expr::message
///*
//        keywords::format = expr::format("%1% [%2%] [%3%] <%4%> %5%")
//            % expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d, %H:%M:%S.%f")
//            % expr::format_date_time< attrs::timer::value_type >("Uptime", "%O:%M:%S")
//            % expr::format_named_scope("Scope", keywords::format = "%n (%f:%l)")
//            % expr::attr< severity_level >("Severity")
//            % expr::message
//*/
//    );
//
//    // Also let's add some commonly used attributes, like timestamp and record counter.
//    logging::add_common_attributes();
//    logging::core::get()->add_thread_attribute("Scope", attrs::named_scope());
//
//    BOOST_LOG_FUNCTION();
	

	
	std::cout << " Debugging file (1/2/3): " << dbgfl << " "<< idbg << "\n";


};
