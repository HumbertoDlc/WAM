#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()
#include <iostream> 
#include <string>

#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>

#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>

#include<ctime>
#include <time.h>

using namespace barrett;
using systems::connect;


bool validate_args(int argc, char** argv) {
	return true;
}


template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	std::string filestr;
	const char *filename;
	
	//time_t _tm=time(NULL );
	//struct tm * curtime=localtime(&_tm);
	printf("Please enter file name to log data to: ");
	std::getline(std::cin,filestr);
	filename = filestr.c_str();
	
	//std::string filename;//char *fileptr;					asctime(curtime)
	//printf("Please enter file name to log data to: ");
	//

	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}


	systems::Ramp time(pm.getExecutionManager(), 1.0);

	systems::TupleGrouper<double, jp_type, jv_type, jt_type, cp_type, Eigen::Quaterniond> tg;
	connect(time.output, tg.template getInput<0>());
	connect(wam.jpOutput, tg.template getInput<1>());
	connect(wam.jvOutput, tg.template getInput<2>());
	connect(wam.jtSum.output, tg.template getInput<3>());
	connect(wam.toolPosition.output, tg.template getInput<4>());
	connect(wam.toolOrientation.output, tg.template getInput<5>());

	typedef boost::tuple<double, jp_type, jv_type, jt_type, cp_type, Eigen::Quaterniond> tuple_type;
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);

	time.start();
	connect(tg.output, logger.input);
	printf("Logging started.\n");


	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);


	logger.closeLog();
	printf("Logging stopped.\n");

	log::Reader<tuple_type> lr(tmpFile);
	lr.exportCSV(filename);
	printf("Output written to %s.\n", filename);
	std::remove(tmpFile);

	return 0;
}
