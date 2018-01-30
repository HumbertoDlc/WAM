// C++
#include <cstdio>  // For remove()
#include <fstream>
#include <math.h>  
#include <iostream>
#include <string>
#include <cstdlib>  // For std::atexit()
#include <algorithm>
// BARRETT
#include <boost/tuple/tuple.hpp>
#include <barrett/log.h>
#include <barrett/units.h>
#include <curses.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/os.h>
#include <barrett/math.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/standard_main_function.h>
// LABJACK
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "ljacklm.h"

using namespace barrett;
using namespace std;
void waitForEnter() 
{
	std::string line;
	std::getline(std::cin, line);
}

/*void handleError(long error, char *functionName)
{
	char errorString[50];

	if(error)
	{
		GetErrorString(error, errorString);
		printf("%s error %ld: %s\n", functionName, error, errorString);
		exit(1);
	}
}*/
int mainLJ (double ENC_1, double ENC_2, double ENC_3, double ENC_4)
{
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////LABJACK/////////////////////////////////////////////////////
	long idnum12 = 0;//100038943;//-1;  //Using first found U12
	long idnum34 = 1;//100039743;	
	long demo = 0;  //Normal operations
	long stateIO = 0;  //Output states for IO0-IO3

	//AOUpdate specific paramters
	long trisD = 0;  //Directions for D0-D15
	long trisIO = 0;  //Directions for IO0-IO3
	long stateD = 0;  //Output states for D0-D15
	long updateDigital = 0;  //Tris and state values are only being read
	long resetCounter = 0;  //Not resetting counter
	unsigned long count = 0;  //Returned current count value
	float analogOut0 = 0;  //Voltage for AO0 DEVICE 1
	float analogOut1 = 0;  //Voltage for AO1 DEVICE 1
	float analogOut2 = 0;  //Voltage for AO0 DEVICE 2
	float analogOut3 = 0;  //Voltage for AO1 DEVICE 2

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////LABJACK/////////////////////////////////////////////////////
analogOut0=5*(ENC_1+2.6695)/(2.7344+2.6695);
analogOut1=5*(ENC_2+1.9888)/(1.9828+1.9888);
analogOut2=5*(ENC_3+2.6203)/(2.9593+2.6203);
analogOut3=5*(ENC_4+0.8231)/(3.1229+0.8231);
	
	if (analogOut0>=5)
	{analogOut0=5;}
	if (analogOut0<=0)
	{analogOut0=0;}
	if (analogOut1>=5)
	{analogOut1=5;}
	if (analogOut1<=0)
	{analogOut1=0;}
	if (analogOut2>=5)
	{analogOut2=5;}
	if (analogOut2<=0)
	{analogOut2=0;}
	if (analogOut3>=5)
	{analogOut3=5;}
	if (analogOut3<=0)
	{analogOut3=0;}
	
	AOUpdate(&idnum12, demo, trisD, trisIO, &stateD,
					&stateIO, updateDigital, resetCounter, &count, analogOut0,
					analogOut1);
	AOUpdate(&idnum34, demo, trisD, trisIO, &stateD,
					&stateIO, updateDigital, resetCounter, &count, analogOut2,
					analogOut3);

return 0;
}



template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) 
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	jp_type zeroMove;

	
	wam.gravityCompensate(true);
	
	bool going = true;
	//int val;
	//wam.moveTo(jp_type(0.0));

	jp_type J_arr;
	
	std::cout << "Entering to the loop.\n";
	//waitForEnter();

	while (going)
	{
		J_arr=wam.getJointPositions();
		mainLJ(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);
		
	}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}
