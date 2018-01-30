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
#include <math.h>   
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>
#include "ljacklm.h"

using namespace barrett;
using namespace std;
void waitForEnter() 
{
	std::string line;
	std::getline(std::cin, line);
}

float Angle_X (double ENC_1, double ENC_2, double ENC_3, double ENC_4)
{
	float Z4[3];
	float angle = 0.0;
	Z4[0]=-cos(ENC_3)*sin(ENC_1)-cos(ENC_1)*cos(ENC_2)*sin(ENC_3);
	Z4[1]=cos(ENC_1)*cos(ENC_3)-cos(ENC_2)*sin(ENC_1)*sin(ENC_3);
	Z4[2]=sin(ENC_2)*sin(ENC_3);
	angle=acos(-Z4[2]/sqrt((Z4[0]*Z4[0])+(Z4[1]*Z4[1])+(Z4[2]*Z4[2])));
	//std::cout<<"Angle X: "<<angle<<"\n\n\n\n\n\n\n\n";
return angle;
}

float Angle_Y (double ENC_1, double ENC_2, double ENC_3, double ENC_4)
{
	float Z4[3];
	float angle = 0.0;
	Z4[0]=-cos(ENC_4)*(sin(ENC_1)*sin(ENC_3)-cos(ENC_1)*cos(ENC_2)*cos(ENC_3))-cos(ENC_1)*sin(ENC_2)*sin(ENC_4);
	Z4[1]=cos(ENC_4)*(cos(ENC_1)*sin(ENC_3)+cos(ENC_2)*cos(ENC_3)*sin(ENC_1))-sin(ENC_1)*sin(ENC_2)*sin(ENC_4);
	Z4[2]=-cos(ENC_2)*sin(ENC_4)-cos(ENC_3)*cos(ENC_4)*sin(ENC_2);
	angle=acos(Z4[2]/sqrt((Z4[0]*Z4[0])+(Z4[1]*Z4[1])+(Z4[2]*Z4[2])));
	//std::cout<<"Angle Y: "<<angle<<"\n\n\n\n\n\n\n\n";
return angle;
}

float Angle_Z (double ENC_1, double ENC_2, double ENC_3, double ENC_4)
{
	float Z4[3];
	float angle = 0.0;
	Z4[0]=cos(ENC_1)*cos(ENC_4)*sin(ENC_2)-sin(ENC_4)*(sin(ENC_1)*sin(ENC_3)-cos(ENC_1)*cos(ENC_2)*cos(ENC_3));
	Z4[1]=sin(ENC_4)*(cos(ENC_1)*sin(ENC_3)+cos(ENC_2)*cos(ENC_3)*sin(ENC_1))+cos(ENC_4)*sin(ENC_1)*sin(ENC_2);
	Z4[2]=cos(ENC_2)*cos(ENC_4)-cos(ENC_3)*sin(ENC_2)*sin(ENC_4);
	angle=acos(Z4[2]/sqrt((Z4[0]*Z4[0])+(Z4[1]*Z4[1])+(Z4[2]*Z4[2])));
	//std::cout<<"Angle Z: "<<angle<<"\n\n\n\n\n\n\n\n";
return angle;
}

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) 
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	jp_type zeroMove;

	
	wam.gravityCompensate(true);
	
	bool going = true;
	//int n;
	//int val;
	//wam.moveTo(jp_type(0.0));

	jp_type J_arr;
	
	std::cout << "Entering to the loop.\n";
	//waitForEnter();

	long idnum0 = 0;
	long idnum1 = 1;
	long overVoltage = 0;  //Returns >0 if overvoltage is detected
	float Sg0 = 0.0;  //Returns the voltage reading
	float Sg1 = 0.0;  //Returns the voltage reading
	float Sg2 = 0.0;  //Returns the voltage reading
	float Sg3 = 0.0;  //Returns the voltage reading
	float Sg4 = 0.0;  //Returns the voltage reading
	float Sg5 = 0.0;  //Returns the voltage reading
	float Fx = 0.0;  
	float Fy = 0.0; 
	float Fz = 0.0;  
	float Tx = 0.0;  
	float Ty = 0.0;  
	float Tz = 0.0;
	float angle_x = 0.0;
	float angle_y = 0.0;
	float angle_z = 0.0;
	J_arr[0]=0;
	J_arr[1]=0;//-3.1415/2;
	J_arr[2]=0;//3.1415/2;
	J_arr[3]=3.1415/2;
	wam.moveTo(J_arr);
	while (going)
	{
		EAnalogIn(&idnum1, 0, 0, 0, &overVoltage, &Sg0);
		EAnalogIn(&idnum1, 0, 2, 0, &overVoltage, &Sg1);
		EAnalogIn(&idnum1, 0, 4, 0, &overVoltage, &Sg2);
		EAnalogIn(&idnum0, 0, 0, 0, &overVoltage, &Sg3);
		EAnalogIn(&idnum0, 0, 2, 0, &overVoltage, &Sg4);
		EAnalogIn(&idnum0, 0, 4, 0, &overVoltage, &Sg5);
				
////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// BIAS ///////////////////////////////////
	J_arr=wam.getJointPositions();
	angle_x=Angle_X(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);
	angle_y=Angle_Y(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);
	angle_z=Angle_Z(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);		
	
	Sg0=Sg0-0.0;Sg1=Sg1-0.0;Sg2=Sg2+0.0;
	Sg3=Sg3-0.0;Sg4=Sg4+0.0;Sg5=Sg5-0.0;
	
////////////////////////////////////////////////////////////////////////////
		//std::cout<<"Sg0: "<<Sg0*1000<<"\n\n\n\nSg1: "<<Sg1*1000<<"\n\n\n\n\nSg2: "<<Sg2*1000<<"\n\n\n\n\nSg3: "<<Sg3*1000<<"\n\n\n\n\nSg4: "<<Sg4*1000<<"\n\n\n\n\nSg5: "<<Sg5*1000<<"\n\n\n\n\n";

		Fx=(-1.81798*Sg0)+(0.13728*Sg1)+(5.00500*Sg2)+(-42.29696*Sg3)+(-0.37203*Sg4)+(47.34741*Sg5);
		Fy=(-1.21862*Sg0)+(56.55767*Sg1)+(-0.53762*Sg2)+(-24.78708*Sg3)+(3.93253*Sg4)+(-27.44557*Sg5);
		Fz=(67.24390*Sg0)+(1.30897*Sg1)+(66.99002*Sg2)+(1.32217*Sg3)+(66.64551*Sg4)+(1.39404*Sg5);
		Tx=(-0.01538*Sg0)+(0.39768*Sg1)+(-1.06401*Sg2)+(-0.19800*Sg3)+(1.07781*Sg4)+(-0.16926*Sg5);
		Ty=(1.28994*Sg0)+(0.02545*Sg1)+(-0.67087*Sg2)+(0.28289*Sg3)+(-0.62034*Sg4)+(-0.34534*Sg5);
		Tz=(0.01601*Sg0)+(-0.70847*Sg1)+(0.08505*Sg2)+(-0.61752*Sg3)+(-0.00010*Sg4)+(-0.70166*Sg5);
		Fx=Fx+5.25+(2.825*cos(angle_x));
		Fy=Fy+0.2+(3*cos(angle_y));
		Fz=Fz+1.5+(3*cos(angle_z));
		Tx=Tx+0.01-(0.18*cos(angle_y));
		Ty=Ty-0.28+(0.19*cos(angle_x));
		Tz=Tz+0.225;
		std::cout<<"Fx: "<<Fx<<"\n\n\n\nFy: "<<Fy<<"\n\n\n\n\nFz: "<<Fz<<"\n\n\n\n\nTx: "<<Tx<<"\n\n\n\n\nTy: "<<Ty<<"\n\n\n\n\nTz: "<<Tz<<"\n\n\n\n\n";
	}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}
