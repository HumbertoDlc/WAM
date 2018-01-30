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

//#include <stdlib.h>
//#define BARRETT_SMF_VALIDATE_ARGS
//#include<ctime>
//#include <time.h>

using namespace barrett;
using namespace std;
using systems::connect;

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
class Set_Torque : public systems::System
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:
	Input<jp_type> jp_in;
	Input<jv_type> jv_in;
	Output<jt_type> jt_out;
	Set_Torque(double *Torques, 
		const std::string& sysName = "Set_Torque") :
		System(sysName), 
		jp_in(this),
		jv_in(this),
		jt_out(this, &jtOutputValue),
		Tor(Torques),
		jt(0.0) {}
	virtual ~Set_Torque() { this->mandatoryCleanUp(); }

protected:
	typename Output<jt_type>::Value* jtOutputValue;
	double *Tor;
	jt_type jt;
	virtual void operate() 
	{
		jt[0]=0;//Tor[0];
		jt[1]=Tor[1];
		jt[2]=0;//Tor[2];
		jt[3]=0;//Tor[3];		
		this->jtOutputValue->setData(&jt);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(Set_Torque);
};

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) 
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	
	wam.gravityCompensate(true);
	
	jp_type J_arr;
	jv_type Jv_arr;
	
	// TXT related to the desired position, velocities and acceleration:
	int txt_long=1598;
	// LabJack IDs:
	long idnum0 = 0;long idnum1 = 1;
	long overVoltage = 0;
	// FT signals:
	float Sg0 = 0.0;float Sg1 = 0.0;float Sg2 = 0.0;
	float Sg3 = 0.0;float Sg4 = 0.0;float Sg5 = 0.0;
	// FT values:
	float Fx = 0.0;float Fy = 0.0;float Fz = 0.0;  
	float Tx = 0.0;float Ty = 0.0;float Tz = 0.0;
	// Orientation:	
	float angle_x = 0.0;float angle_y = 0.0;float angle_z = 0.0;
	// New position:	
	J_arr[0]=0;J_arr[1]=0;J_arr[2]=0;J_arr[3]=0;wam.moveTo(J_arr);

	//jp_type zeroMove;
	jp_type jp_in;
	
	// Desired position, velocity and acceleration:
	double Des[txt_long][12];	
	double x_des[txt_long][4];
	double xdot_des[txt_long][4];
	double xddot_des[txt_long][4];
	
	ifstream file;
	file.open("Reference.txt");
	for(int i=0;i<txt_long;++i){
	for(int j=0;j<12;++j){
	file>>Des[i][j];}}
	file.close();
	
	for(int i=0;i<txt_long;++i){
	for(int j=0;j<4;++j){
	x_des[i][j]=Des[i][j];
	xdot_des[i][j]=0;//Des[i][j+4];
	xddot_des[i][j]=0;}}//Des[i][j+8];}}
	
	// Desired impedance:
	double M_des[4]={5,5,2,2};					//{5,5,2,2};
	double K_des[4]={1600,60,1000,90};				//{1600,50,1000,90};
	double B_des[4]={100,30,80,50};					//{100,30,80,50};
	double Kf_des=1;
	
	// Estimated and SM parameters:
	double M[4][4];
	double C[4][4]={0};
	double G_v[4];
	double F1=1;
	double F2=1;
	double A=-10;
	double D=1000;
	double ep=0.001;

	// Gains:
	double K_vz[4],K_pz[4],K_fz[4];
	for(int i=0;i<4;i++){	
	K_vz[i]=(B_des[i]/M_des[i]-F1+A)/F2;
	K_pz[i]=(K_des[i]/M_des[i]+F1*A)/F2;
	K_fz[i]=(Kf_des/M_des[i])/F2;}
	
	// Variables declaration:
	double q1, q2, q3, q4;
	double J[4][6];
	double J_FT[4];
	double e[4];
	double edot[4];
	double z[4];
	double s[4];
	double sgn[4];
	double zdot[4];
	double Torques[4]={0};
	double T_1, T_2, T_3, T_4;
	double T_1_a, T_1_b, T_1_c, T_1_d, T_1_e, T_1_f, T_1_g, T_1_h;
	double T_2_a, T_2_b, T_2_c, T_2_d, T_2_e, T_2_f, T_2_g, T_2_h;
	double T_3_a, T_3_b, T_3_c, T_3_d, T_3_e, T_3_f, T_3_g, T_3_h;
	double T_4_a, T_4_b, T_4_c, T_4_d, T_4_e, T_4_f, T_4_g, T_4_h;

	std::string filestr;
	//const char *filename;
	
	printf("Please enter file name to log data to: ");
	std::getline(std::cin,filestr);
	//filename = filestr.c_str();
	
	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;}

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

	J_arr[0]=0;//x_des[0][0];
	J_arr[1]=x_des[0][1];
	J_arr[2]=0;//x_des[0][2];
	J_arr[3]=0;//x_des[0][3];
	wam.moveTo(J_arr);
	
	time.start();
	connect(tg.output, logger.input);
	printf("Logging started.\n");
	
	for (int i=0;i<txt_long;++i){

		Set_Torque<DOF> jsd(Torques);
		systems::connect(wam.jpOutput, jsd.jp_in);
		systems::connect(wam.jvOutput, jsd.jv_in);
		wam.trackReferenceSignal(jsd.jt_out);
		//usleep(500);

		// FT signals:		
		EAnalogIn(&idnum1, 0, 0, 0, &overVoltage, &Sg0);
		EAnalogIn(&idnum1, 0, 2, 0, &overVoltage, &Sg1);
		EAnalogIn(&idnum1, 0, 4, 0, &overVoltage, &Sg2);
		EAnalogIn(&idnum0, 0, 0, 0, &overVoltage, &Sg3);
		EAnalogIn(&idnum0, 0, 2, 0, &overVoltage, &Sg4);
		EAnalogIn(&idnum0, 0, 4, 0, &overVoltage, &Sg5);
				
		// Orientation:
		J_arr=wam.getJointPositions();
		q1=J_arr[0];q2=J_arr[1];q3=J_arr[2];q4=J_arr[3];
		Jv_arr=wam.getJointVelocities();
		angle_x=Angle_X(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);
		angle_y=Angle_Y(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);
		angle_z=Angle_Z(J_arr[0],J_arr[1],J_arr[2],J_arr[3]);		
		
		// Bias:		
		Sg0=Sg0-0.0;Sg1=Sg1-0.0;Sg2=Sg2+0.0;
		Sg3=Sg3-0.0;Sg4=Sg4+0.0;Sg5=Sg5-0.0;
		
		// FT values:
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
		//std::cout<<"Fx: "<<Fx<<"\n\n\n\nFy: "<<Fy<<"\n\n\n\n\nFz: "<<Fz<<"\n\n\n\n\nTx: "<<Tx<<"\n\n\n\n\nTy: "<<Ty<<"\n\n\n\n\nTz: "<<Tz<<"\n\n\n\n\n";

		// Jacobian:
		J[0][0]=(9*cos(q4)*(cos(q1)*sin(q3)+cos(q2)*cos(q3)*sin(q1)))/200-(11*sin(q1)*sin(q2))/20-(9*cos(q1)*sin(q3))/200-(sin(q4)*(cos(q1)*sin(q3)+
		cos(q2)*cos(q3)*sin(q1)))/2-(9*sin(q1)*sin(q2)*sin(q4))/200-(9*cos(q2)*cos(q3)*sin(q1))/200-(cos(q4)*sin(q1)*sin(q2))/2;
		J[1][0]=cos(q1)*((11*cos(q2))/20 + (cos(q2)*cos(q4))/2-(9*cos(q3)*sin(q2))/200+(9*cos(q2)*sin(q4))/200+(9*cos(q3)*cos(q4)*sin(q2))/200-
		(cos(q3)*sin(q2)*sin(q4))/2);
		J[2][0]=-((cos(q3)*sin(q1)+cos(q1)*cos(q2)*sin(q3))*(100*sin(q4)-9*cos(q4)+9))/200;
		J[3][0]=(cos(q1)*cos(q3)-cos(q2)*sin(q1)*sin(q3))*((cos(q2)*cos(q4))/2+(9*cos(q2)*sin(q4))/200+(9*cos(q3)*cos(q4)*sin(q2))/200-
		(cos(q3)*sin(q2)*sin(q4))/2)-sin(q2)*sin(q3)*((sin(q4)*(cos(q1)*sin(q3)+cos(q2)*cos(q3)*sin(q1)))/2-(9*cos(q4)*(cos(q1)*sin(q3)+
		cos(q2)*cos(q3)*sin(q1)))/200+(9*sin(q1)*sin(q2)*sin(q4))/200+(cos(q4)*sin(q1)*sin(q2))/2);
		//
		J[0][1]=(11*cos(q1)*sin(q2))/20-(9*sin(q1)*sin(q3))/200+(9*cos(q4)*(sin(q1)*sin(q3)-cos(q1)*cos(q2)*cos(q3)))/200-(sin(q4)*(sin(q1)*sin(q3)-
		cos(q1)*cos(q2)*cos(q3)))/2+(9*cos(q1)*cos(q2)*cos(q3))/200+(cos(q1)*cos(q4)*sin(q2))/2+(9*cos(q1)*sin(q2)*sin(q4))/200;
		J[1][1]=sin(q1)*((11*cos(q2))/20+(cos(q2)*cos(q4))/2-(9*cos(q3)*sin(q2))/200+(9*cos(q2)*sin(q4))/200+(9*cos(q3)*cos(q4)*sin(q2))/200-(cos(q3)*sin(q2)*sin(q4))/2);
		J[2][1]=((cos(q1)*cos(q3)-cos(q2)*sin(q1)*sin(q3))*(100*sin(q4)-9*cos(q4)+9))/200;
		J[3][1]=(cos(q3)*sin(q1)+cos(q1)*cos(q2)*sin(q3))*((cos(q2)*cos(q4))/2+(9*cos(q2)*sin(q4))/200+(9*cos(q3)*cos(q4)*sin(q2))/200-
		(cos(q3)*sin(q2)*sin(q4))/2)+sin(q2)*sin(q3)*((9*cos(q4)*(sin(q1)*sin(q3)-cos(q1)*cos(q2)*cos(q3)))/200-(sin(q4)*(sin(q1)*sin(q3)-
		cos(q1)*cos(q2)*cos(q3)))/2+(cos(q1)*cos(q4)*sin(q2))/2+(9*cos(q1)*sin(q2)*sin(q4))/200);
		//
		J[0][2]=0;
		J[1][2]=(9*cos(q2)*cos(q3)*cos(q4))/200-(9*cos(q2)*cos(q3))/200-(cos(q4)*sin(q2))/2-(9*sin(q2)*sin(q4))/200-(11*sin(q2))/20-(cos(q2)*cos(q3)*sin(q4))/2;
		J[2][2]=(sin(q2)*sin(q3)*(100*sin(q4)-9*cos(q4) +9))/200;
		J[3][2]=(9*cos(q2)*cos(q4))/200-(cos(q2)*sin(q4))/2-(cos(q3)*cos(q4)*sin(q2))/2-(9*cos(q3)*sin(q2)*sin(q4))/200;
		//
		J[0][3]=0;                                                      
		J[1][3]=-sin(q1);
		J[2][3]=cos(q1)*sin(q2);
		J[3][3]=-cos(q3)*sin(q1)-cos(q1)*cos(q2)*sin(q3);
		//
		J[0][4]=0;
		J[1][4]=cos(q1);
		J[2][4]=sin(q1)*sin(q2);
		J[3][4]=cos(q1)*cos(q3)-cos(q2)*sin(q1)*sin(q3);
		//
		J[0][5]=1;
		J[1][5]=0;
		J[2][5]=cos(q2);
		J[3][5]=sin(q2)*sin(q3);
		
		//Mass Matrix
		M[0][0]=cos(q2)*(3.2e-3*cos(q2) - 2.5e-3*cos(q3)*sin(q2) + 1.8e-5*sin(q2)*sin(q3)) + (3.3e-5*cos(q3)*sin(q1) - 0.34*cos(q1)*sin(q2) + 6.7e-3*sin(q1)*sin(q3) - 			6.7e-3*cos(q1)*cos(q2)*cos(q3) + 3.3e-5*cos(q1)*cos(q2)*sin(q3))*(6.0e-5*cos(q3)*sin(q1) - 0.62*cos(q1)*sin(q2) + 0.012*sin(q1)*sin(q3) - 
		0.012*cos(q1)*cos(q2)*cos(q3) + 6.0e-5*cos(q1)*cos(q2)*sin(q3)) + (3.3e-5*cos(q1)*cos(q3) + 6.7e-3*cos(q1)*sin(q3) + 0.34*sin(q1)*sin(q2) +
		6.7e-3*cos(q2)*cos(q3)*sin(q1) - 3.3e-5*cos(q2)*sin(q1)*sin(q3))*(6.0e-5*cos(q1)*cos(q3) + 0.012*cos(q1)*sin(q3) + 0.62*sin(q1)*sin(q2) + 
		0.012*cos(q2)*cos(q3)*sin(q1) - 6.0e-5*cos(q2)*sin(q1)*sin(q3)) + (cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2))*(1.5e-3*cos(q2)*cos(q4) + 
		0.015*cos(q2)*sin(q4) - 1.7e-5*sin(q2)*sin(q3) + 0.015*cos(q3)*cos(q4)*sin(q2) - 1.5e-3*cos(q3)*sin(q2)*sin(q4)) + (0.11*cos(q1)*sin(q3) -
		5.5e-4*cos(q1)*cos(q3) + 1.3*sin(q1)*sin(q2) - 0.096*cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) + 0.32*sin(q4)*(cos(q1)*sin(q3) + 
		cos(q2)*cos(q3)*sin(q1)) + 0.096*sin(q1)*sin(q2)*sin(q4) + 0.11*cos(q2)*cos(q3)*sin(q1) + 5.5e-4*cos(q2)*sin(q1)*sin(q3) + 
		0.32*cos(q4)*sin(q1)*sin(q2))*(0.045*cos(q1)*sin(q3) - 2.3e-4*cos(q1)*cos(q3) + 0.55*sin(q1)*sin(q2) - 0.04*cos(q4)*(cos(q1)*sin(q3) + 
		cos(q2)*cos(q3)*sin(q1)) + 0.13*sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) + 0.04*sin(q1)*sin(q2)*sin(q4) + 0.045*cos(q2)*cos(q3)*sin(q1) +
		2.3e-4*cos(q2)*sin(q1)*sin(q3) + 0.13*cos(q4)*sin(q1)*sin(q2)) + (0.12*sin(q1) + 9.2e-3*cos(q1)*cos(q2) - 0.06*cos(q1)*sin(q2))*(0.031*sin(q1) + 
		2.4e-3*cos(q1)*cos(q2) - 0.015*cos(q1)*sin(q2)) + (cos(q2)*cos(q4) - 1.0*cos(q3)*sin(q2)*sin(q4))*(2.9e-3*cos(q2)*cos(q4) + 1.5e-3*cos(q2)*sin(q4) - 
		2.1e-5*sin(q2)*sin(q3) + 1.5e-3*cos(q3)*cos(q4)*sin(q2) - 2.9e-3*cos(q3)*sin(q2)*sin(q4)) + (4.4e-3*cos(q1) - 6.6e-4*sin(q1))*(0.048*cos(q1) - 
		7.2e-3*sin(q1)) + (6.6e-4*cos(q1) + 4.4e-3*sin(q1))*(7.2e-3*cos(q1) + 0.048*sin(q1)) + (0.12*cos(q1) - 9.2e-3*cos(q2)*sin(q1) + 
		0.06*sin(q1)*sin(q2))*(0.031*cos(q1) - 2.4e-3*cos(q2)*sin(q1) + 0.015*sin(q1)*sin(q2)) - 3.0e-20*sin(q2)*(8.1e14*cos(q2) - 7.0e17*sin(q2)) + 
		cos(q2)*(0.016*cos(q2) - 2.5e-5*sin(q2)) + (0.55*cos(q1)*sin(q2) + 2.3e-4*cos(q3)*sin(q1) - 0.045*sin(q1)*sin(q3) + 0.04*cos(q4)*(sin(q1)*sin(q3) - 
		1.0*cos(q1)*cos(q2)*cos(q3)) - 0.13*sin(q4)*(sin(q1)*sin(q3) - 1.0*cos(q1)*cos(q2)*cos(q3)) + 0.045*cos(q1)*cos(q2)*cos(q3) + 
		2.3e-4*cos(q1)*cos(q2)*sin(q3) + 0.13*cos(q1)*cos(q4)*sin(q2) + 0.04*cos(q1)*sin(q2)*sin(q4))*(1.3*cos(q1)*sin(q2) + 5.5e-4*cos(q3)*sin(q1) - 
		0.11*sin(q1)*sin(q3) + 0.096*cos(q4)*(sin(q1)*sin(q3) - 1.0*cos(q1)*cos(q2)*cos(q3)) - 0.32*sin(q4)*(sin(q1)*sin(q3) - 1.0*cos(q1)*cos(q2)*cos(q3)) + 
		0.11*cos(q1)*cos(q2)*cos(q3) + 5.5e-4*cos(q1)*cos(q2)*sin(q3) + 0.32*cos(q1)*cos(q4)*sin(q2) + 0.096*cos(q1)*sin(q2)*sin(q4)) - 
		1.0*cos(q3)*sin(q2)*(2.5e-3*cos(q2) - 0.059*cos(q3)*sin(q2) + 7.4e-6*sin(q2)*sin(q3)) + sin(q2)*sin(q3)*(1.8e-5*cos(q2) - 7.4e-6*cos(q3)*sin(q2) + 
		0.059*sin(q2)*sin(q3)) - 1.0*sin(q2)*sin(q3)*(2.1e-5*cos(q2)*cos(q4) + 1.7e-5*cos(q2)*sin(q4) - 0.015*sin(q2)*sin(q3) + 1.7e-5*cos(q3)*cos(q4)*sin(q2) -
		2.1e-5*cos(q3)*sin(q2)*sin(q4)) + 0.11;
		M[0][1]=3.0e-4*cos(q2)*cos(q3) - 5.3e-4*sin(q2) - 3.7e-3*cos(q2) - 0.072*cos(q2)*sin(q3) - 4.6e-6*cos(q4)*sin(q2) + 5.2e-5*sin(q2)*sin(q4) - 
		6.4e-5*cos(q3)*cos(q3)*sin(q2) + 0.022*cos(q2)*cos(q4)*cos(q4)*sin(q3) + 9.2e-6*cos(q3)*cos(q3)*cos(q4)*sin(q2) - 1.0e-4*cos(q3)*cos(q3)*sin(q2)*sin(q4) + 
		5.2e-5*cos(q2)*cos(q3)*cos(q4) + 4.6e-6*cos(q2)*cos(q3)*sin(q4) + 0.038*cos(q2)*cos(q4)*sin(q3) + 0.059*cos(q3)*sin(q2)*sin(q3) - 
		0.18*cos(q2)*sin(q3)*sin(q4) - 0.05*cos(q3)*cos(q4)*cos(q4)*sin(q2)*sin(q3) - 8.6e-3*cos(q3)*cos(q4)*sin(q2)*sin(q3) - 0.05*cos(q2)*cos(q4)*sin(q3)*sin(q4) +
		0.029*cos(q3)*sin(q2)*sin(q3)*sin(q4) - 0.022*cos(q3)*cos(q4)*sin(q2)*sin(q3)*sin(q4);
		M[0][2]=0.065*cos(q2) - 8.6e-3*cos(q2)*cos(q4) + 0.072*cos(q3)*sin(q2) + 0.029*cos(q2)*sin(q4) + 3.0e-4*sin(q2)*sin(q3) - 0.05*cos(q2)*cos(q4)*cos(q4) +
		4.6e-6*sin(q2)*sin(q3)*sin(q4) - 0.022*cos(q3)*cos(q4)*cos(q4)*sin(q2) - 0.038*cos(q3)*cos(q4)*sin(q2) - 0.022*cos(q2)*cos(q4)*sin(q4) + 
		0.18*cos(q3)*sin(q2)*sin(q4) + 5.2e-5*cos(q4)*sin(q2)*sin(q3) + 0.05*cos(q3)*cos(q4)*sin(q2)*sin(q4);
		M[0][3]=5.2e-5*cos(q2)*cos(q4) + 4.6e-6*cos(q2)*sin(q4) + 0.061*sin(q2)*sin(q3) + 0.067*sin(q2)*sin(q3)*sin(q4) + 4.6e-6*cos(q3)*cos(q4)*sin(q2) - 
		5.2e-5*cos(q3)*sin(q2)*sin(q4) + 0.17*cos(q4)*sin(q2)*sin(q3);

		M[1][0]=3.0e-4*cos(q2)*cos(q3) - 5.3e-4*sin(q2) - 3.7e-3*cos(q2) - 0.072*cos(q2)*sin(q3) - 4.6e-6*cos(q4)*sin(q2) + 5.2e-5*sin(q2)*sin(q4) - 
		6.4e-5*cos(q3)*cos(q3)*sin(q2) + 0.022*cos(q2)*cos(q4)*cos(q4)*sin(q3) + 9.2e-6*cos(q3)*cos(q3)*cos(q4)*sin(q2) - 1.0e-4*cos(q3)*cos(q3)*sin(q2)*sin(q4) + 
		5.2e-5*cos(q2)*cos(q3)*cos(q4) + 4.6e-6*cos(q2)*cos(q3)*sin(q4) + 0.038*cos(q2)*cos(q4)*sin(q3) + 0.059*cos(q3)*sin(q2)*sin(q3) - 
		0.18*cos(q2)*sin(q3)*sin(q4) - 0.05*cos(q3)*cos(q4)*cos(q4)*sin(q2)*sin(q3) - 8.6e-3*cos(q3)*cos(q4)*sin(q2)*sin(q3) - 0.05*cos(q2)*cos(q4)*sin(q3)*sin(q4) +
		0.029*cos(q3)*sin(q2)*sin(q3)*sin(q4) - 0.022*cos(q3)*cos(q4)*sin(q2)*sin(q3)*sin(q4);
		M[1][1]=0.017*cos(2.0*q3) + 0.013*cos(2.0*q4) + 3.2e-5*sin(2.0*q3) + 5.6e-3*sin(2.0*q4) + 0.35*cos(q4) + 0.12*sin(q4) - 0.013*cos(2.0*q3)*cos(2.0*q4) -
		5.6e-3*cos(2.0*q3)*sin(2.0*q4) - 4.3e-3*cos(2.0*q3)*cos(q4) + 0.014*cos(2.0*q3)*sin(q4) - 4.6e-6*sin(2.0*q3)*cos(q4) + 5.2e-5*sin(2.0*q3)*sin(q4) + 1.1;
		M[1][2]=3.0e-4*cos(q3) - 0.072*sin(q3) + 5.2e-5*cos(q3)*cos(q4) + 4.6e-6*cos(q3)*sin(q4) + 0.038*cos(q4)*sin(q3) - 0.18*sin(q3)*sin(q4) + 
		0.022*cos(q4)*cos(q4)*sin(q3) - 0.05*cos(q4)*sin(q3)*sin(q4);
		M[1][3]=0.061*cos(q3) + 0.17*cos(q3)*cos(q4) + 0.067*cos(q3)*sin(q4) - 4.6e-6*cos(q4)*sin(q3) + 5.2e-5*sin(q3)*sin(q4);

		M[2][0]=0.04*cos(q2) - 8.6e-3*cos(q2)*cos(q4) + 0.061*cos(q3)*sin(q2) + 0.029*cos(q2)*sin(q4) + 3.0e-4*sin(q2)*sin(q3) - 
		0.025*cos(2.0*q4)*cos(q2) - 0.011*sin(2.0*q4)*cos(q2) + 4.6e-6*sin(q2)*sin(q3)*sin(q4) - 0.011*cos(2.0*q4)*cos(q3)*sin(q2) + 
		0.025*sin(2.0*q4)*cos(q3)*sin(q2) - 0.038*cos(q3)*cos(q4)*sin(q2) + 0.18*cos(q3)*sin(q2)*sin(q4) + 5.2e-5*cos(q4)*sin(q2)*sin(q3);
		M[2][1]=3.0e-4*cos(q3) - 0.072*sin(q3) + 5.2e-5*cos(q3)*cos(q4) + 4.6e-6*cos(q3)*sin(q4) + 0.038*cos(q4)*sin(q3) - 0.18*sin(q3)*sin(q4) + 
		0.022*cos(q4)*cos(q4)*sin(q3) - 0.05*cos(q4)*sin(q3)*sin(q4);
		M[2][2]=0.029*sin(q4) - 8.6e-3*cos(q4) - 0.022*cos(q4)*sin(q4) - 0.05*cos(q4)*cos(q4) + 0.065;
		M[2][3]=(0.04*sin(q1)*sin(q3)*sin(q4) - 0.04*cos(q1)*cos(q4)*sin(q2) + 0.13*cos(q1)*sin(q2)*sin(q4) + 0.13*cos(q4)*sin(q1)*sin(q3) - 
		0.13*cos(q1)*cos(q2)*cos(q3)*cos(q4) - 0.04*cos(q1)*cos(q2)*cos(q3)*sin(q4))*(0.11*cos(q3)*sin(q1) + 5.5e-4*sin(q1)*sin(q3) - 
		5.5e-4*cos(q1)*cos(q2)*cos(q3) + 0.11*cos(q1)*cos(q2)*sin(q3) - 0.096*cos(q3)*cos(q4)*sin(q1) + 0.32*cos(q3)*sin(q1)*sin(q4) - 
		0.096*cos(q1)*cos(q2)*cos(q4)*sin(q3) + 0.32*cos(q1)*cos(q2)*sin(q3)*sin(q4)) - 1.0*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3))*((cos(q3)*sin(q1) + 
		cos(q1)*cos(q2)*sin(q3))*(2.1e-5*cos(q4) + 1.7e-5*sin(q4)) + (2.9e-3*cos(q4) + 1.5e-3*sin(q4))*(cos(q1)*cos(q4)*sin(q2) - 1.0*sin(q1)*sin(q3)*sin(q4) + 
		cos(q1)*cos(q2)*cos(q3)*sin(q4)) + (1.5e-3*cos(q4) + 0.015*sin(q4))*(cos(q1)*sin(q2)*sin(q4) + cos(q4)*sin(q1)*sin(q3) - 1.0*cos(q1)*cos(q2)*cos(q3)*cos(q4))) -
		1.0*(cos(q1)*cos(q3) - 1.0*cos(q2)*sin(q1)*sin(q3))*((1.5e-3*cos(q4) + 0.015*sin(q4))*(cos(q1)*cos(q4)*sin(q3) - 1.0*sin(q1)*sin(q2)*sin(q4) + 
		cos(q2)*cos(q3)*cos(q4)*sin(q1)) - 1.0*(2.9e-3*cos(q4) + 1.5e-3*sin(q4))*(cos(q4)*sin(q1)*sin(q2) + cos(q1)*sin(q3)*sin(q4) + cos(q2)*cos(q3)*sin(q1)*sin(q4)) +
		(2.1e-5*cos(q4) + 1.7e-5*sin(q4))*(cos(q1)*cos(q3) - 1.0*cos(q2)*sin(q1)*sin(q3))) + (0.13*cos(q1)*cos(q4)*sin(q3) - 0.13*sin(q1)*sin(q2)*sin(q4) + 
		0.04*cos(q4)*sin(q1)*sin(q2) + 0.04*cos(q1)*sin(q3)*sin(q4) + 0.13*cos(q2)*cos(q3)*cos(q4)*sin(q1) + 
		0.04*cos(q2)*cos(q3)*sin(q1)*sin(q4))*(0.11*cos(q1)*cos(q3) + 5.5e-4*cos(q1)*sin(q3) - 0.096*cos(q1)*cos(q3)*cos(q4) + 5.5e-4*cos(q2)*cos(q3)*sin(q1) + 
		0.32*cos(q1)*cos(q3)*sin(q4) - 0.11*cos(q2)*sin(q1)*sin(q3) + 0.096*cos(q2)*cos(q4)*sin(q1)*sin(q3) - 0.32*cos(q2)*sin(q1)*sin(q3)*sin(q4)) + 
		2.6e-21*sin(q2)*(2.1e17*cos(q3) - 4.2e19*sin(q3) + 3.7e19*cos(q4)*sin(q3) - 1.2e20*sin(q3)*sin(q4))*(0.13*cos(q2)*sin(q4) - 0.04*cos(q2)*cos(q4) + 
		0.13*cos(q3)*cos(q4)*sin(q2) + 0.04*cos(q3)*sin(q2)*sin(q4)) + sin(q2)*sin(q3)*((cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2))*(1.5e-3*cos(q4) + 0.015*sin(q4)) + 
		(2.9e-3*cos(q4) + 1.5e-3*sin(q4))*(cos(q2)*cos(q4) - 1.0*cos(q3)*sin(q2)*sin(q4)) - 1.0*sin(q2)*sin(q3)*(2.1e-5*cos(q4) + 1.7e-5*sin(q4)));


		M[3][0]=5.2e-5*cos(q2)*cos(q4) + 4.6e-6*cos(q2)*sin(q4) + 0.061*sin(q2)*sin(q3) + 0.067*sin(q2)*sin(q3)*sin(q4) + 4.6e-6*cos(q3)*cos(q4)*sin(q2) - 
		5.2e-5*cos(q3)*sin(q2)*sin(q4) + 0.17*cos(q4)*sin(q2)*sin(q3);
		M[3][1]=0.061*cos(q3) + 0.17*cos(q3)*cos(q4) + 0.067*cos(q3)*sin(q4) - 4.6e-6*cos(q4)*sin(q3) + 5.2e-5*sin(q3)*sin(q4);
		M[3][2]=5.2e-5*cos(q4 - 0.089);
		M[3][3]=0.061;

		//Gravity Compensator
		G_v[0]=0;
		G_v[1]=0.09*cos(q2)-20.0*sin(q2)-1.2*cos(q2)*cos(q3)-4.8e-3*cos(q2)*sin(q3)-3.1*cos(q4)*sin(q2)-0.94*sin(q2)*sin(q4)+0.94*cos(q2)*cos(q3)*cos(q4)-
		3.1*cos(q2)*cos(q3)*sin(q4);
                G_v[2]=-9.6e-37*sin(q2)*(5.0e33*cos(q3)-1.2e36*sin(q3)+9.8e35*cos(q4)*sin(q3)-3.2e36*sin(q3)*sin(q4));
		G_v[3]=0.94*cos(q2)*cos(q4)-3.1*cos(q2)*sin(q4)-3.1*cos(q3)*cos(q4)*sin(q2)-0.94*cos(q3)*sin(q2)*sin(q4);

		for(int j = 0; j < 4; ++j){
		e[j]=J_arr[j]-x_des[i][j];
		edot[j]=Jv_arr[j]-xdot_des[i][j];
		J_FT[j] =(J[j][0]*Fx+J[j][1]*Fy+J[j][2]*Fz+J[j][3]*Tx+J[j][4]*Ty+J[j][5]*Tz);
		z[j]=-(edot[j]+F1*e[j])/F2;
		s[j]=edot[j]+F1*e[j]+F2*z[j];
		sgn[j]=2/(1+exp(-s[j]))-1;
		zdot[j]=(A*z[j]+K_pz[j]*e[j]+K_vz[j]*edot[j]+K_fz[j]*J_FT[j]);}

		T_1_a=0;T_1_b=0,T_1_d=0;T_1_e=0;
		T_2_a=0;T_2_b=0,T_2_d=0;T_2_e=0;
		T_3_a=0;T_3_b=0,T_3_d=0;T_3_e=0;
		T_4_a=0;T_4_b=0,T_4_d=0;T_4_e=0;				
		for(int j = 0; j < 4; ++j){
		T_1_a=T_1_a+M[0][j]*(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_1_b=T_1_b+C[0][j]*(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_1_d=T_1_d-M[0][j]*0.25*abs(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_1_e=T_1_e-C[0][j]*0.25*abs(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_2_a=T_2_a+M[1][j]*(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_2_b=T_2_b+C[1][j]*(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_2_d=T_2_d-M[1][j]*0.25*abs(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_2_e=T_2_e-C[1][j]*0.25*abs(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_3_a=T_3_a+M[2][j]*(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_3_b=T_3_b+C[2][j]*(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_3_d=T_3_d-M[2][j]*0.25*abs(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_3_e=T_3_e-C[2][j]*0.25*abs(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_4_a=T_4_a+M[3][j]*(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_4_b=T_4_b+C[3][j]*(xdot_des[i][j]-F1*e[j]-F2*z[j]);
		T_4_d=T_4_d-M[3][j]*0.25*abs(xddot_des[i][j]-F1*edot[j]-F2*zdot[j]);
		T_4_e=T_4_e-C[3][j]*0.25*abs(xdot_des[i][j]-F1*e[j]-F2*z[j]);}
		T_1_c=G_v[0];		
		T_1_f=-G_v[0]*0.25-J_FT[0]*0.1;
		T_1_g=-D*s[0]-ep*sgn[0];
		T_1_h=J_FT[0];
		T_2_c=G_v[1];		
		T_2_f=-G_v[1]*0.25-J_FT[1]*0.1;
		T_2_g=-D*s[1]-ep*sgn[1];
		T_2_h=J_FT[1];
		T_3_c=G_v[2];		
		T_3_f=-G_v[2]*0.25-J_FT[2]*0.1;
		T_3_g=-D*s[2]-ep*sgn[2];
		T_3_h=J_FT[2];
		T_4_c=G_v[3];		
		T_4_f=-G_v[3]*0.25-J_FT[3]*0.1;
		T_4_g=-D*s[3]-ep*sgn[3];
		T_4_h=J_FT[3];		
		T_1=T_1_a+T_1_b+T_1_c+(T_1_d+T_1_e+T_1_f)*sgn[0]+T_1_g+T_1_h;
		T_2=T_2_a+T_2_b+T_2_c+(T_2_d+T_2_e+T_2_f)*sgn[1]+T_2_g+T_2_h;
		T_3=T_3_a+T_3_b+T_3_c+(T_3_d+T_3_e+T_3_f)*sgn[2]+T_3_g+T_3_h;
		T_4=T_4_a+T_4_b+T_4_c+(T_4_d+T_4_e+T_4_f)*sgn[3]+T_4_g+T_4_h;

		//cout<<"T_1_a: "<<T_1_a<<"\nT_1_b: "<<T_1_b<<"\nT_1_c: "<<T_1_c<<"\nT_1_d: "<<T_1_d<<"\nT_1_e: "<<T_1_e<<"\nT_1_f: "<<T_1_f<<"\nT_1_g: "<<T_1_g<<"\nT_1_h: "<<T_1_h<<"\ni: "<<i<<"\n\n\n\n\n\n";
		//cout<<"T_2_a: "<<T_2_a<<"\nT_2_b: "<<T_2_b<<"\nT_2_c: "<<T_2_c<<"\nT_2_d: "<<T_2_d<<"\nT_2_e: "<<T_2_e<<"\nT_2_f: "<<T_2_f<<"\nT_2_g: "<<T_2_g<<"\nT_2_h: "<<T_2_h<<"\ni: "<<i<<"\n\n\n\n\n\n";
		cout<<"T_4_a: "<<T_4_a<<"\nT_4_b: "<<T_4_b<<"\nT_4_c: "<<T_4_c<<"\nT_4_d: "<<T_4_d<<"\nT_4_e: "<<T_4_e<<"\nT_4_f: "<<T_4_f<<"\nT_4_g: "<<T_4_g<<"\nT_4_h: "<<T_4_h<<"\ni: "<<i<<"\n\n\n\n\n\n";	
			
			Torques[0]=T_1;
			Torques[1]=T_2;
			Torques[2]=T_3;
			Torques[3]=T_4;
}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}
