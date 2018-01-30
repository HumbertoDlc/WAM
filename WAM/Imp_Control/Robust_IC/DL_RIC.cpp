//#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()
#include <iostream> 
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>  
#include <stdlib.h>

#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>

#include <curses.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/os.h>  // For btsleep()
#include <barrett/math.h>  // For barrett::math::saturate()

#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>

#include<ctime>
#include <time.h>

using namespace barrett;
using namespace std;
using systems::connect;

bool validate_args(int argc, char** argv) {
	return true;
}

void waitForEnter() 
{
	std::string line;
	std::getline(std::cin, line);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
		jt[0]=Tor[0];
		jt[1]=Tor[1];
		jt[2]=Tor[2];
		jt[3]=Tor[3];		
		this->jtOutputValue->setData(&jt);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(Set_Torque);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	
	int txt_long=6600;

	jp_type zeroMove;
	jp_type jp_in;

	std::cout << "Press [Enter].\n";
	waitForEnter();
	wam.gravityCompensate(true);
	int going = txt_long;
	double x_des[txt_long][4];

	//////////////////
	double xdot=0;
	double xdot_des=0;
	double xddot_des=0;
	//////////////////
	
	// Desired Impedance
	double M_des=5;
	double K_des=30;
	double B_des=1;
	double Kf_des=1;
	
	// Constants
	double M_est=1.6;
	double C_est=0.2;
	double G_v=0;
	double F_Delta=0;
	double F1=1;
	double F2=1;
	double A=-10;
	double D=500;
	double ep=0.001;

	// Gains
	double K_vz=(B_des/M_des-F1+A)/F2;
	double K_pz=(K_des/M_des+F1*A)/F2;
	double K_fz=(Kf_des/M_des)/F2;


	double e;
	double edot;
	double ef;
	double z;
	double s;
	double sgn;
	double zdot;


	double T_1;
	double T_2;
	double T_3;
	double T_4_1;
	double T_4_2;
	double T_4;
	double T_5;
	double T_6;
	double Torque_1;


	ifstream file;
	file.open("exp.txt");
	
	for(int i=0;i<txt_long;++i){
	for(int j=0;j<4;++j){
	file>>x_des[i][j];}}
	file.close();

	std::string filestr;
	const char *filename;
	
	printf("Please enter file name to log data to: ");
	std::getline(std::cin,filestr);
	filename = filestr.c_str();
	
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

	double Torques[] = {0, 0, 0, 0};
	int k=0;
	jp_type x;
	
	jp_type J_arr;
	J_arr[0]=x_des[0][0];
	wam.moveTo(J_arr);	
	
	while (going>k)
	{
	k=k+1;
		x=wam.getJointPositions();
		Set_Torque<DOF> jsd(Torques);
		systems::connect(wam.jpOutput, jsd.jp_in);
		systems::connect(wam.jvOutput, jsd.jv_in);
		wam.trackReferenceSignal(jsd.jt_out);
		usleep(1000);
		
	//Variables
	e=x[0]-x_des[k][0];
	edot=xdot-xdot_des;
	ef=J_FT;                    // F??
	z=-(edot+F1*e)/F2;
	s=edot+F1*e+F2*z;
	sgn=2/(1+exp(-s))-1;
	zdot=(A*z+K_pz*e+K_vz*edot+K_fz*ef);


	T_1=M_est*(xddot_des-F1*edot-F2*zdot);
	T_2=C_est*(xdot_des-F1*e-F2*z);
	T_3=G_v;
	T_4_1=(M_est*0.25)*abs(xddot_des-F1*edot-F2*zdot);
	T_4_2=(C_est*0.25)*abs(xdot_des-F1*e-F2*z);
	T_4=(T_4_1+T_4_2+(G_v*0.25)+(J_FT*0.1))*sgn;
	T_5=D*s+ep*sgn;
	T_6=J_FT;
	Torque_1=T_1+T_2+T_3-T_4-T_5+T_6;
	
	cout<<Torque_1;
			Torques[0]=Torque_1;
			Torques[1]=0;
			Torques[2]=0;
			Torques[3]=0;
					}

	wam.moveHome();
	
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
