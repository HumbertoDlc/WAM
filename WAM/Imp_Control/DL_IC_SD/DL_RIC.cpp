#include <iostream>
#include <string>
#include <cstdlib>  // For std::atexit()
#include <curses.h>
#include <cstdio>  // For remove()
#include <fstream>
#include <algorithm>
#include <math.h>  
#include <stdlib.h>

#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/os.h>  // For btsleep()
#include <barrett/math.h>  // For barrett::math::saturate()
#include <barrett/standard_main_function.h>

using namespace barrett;
using namespace std;
using detail::waitForEnter;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t DOF>
class JSpringDamper : public systems::System
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:
	Input<jp_type> jp_in;
	Input<jv_type> jv_in;
	Output<jt_type> jt_out;
	JSpringDamper(double *centerAngle,
		double *springConstant,
		double *dampingConstant, 
		const std::string& sysName = "JSpringDamper") :
		System(sysName), 
		jp_in(this),
		jv_in(this),
		jt_out(this, &jtOutputValue),
		center(centerAngle),
		stiffness(springConstant),
		damping(dampingConstant), 
		jt(0.0) {}
	virtual ~JSpringDamper() { this->mandatoryCleanUp(); }

protected:
	typename Output<jt_type>::Value* jtOutputValue;
	double *center;
	double *stiffness, *damping;
	jt_type jt;
	jp_type theta;
	jv_type thetaD; 

	virtual void operate() 
	{
		theta = jp_in.getValue();
		thetaD = jv_in.getValue();

		for (int i = 0; i < 4; i++)
		{
			jt[i] = (center[i] - theta[i]) * (*(stiffness+i)) - thetaD[i] * (*(damping+i));
		}

		this->jtOutputValue->setData(&jt);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(JSpringDamper);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) 
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	//wam.moveHome();
	wam.gravityCompensate(true);
	
	double stiffness[] = {800, 200, 10, 30};
	double stiffness_up_2 = 100;//800
	double stiffness_down_2 = 100;
	double stiffness_up_4 = 70;
	double stiffness_down_4 = 70;
	double damping[] = {20, 50, 5, 10};
	int txt_long=15991;
	double theta_des[txt_long][4];
	double theta_des_in[4];

	ifstream file;
	file.open("exp.txt");
	
	for(int i=0;i<txt_long;++i){
	for(int j=0;j<4;++j){
	file>>theta_des[i][j];}}
	file.close();

	jp_type J_arr;
	J_arr[0]=theta_des[0][0];
	J_arr[1]=theta_des[0][1];
	J_arr[2]=theta_des[0][2];
	J_arr[3]=theta_des[0][3]+0.6;
	wam.moveTo(J_arr);
	printf("Press [Enter] to add a virtual spring-damper to all joints.");
	waitForEnter();
	for(int times=1;times<4;times++){
	for(int k=1;k<txt_long+1;k++){
	theta_des_in[0]=theta_des[k][0];
	theta_des_in[1]=theta_des[k][1];
	theta_des_in[2]=theta_des[k][2];
	theta_des_in[3]=theta_des[0][3]+0.6;

	if (theta_des[k][1]<theta_des[k-1][1])
	{stiffness[1]=stiffness_down_2;}
	else{stiffness[1]=stiffness_up_2;}
	if (theta_des[k][3]<theta_des[k-1][3])
	{stiffness[3]=stiffness_down_4;}
	else{stiffness[3]=stiffness_up_4;}
	if (theta_des[k][3]<theta_des[0][3])
	{theta_des_in[3]=theta_des[0][3]+0.6;}

	//while (abs(theta_des_in[0]-wam.getJointPositions()[0])>0.01){
	JSpringDamper<DOF> jsd(theta_des_in,stiffness,damping);
	systems::connect(wam.jpOutput, jsd.jp_in);
	systems::connect(wam.jvOutput, jsd.jv_in);

	wam.trackReferenceSignal(jsd.jt_out);
	usleep(500);}}//}
	
	wam.moveHome();
	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}

