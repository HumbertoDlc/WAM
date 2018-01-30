#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>  // For std::atexit()
#include <algorithm>

#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/standard_main_function.h>

using namespace barrett;
using namespace std;

void waitForEnter() 
{
	std::string line;
	std::getline(std::cin, line);
}

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) 
{
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);
	jp_type zeroMove;

	std::cout << "Press [Enter].\n";
	waitForEnter();
	wam.gravityCompensate(false);
	
	wam.moveTo(jp_type(0.0));

	jp_type J_arr;

	float theta[160][4];
	ifstream file;
	file.open("exp.txt");
	
	for(int i=0;i<160;++i){
	for(int j=0;j<4;++j){
	file>>theta[i][j];}}

	file.close();

	for(int k=0;k<160;++k){
	{
	J_arr[0]=theta[k][0];
	J_arr[1]=theta[k][1];
	J_arr[2]=theta[k][2];
	J_arr[3]=theta[k][3];
	wam.moveTo(J_arr);					}
	}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}