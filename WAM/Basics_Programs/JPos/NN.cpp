// Lock joints together
#include <iostream>
#include <string>
#include <cstdlib>  // For std::atexit()
#include <algorithm>

#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/standard_main_function.h>

using namespace barrett;

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
	wam.gravityCompensate(true);
	
	bool going = true;
	int val;
	wam.moveTo(jp_type(0.0));

	jp_type J_arr;

	while (going)
	{
		std::cout <<">>> ";
		std::cin >> val;
		switch (val) {
		case 0:
			going = false;
			break;
		case 8:
			std::cout << "J1 = ";
			std::cin >> J_arr[0];
			std::cout << "J2 = ";
			std::cin >> J_arr[1];
			std::cout << "J3 = ";
			std::cin >> J_arr[2];
			std::cout << "J4 = ";
			std::cin >> J_arr[3];
			wam.moveTo(J_arr);
			break;
					}
	}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}