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
	float pos=0.0;
	wam.moveTo(jp_type(pos));

	cp_type cartSetPoint;

	while (going)
	{
		std::cout <<pos<<">>> ";
		std::cin >> val;
		switch (val) {
		case 0:
			going = false;
			break;
		case 8:
			std::cout << "x = ";
			std::cin >> cartSetPoint[0];
			std::cout << "y = ";
			std::cin >> cartSetPoint[1];
			std::cout << "z = ";
			std::cin >> cartSetPoint[2];
			wam.moveTo(cartSetPoint);
			break;
					}
	}
	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}