#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <cstdlib>

using namespace std;

class integrator {

	virtual double run_dynamics(double moves) = 0;

};

#endif
