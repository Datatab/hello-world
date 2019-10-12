/*
 * Solver.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Method.h"

class Solver {
public:
	static Method* initMethod( char* fileName );
	static void runMethod( Method* m );
	static void destroyMethod( Method* m );
};

#endif /* SOLVER_H_ */
