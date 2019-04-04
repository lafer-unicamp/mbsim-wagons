/* 
 * File:   main.cpp
 * Author: lbb
 *
 * Created on March 30, 2013, 2:20 PM
 */

#ifndef _OPENMP
#define _OPENMP
#endif

#include <cstdlib>
#include <stdio.h> /* needed this include just because Eclipse can't
                       track down namespace std from cstdlib */
#include <iostream>
#include <string>
#include <sstream>
#include "mbsim/integrators/integrators.h"
#include "inputTools.h"
#include "system.h"
#include "profilecontour.h"


using namespace std;
using namespace MBSim;

/* 
 * 
 */
int main(int argc, char** argv)
{

	// read input file
	vector<string> inputFileNames;
	if (argc >= 2){
		for (int i=1; i<argc; i++) inputFileNames.push_back(argv[i]);
	}
	else inputFileNames.push_back("input2.in");

	// settings
	double endTime = atof(searchParameter(inputFileNames[0],"SIMULATION_TIME").c_str());
	double timeStep = atof(searchParameter(inputFileNames[0],"TIME_STEP").c_str());
	double outputFrequencyInHertz =
			atof(searchParameter(inputFileNames[0],"OUTPUT_FREQHZ").c_str());

	// build single modules
	DynamicSystemSolver *sys1 = new System(searchParameter(inputFileNames[0],"SYSTEM_NAME"),inputFileNames[0]);

	// add modules to overall dynamic system
	sys1->setConstraintSolver(DynamicSystemSolver::FixedPointSingle);
	sys1->setImpactSolver(DynamicSystemSolver::FixedPointSingle);
	sys1->setFlushEvery(250);
	sys1->initialize();
	sys1->setStopIfNoConvergence(false,true);
	sys1->setMaxIter(1e5);
	sys1->setGeneralizedImpulseTolerance(1e-6);
	sys1->setGeneralizedForceTolerance(1e-6);
	sys1->setGeneralizedRelativePositionTolerance(5e-6);
	sys1->setGeneralizedRelativeVelocityTolerance(8e-2);

	MBSimIntegrator::TimeSteppingSSCIntegrator integrator;
	integrator.setRelativeTolerance(1e-7);
	integrator.setAbsoluteTolerance(1e-6);
	integrator.setInitialStepSize(timeStep);
	integrator.setEndTime(endTime);
	integrator.setPlotStepSize(1/outputFrequencyInHertz);
	// integrator.setStepSize(timeStep);
	integrator.setStepSizeMax(5e-3);
	integrator.setStepSizeMin(1e-10);
	integrator.integrate(*sys1);
	cout << "finished" << endl;
	cout << timeStep << endl;

	return 0;
}
