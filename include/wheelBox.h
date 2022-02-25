/*
 * wheelBox.h
 *
 *  Created on: 11 de set de 2018
 *      Author: leonardo
 *
 *      This class implements a generic wheelbox for freight wagons. Derived classes will hold the effective values for each type of wagon
 *		z axis is the axial direction
 */

#ifndef WHEELBOX_H_
#define WHEELBOX_H_

#include <mbsim/links/elastic_joint.h>
#include <mbsim/functions/kinetics/linear_elastic_function.h>

using namespace MBSim;
using namespace fmatvec;

class WheelBox: public ElasticJoint {
public:
	WheelBox(const std::string &name, Frame* sideFrameFrame, Frame* wheelSetFrame, SymMat stiffnessMatrix);
	void setStiffnessMatrix(SymMatV _stiffness);

private:
	LinearElasticFunction* stiffnessFunction;
	SymMatV myStiffnessMatrix;
};

#endif /* WHEELBOX_H_ */
