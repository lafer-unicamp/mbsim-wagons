/*
 * wheelBox.cpp
 *
 *  Created on: 11 de set de 2018
 *      Author: leonardo
 */

#include "wheelBox.h"

WheelBox::WheelBox(const std::string &name, Frame* sideFrameFrame,
		Frame* wheelSetFrame, SymMat stiffnessMatrix) : ElasticJoint(name) {

	this->stiffnessFunction = new LinearElasticFunction();

	this->setStiffnessMatrix(stiffnessMatrix);

	this->connect(sideFrameFrame,wheelSetFrame);
	this->setGeneralizedForceFunction(stiffnessFunction);
	this->setForceDirection("[1,0,0;0,1,0;0,0,1]");
	this->setMomentDirection("[1,0,0;0,1,0;0,0,1]");

	this->setPlotFeature(force,true);
	this->setPlotFeature(MBSim::moment,true);
}

void WheelBox::setStiffnessMatrix(SymMatV _stiffness){

	this->myStiffnessMatrix <<= _stiffness;

	this->stiffnessFunction->setStiffnessMatrix(this->myStiffnessMatrix);
	this->stiffnessFunction->setDampingMatrix(this->myStiffnessMatrix * 0.2);

}
