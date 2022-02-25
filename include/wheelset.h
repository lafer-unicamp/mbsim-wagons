/*
 * wheelset.h
 *
 *  Created on: Nov 27, 2013
 *      Author: leonardo
 */

#ifndef WHEELSET_H_
#define WHEELSET_H_

#include <string>
#include "mbsim/objects/rigid_body.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "openmbvcppinterface/rotation.h"
#include "openmbvcppinterface/frustum.h"
#include "openmbvcppinterface/compoundrigidbody.h"
#include "openmbvcppinterface/ivbody.h"
#include "wheel_profile.h"

class Wheelset : public MBSim::RigidBody
{
public:
	// Default constructor
	Wheelset(const std::string& name);
	Wheelset(const std::string& name, double _track, double _sfTrack);

	/*
	 * Enable graphical representation
	 * To describe the wheel profile, the OpenMBV::Rotation class is used
	 * This class performs a rotation around the Y axis
	 */
	void enableOpenMBV();

	/*
	 * GETTERS AND SETTERS
	 */
	MBSim::WheelProfile* getLeftWheel() { return wheelLeft; };
	MBSim::WheelProfile* getRightWheel() { return wheelRight; };

protected:
	double track;
	double sideframeTrack;
	MBSim::WheelProfile *wheelRight;
	MBSim::WheelProfile *wheelLeft;
};

#endif /* WHEELSET_H_ */
