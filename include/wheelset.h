/*
 * wheelset.h
 *
 *  Created on: Nov 27, 2013
 *      Author: leonardo
 */

#ifndef WHEELSET_H_
#define WHEELSET_H_

#include <string>
#include <mbsim/objects/rigid_body.h>
#include "mbsim/frames/fixed_relative_frame.h"
#include "openmbvcppinterface/rotation.h"
#include "openmbvcppinterface/compoundrigidbody.h"

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

protected:
	double track;
	double sideframeTrack;
};

#endif /* WHEELSET_H_ */
