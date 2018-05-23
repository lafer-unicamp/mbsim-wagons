/*
    Implements a truck with "variable" damping
    Copyright (C) 2013  Leonardo Baruffaldi leobart@fem.unicamp.br

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef BARBERTRUCK_H
#define BARBERTRUCK_H

#include "truck.h"
#include "wedge.h"
#include "bolster.h"
#include "sideframe.h"
#include "sinusoidalmovement.h"
#include "wheelset.h"
#include "jointFactory.h"
#include "mbsim/objects/rigid_body.h"
#include "mbsim/links/joint.h"
#include "mbsim/links/spring_damper.h"
#include "mbsim/links/elastic_joint.h"
#include "mbsim/environment.h"
#include "mbsim/contours/sphere.h"
#include "mbsim/contours/plane.h"
#include "mbsim/contours/edge.h"
#include "mbsim/contours/compound_contour.h"
#include "mbsim/links/contact.h"
#include "mbsim/functions/kinetics/kinetics.h"
#include "mbsim/constitutive_laws/constitutive_laws.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "mbsim/utils/rotarymatrices.h"
#include "mbsim/observers/contact_observer.h"
#include "mbsim/observers/mechanical_link_observer.h"
#include "mbsim/functions/kinetics/linear_elastic_function.h"
#include "mbsim/links/generalized_elastic_connection.h"

/*
 * This object is not included in the mbsim 11.0 official distribution,
 * so it was added individually
 */
#include "mbsim/links/isotropic_rotational_spring_damper.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include "openmbvcppinterface/coilspring.h"
#include "openmbvcppinterface/cuboid.h"
#include "openmbvcppinterface/extrusion.h"
#include "openmbvcppinterface/sphere.h"
#endif

#include <iostream>

using namespace MBSim;
using namespace fmatvec;
using namespace std;



class BarberTruck : public Truck
{
public:
	BarberTruck ( const std::string& projectName );
	BarberTruck ( const std::string& projectName, bool withBushings, double _wBase );

	Wedge* wedge1;
	Wedge* wedge2;
	Wedge* wedge3;
	Wedge* wedge4;
	Bolster* bolster;
	Sideframe* sideFrameLeft;
	Sideframe* sideFrameRight;
	Wheelset* wheelRear;
	Wheelset* wheelFront;

	void setWheelBase(double wB_);
	double getWheelBase();

private:
	double wheelBase;
	bool bolsterWithBushings;

	/// \brief Sets up the position of the sideframe springs connection points
	/// \param sframe, the sideframe rigid body to project the frames
	/// \param _distanceBetweenSprings distance between spring packs (should be the same along x and z)
	/// \param _lateralOffset lateral distance between the central spring packs of the two sideframes
	/// \param _designLength vertical distance between spring ends.
	/// TODO the parameter _distanceBetweenSprings sets a property from the bolster class
	/// called springOffsetDistace. Maybe it is better to configure it with a setter?
	void setSpringConnectionPointsFrames(MBSim::RigidBody *body, double _distanceBetweenSprings,
			double _lateralOffset, double _height, int lado);

	/// \brief configures wedge contacts
	/// \param Name of the wedge face TODO improve this
	/// \param Other component face contour
	/// \param Friction coefficient
	/// \param Restitution coefficient
	/// \param Whether there is an observer (default = false, no observer)
	void setWedgeContacts(Contour *,Contour *,double ,double, bool observerActive=false );

};

#endif // BARBERTRUCK_H
