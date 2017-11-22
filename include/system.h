/*
    System: instances a friction wedge 2D dynamical system
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


#ifndef _SYSTEM_H
#define _SYSTEM_H

#ifndef PI
#define PI 3.1416
#endif

#include "mbsim/dynamic_system.h"
#include "mbsim/objects/rigid_body.h"
#include "mbsim/links/spring_damper.h"
#include "mbsim/links/isotropic_rotational_spring_damper.h"
#include "mbsim/environment.h"
#include "mbsim/contours/sphere.h"
#include "mbsim/contours/plane.h"
#include "mbsim/contours/edge.h"
#include <mbsim/contours/compound_contour.h>
#include "mbsim/links/contact.h"
#include "mbsim/constitutive_laws/constitutive_laws.h"
#include "mbsim/functions/sinusoidal_function.h"

#include "openmbvcppinterface/coilspring.h"
#include "openmbvcppinterface/cuboid.h"
#include "openmbvcppinterface/extrusion.h"
#include "openmbvcppinterface/sphere.h"

#include <iostream>

using namespace MBSim;
using namespace fmatvec;
using namespace std;

#include <string>

#include "mbsim/dynamic_system_solver.h"
#include <barbertruck.h>
#include "wagonbox.h"
#include "wagonmovmass.h"
#include "wagonsimple.h"
#include "wagonsloshing.h"
#include "wheelset.h"
#include "sinusoidalmovement.h"
#include "inputTools.h"
/*
 * This object is not included in the mbsim 11.0 official distribution,
 * so it was added individually
 */
#include <isotropic_rotational_spring_damper.hx>

class System : public MBSim::DynamicSystemSolver
{
public:
  System(const std::string &projectName, const std::string &inputFileName);
  /// GETTERS AND SETTERS
  
  
private:
  // Sinusoidal movement inputs
  SinusoidalMovement *wheel1; // front wheel, front truck
  SinusoidalMovement *wheel2; // rear wheel, front truck
  SinusoidalMovement *wheel3; // front wheel, rear truck
  SinusoidalMovement *wheel4; // rear wheel, rear truck
  double amplitude;	// movement amplitude [m]
  double freq;	// movement angular speed [rad/s]
  double t0;	// movement delay to enter [s]
  double truckBaseDistance; // car wheel base [m]
  double truckWheelBase; // truck wheel base [m]
  double wagonMass; // wagon box mass [kg]
  double fillRatio; // amount of fluid in the box for the case of liquid cargo [0-1]
  fmatvec::SymMat wagonInertiaTensor; // [kg.m^2]
  bool bolsterBushing; // whether the bolster connections are modelled using 3d-stiffness
  
  /// \brief Get input data from text file
  /// \param inputFileName string containing the name of the input file
  /// \return 0 if everything is fine
  int initializeFromFile(const std::string &inputFileName);
};

#endif // _SYSTEM_H
