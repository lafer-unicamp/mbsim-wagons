/*
 Wagon with moving mass
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

#ifndef WAGONMOVMASS_H
#define WAGONMOVMASS_H

#include "wagonbox.h"
#include "wagonsimple.h"
#include "mbsim/links/spring_damper.h"
#include "mbsim/links/joint.h"
#include "mbsim/constitutive_laws/constitutive_laws.h"
#include "mbsim/functions/kinetics/linear_spring_damper_force.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/coilspring.h>
#endif

class WagonMovMass : public WagonSimple
{
public:
  /// Wagon with moving mass
  /// this implementation supposes there is a factor that splits the overall
  /// mass of the wagon between the free mass and the overall wagon body
  /// \brief implements a WagonBox with a moving mass on its center
  WagonMovMass(const std::string& name, double sFactor_);

  /// Overload setter for mass
  /// \param mass
  void
  setTotalMass(double m_);

  double
  getSplitFactor() const
  {
    return splitFactor;
  }

  void
  setSplitFactor(double splitFactor)
  {
    this->splitFactor = splitFactor;
  }

  /// Getters for moving mass, spring, and joint
  MBSim::RigidBody*
  getMovingMass()
  {
    return movingMass;
  }
  ;
  MBSim::SpringDamper*
  getMovingMassSpring()
  {
    return spring;
  }
  ;
  MBSim::Joint*
  getTransJoint()
  {
    return translational;
  }
  ;

private:
  /**
   * Percentage of the total mass that is assigned to the wagon box
   */
  double splitFactor;

  // mass attached to the wagons center of gravity
  MBSim::RigidBody* movingMass;
  // spring that connects the moving mass and the wagon
  MBSim::SpringDamper* spring;
  // joint to constrain the movement of the free mass to the wagon
  MBSim::Joint* translational;
};

#endif // WAGONMOVMASS_H
