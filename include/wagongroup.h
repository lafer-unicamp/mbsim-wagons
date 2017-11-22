/*
 * wagongroup.h
 * This class implements the wagon as a group, instead of a RigidBody
 * allowing the addition of features and simplifying its interface with
 * the complete system
 *
 * The present class is meant to be an abstract class that should be
 * specialized
 *
 *  Created on: Nov 6, 2013
 *      Author: leonardo
 */

#ifndef WAGONGROUP_H
#define WAGONGROUP_H

#include "mbsim/group.h"

class WagonGroup : public MBSim::Group
{
public:
  WagonGroup(const std::string& name);

  void
  setInertiaTensor(fmatvec::SymMat I_)
  {
    inertiaTensor = I_;
  }

  fmatvec::SymMat
  getInertiaTensor() const
  {
    return inertiaTensor;
  }

  double
  getTotalMass() const
  {
    return totalMass;
  }

  void
  setTotalMass(double totalMass)
  {
    this->totalMass = totalMass;
  }

private:
  /**
   * Total mass of the wagon box
   */
  double totalMass;
  /**
   * Total inertia tensor of the wagon box
   */
  fmatvec::SymMat inertiaTensor;

};
#endif /* WAGONGROUP_H */
