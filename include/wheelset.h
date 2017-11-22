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
#include "openmbvcppinterface/rotation.h"
#include "openmbvcppinterface/compoundrigidbody.h"

class Wheelset : public MBSim::RigidBody
{
public:
  Wheelset(const std::string& name);

  /*
   * Enable graphical representation
   * To describe the wheel profile, the OpenMBV::Rotation class is used
   * This class performs a rotation around the Y axis
   */
  void enableOpenMBV();

  double
  getTrack() const
  {
    return track;
  }

  void
  setTrack(double track)
  {
    this->track = track;
  }

protected:
  double track;
};

#endif /* WHEELSET_H_ */
