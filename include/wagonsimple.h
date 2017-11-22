/*
 * wagonsimple.h
 *
 *  Created on: Nov 7, 2013
 *      Author: leonardo
 */

#ifndef WAGONSIMPLE_H_
#define WAGONSIMPLE_H_

#include "wagongroup.h"
#include "mbsim/objects/rigid_body.h"
#include "mbsim/frames/frame.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "mbsim/functions/kinematics/kinematics.h"

#include "openmbvcppinterface/cuboid.h"

using namespace MBSim;

class WagonSimple : public WagonGroup
{
public:
  WagonSimple(const std::string& name);

  /**
   * Overloaded function to set the total mass.
   * In this case, the whole mass is attributed to the box
   * @param totalMass
   */
  void
  setTotalMass(double totalMass);

  /**
   * Sets the geometryReferenceFrame relative position
   * @param x_, y_, z_ position components on a Cartesian frame
   * @param pos_ position vector on a Cartesian frame
   */
  void
  setGeometryReferenceFramePosition(double x_, double y_, double z_);
  void
  setGeometryReferenceFramePosition(fmatvec::Vec3 pos_);

  /**
   * Sets the frontBolsterConnection relative position
   * @param x_, y_, z_ position components on a Cartesian frame
   * @param pos_ position vector on a Cartesian frame
   */
  void
  setFrontBolsterConnectionPosition(double x_, double y_, double z_);
  void
  setFrontBolsterConnectionPosition(fmatvec::Vec3 pos_);

  /**
   * Sets the rearBolsterConnection relative position
   * @param x_, y_, z_ position components on a Cartesian frame
   * @param pos_ position vector on a Cartesian frame
   */
  void
  setRearBolsterConnectionPosition(double x_, double y_, double z_);
  void
  setRearBolsterConnectionPosition(fmatvec::Vec3 pos_);

  /**
   * Enables OpenMBV representation
   * @param enable Boolean flag
   */
  void enableOpenMBV(bool enable);

  double
  getHeight() const
  {
    return height;
  }

  void
  setHeight(double height)
  {
    this->height = height;
  }

  double
  getLength() const
  {
    return length;
  }

  void
  setLength(double length)
  {
    this->length = length;
  }

  double
  getWidth() const
  {
    return width;
  }

  void
  setWidth(double width)
  {
    this->width = width;
  }

  MBSim::RigidBody*
  getWagonBox()
  {
    return wagonbox;
  }

protected:
  MBSim::RigidBody* wagonbox;
  MBSim::FixedRelativeFrame* geometryReferenceFrame;
  MBSim::FixedRelativeFrame* frontBolsterConnection;
  MBSim::FixedRelativeFrame* rearBolsterConnection;
  /**
   * Wagon box dimensions
   */
  double length;
  double height;
  double width;
};

#endif /* WAGONSIMPLE_H_ */
