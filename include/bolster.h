/*
    Bolster. Implements a three-piece-truck bolster
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


#ifndef BOLSTER_H
#define BOLSTER_H

#include <iostream>

#include <mbsim/objects/rigid_body.h>
#include <mbsim/contours/plane.h>
#include "mbsim/contours/plate.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "mbsim/frames/floating_relative_frame.h"
#include "mbsim/frames/contour_frame.h"
#include "mbsim/utils/rotarymatrices.h"
#include <fmatvec/fmatvec.h>

#include "openmbvcppinterface/extrusion.h"

using namespace MBSim;

class Bolster : public MBSim::RigidBody
{
public:
  /// Default constructor
  /// \brief default class constructor
  /// \param name Name of the object
  Bolster (const std::string& name);
  
  // Setters | Getters
  void setAngle(double ang_); /// \param ang_ angle of the bolster face
  void setHeight(double height_); /// \param height_ bolster height
  void setWidth(double width_); /// \param width_ bolster width at the lowest plane
  void setLength(double length_); /// \param length_ bolter length in transversal direction
  void setWagonConnectionPointPosition(fmatvec::Vec pos_); /// \param pos_ Vec 3x3 with position
  void setWagonConnectionPointPosition(double x_, double y_, double z_); /// \param x_ \param y_ \param z_ position coordinates referenced on the center of mass
  void setGeometryReferenceFramePosition(const fmatvec::Vec3 *pos_); /// \param pos_ Vec 3x3 with position referenced on the center of mass
  void setGeometryReferenceFramePosition(double x_, double y_, double z_); /// \param x_ \param y_ \param z_ position coordinates referenced on the center of mass
  /// This method should be called after all body's parameters had already been defined
  /// In case there is any changes in parameters, a new body should be created
  /// \brief adjusts the position of the contact plane based on geometric parameters
  void setContactPlanes();

  /// OpenMBV interface
  /// \brief enables OpenMBV visualization as an extrusion
  void enableOpenMBV(bool enable);
  
 private:
  // private parameters
  double angle1Radians;		/// angle between vertical plane and wedge contact face
  double height;
  double width;			/// width at the lowest horizontal plane
  double length;
  fmatvec::Vec wagonConnectionPointPosition;	/// three dimensional vector with the position of the connection point with wagon
  MBSim::FixedRelativeFrame *wagonConnectionFrame;	/// frame at the connection with the wagon
  /// Frame for geometry construction.
  /// Initialy it is placed at the center of the bounding box (mid-height, -width, and - length)
  MBSim::FixedRelativeFrame *geometryReferenceFrame;
  MBSim::Plate *leftWedgePlane;		/// contact plane for the left wedge (-X)
  MBSim::Plate *rightWedgePlane;	/// contact plane for the right wedge (+X)
  
  // private methods
  /// \brief Updates the position of the connection point
  /// \param none Uses the values already available from the variables
  void updateWagonConnectionPointPosition();
};

#endif // BOLSTER_H
