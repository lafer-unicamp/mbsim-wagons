/*
    Creates a wedged rigid body
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

///TODO include center of mass displacement
///TODO include reference frame for the geometry

#ifndef WEDGE_H
#define WEDGE_H

#include "mbsim/objects/rigid_body.h"
#include "mbsim/frames/frame.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "mbsim/contours/contour.h"
#include "mbsim/contours/point.h"
#include "mbsim/contours/compound_contour.h"
#include "fmatvec/vector.h"
#include "fmatvec/matrix.h"
#include <math.h>


#include "openmbvcppinterface/extrusion.h"

class Wedge : public MBSim::RigidBody
{
public:
  /// Default constructor
  /// \brief default class contructor
  /// \param name Name of the wedge object
  Wedge(const std::string& name);
  
  /// Setters | Getters
  
  void setHeight(double h_);
  void setDepth(double d_);
  void setAngles(fmatvec::Vec angs_); /// \param angs_ vector with two elements for angles 1 and 2
  
  double getHeight(void){return height;};
  
  fmatvec::Mat* getPointTable(void) { return pointTable; }
  fmatvec::Vec* getGrefFramePosition(void) { return gRefFramePosition; };
  
  /// Uses geometric information to form the compound contour
  /// Fills point table matrix and the compound contour
  void buildContour(void);
  
  /// \param enable Boolean to activate OpenMBV
  void enableOpenMBV(bool enable);

private:
  double height; /// vertical height of the wedge
  double depth; /// depth for extrusion
  double angle1Radians; /// angle between vertical plane and rigth face
  double angle2Radians; /// angle between vertical plane and left face
  MBSim::FixedRelativeFrame *geometryReferenceFrame;
  fmatvec::Vec *gRefFramePosition; /// position vector of the reference frame with respect with the center of mass
  fmatvec::Mat *pointTable; /// matrix Nx3 with points coordinates
  void updateGrefFramePosition(); /// update geometry reference frame position
};
#endif // WEDGE_H
