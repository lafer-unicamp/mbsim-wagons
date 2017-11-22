/*
    Implements wagon freight box object
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


#ifndef WAGONBOX_H
#define WAGONBOX_H

#include <mbsim/objects/rigid_body.h>

#include "fmatvec/fmatvec.h"
#include "mbsim/frames/frame.h"
#include "mbsim/frames/fixed_relative_frame.h"

#include "openmbvcppinterface/cuboid.h"


class WagonBox : public MBSim::RigidBody
{
public:
  /// Default constructor
  WagonBox(const std::string& name);
  
  /// \brief sets geometry reference frame position at the center of gravity frame
  /// \param r_ 3x3 position vector
  void setGeometryReferenceFramePosition(fmatvec::Vec r_);
  /// \brief alternate implementation using double
  /// \param x_
  /// \param y_
  /// \param z_
  void setGeometryReferenceFramePosition(double x_, double y_, double z_);
  
  /// \brief sets front bolster connection frame position with respect to the center of gravity
  /// \param r_ 3x3 position vector
  void setFrontBolsterConnectionPosition(fmatvec::Vec r_);
  /// \brief alternate implementation using double
  /// \param x_
  /// \param y_
  /// \param z_
  void setFrontBolsterConnectionPosition(double x_, double y_, double z_);
  
  /// \brief sets rear bolster connection frame position with respect to the center of gravity
  /// \param r_ 3x3 position vector
  void setRearBolsterConnectionPosition(fmatvec::Vec r_);
  /// \brief alternate implementation using double
  /// \param x_
  /// \param y_
  /// \param z_
  void setRearBolsterConnectionPosition(double x_, double y_, double z_);
  
  /// \brief activates OpenMBV representation
  /// \param enable Boolean flag
  void enableOpenMBV(bool enable);
  
  void setLength(double l_){ length = l_; }
  void setWidth(double w_){ width = w_; }
  void setHeight(double h_){ height = h_; }
  
private:
  double length; /// \brief X dimension
  double height; /// \brief Y dimension
  double width; /// \brief Z dimension
  /// \brief frame to reference geometry
  MBSim::FixedRelativeFrame *geometryReferenceFrame;
  /// \brief connection with front bolster frame
  MBSim::FixedRelativeFrame* frontBolsterConnectionFrame;
  /// \brief connection with rear bolster frame
  MBSim::FixedRelativeFrame* rearBolsterConnectionFrame;
};

#endif // WAGONBOX_H
