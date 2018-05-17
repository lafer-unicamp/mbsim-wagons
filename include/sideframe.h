/*
    Implements a sideframe
    Copyright (C) 2018 Leonardo Baruffaldi leonardo.baruffaldi@ifsp.edu.br

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


#ifndef SFRAME_H_
#define SFRAME_H_

#include "mbsim/objects/rigid_body.h"
#include "mbsim/frames/frame.h"
#include "mbsim/frames/fixed_relative_frame.h"
#include "mbsim/contours/contour.h"
#include "mbsim/contours/plane.h"
#include "mbsim/utils/rotarymatrices.h"
#include "openmbvcppinterface/ivbody.h"

class Sideframe : public MBSim::RigidBody
{
public:
	/// Default constructor
	/// \brief Default sideframe constructor
	/// \param Name of the sideframe
	Sideframe( const std::string &name);

	/// Sets bolster window width
	void setWindowWidth ( double _w ){ sideFrameWidth = _w;};
	void setWindowHeight ( double _h ){ sideFrameHeight = _h;};
	void setWindowAngle ( double _a ){ angleSideFrame = _a;};
	void setWheelBase (double _b){ wheelBase = _b;};

	/// Configures contact contours
	void buildContour (void);

	/// Returns selected contour
	///
	MBSim::Contour* getSideframeContour(int i);

private:
	double sideFrameWidth;
	double sideFrameHeight;
	double angleSideFrame;
	double wheelBase;

	MBSim::Plane *sideframeL;
	MBSim::Plane *sideframeR;
	MBSim::Plane *ground;
};

#endif
