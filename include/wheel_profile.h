/*
 * WheelProfile.h
 * Copyright (C) 2019 MBSim Development Team
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 *  Created on: 28 de mar de 2019
 *      Author: leonardo
 *	   Contact: <e-mail>
 */

#ifndef INCLUDE_WHEEL_PROFILE_H_
#define INCLUDE_WHEEL_PROFILE_H_

#include "profilecontour.h"
#include "openmbvcppinterface/rotation.h"

namespace MBSim {

/**
 * \brief Implements a class to create rotations from a profile given by a set of points
 * Class OpenMBVRotation
 * Implemented because the openmbv_util.h did not implement this
 */
class OpenMBVRotation : public OpenMBVColoredBody {
protected:
	/**
	 * \brief Points in format (x,y).
	 */
	fmatvec::MatVx2 points2d;
	std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>> points;
public:
	OpenMBVRotation(fmatvec::MatVx2 points_,const fmatvec::Vec3 &dc="[-1;1;1]", double tp=0) :
		OpenMBVColoredBody(dc,tp), points2d(points_) {}
	void initializeUsingXML(xercesc::DOMElement *element);
	void initializeObject(const std::shared_ptr<OpenMBV::Rotation> &object);
	/**
	 * \brief Configures the points, which are input as a two column matrix, to OpenMBVPolygonPoint
	 */
	void configurePoints(fmatvec::MatVx2 points_);
	std::shared_ptr<OpenMBV::Rotation> createOpenMBV(xercesc::DOMElement* e=0);
};



class WheelProfile: public ProfileContour {
public:
	WheelProfile(const std::string &name="", const std::string &file_="", Frame *R=nullptr) :
		ProfileContour(name,file_,false,true,R) {  }

	/* INHERITED INTERFACE OF ELEMENT */
	virtual void init(InitStage stage, const InitConfigSet &config);
	/***************************************************/

	/* INHERITED INTERFACE OF CONTOUR */
	     // virtual fmatvec::Vec3 evalKs(const fmatvec::Vec2 &zeta);
	     // virtual fmatvec::Vec3 evalKt(const fmatvec::Vec2 &zeta);
	     // virtual fmatvec::Vec3 evalParDer1Kn(const fmatvec::Vec2 &zeta);
	     // virtual fmatvec::Vec3 evalParDer1Ku(const fmatvec::Vec2 &zeta);
	     // virtual fmatvec::Vec2 evalZeta(const fmatvec::Vec3& WrPoint);
	/***************************************************/

	BOOST_PARAMETER_MEMBER_FUNCTION( (void), enableOpenMBV, tag, (optional (diffuseColor,(const fmatvec::Vec3&),"[-1;1;1]")(transparency,(double),0))) {
		OpenMBVRotation ombv(points,diffuseColor,transparency);
		openMBVRigidBody=ombv.createOpenMBV();
	}

protected:


};

} /* namespace MBSim */

#endif /* INCLUDE_WHEEL_PROFILE_H_ */
