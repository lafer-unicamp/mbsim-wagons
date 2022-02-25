/*
 * WheelProfile.cpp
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

#include "../include/wheel_profile.h"

namespace MBSim {

	MBSIM_OBJECTFACTORY_REGISTERCLASS(MBSIM, WheelProfile)

void WheelProfile::init(InitStage stage, const InitConfigSet &config) {
	if(stage==preInit) {
		//sign = solid?1:-1;
		readInputFile();
	}
	else if(stage==plotting) {
		if(plotFeature[openMBV] && openMBVRigidBody) {

		}
	}
	ProfileContour::init(stage, config);
}

void OpenMBVRotation::initializeUsingXML(xercesc::DOMElement *e) {
	OpenMBVColoredBody::initializeUsingXML(e);
	// TODO finish implementation
}

void OpenMBVRotation::initializeObject(
		const std::shared_ptr<OpenMBV::Rotation>& object) {
	OpenMBVColoredBody::initializeObject(object);
	std::cout << points2d << std::endl;
	configurePoints(points2d);
	object->setContour(points);
}

std::shared_ptr<OpenMBV::Rotation> OpenMBVRotation::createOpenMBV(
		xercesc::DOMElement* e) {
	std::shared_ptr<OpenMBV::Rotation> object = OpenMBV::ObjectFactory::create<OpenMBV::Rotation>();
	if(e) initializeUsingXML(e);
	initializeObject(object);
	return object;
}

void OpenMBVRotation::configurePoints(fmatvec::MatVx2 points_) {
	points = std::make_shared<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>>();

	int n = points_.rows();

	points->push_back(OpenMBV::PolygonPoint::create(0.0,points_(0,0),0.));
	for (int i = 0;i < n; i++)
	{
		points->push_back(OpenMBV::PolygonPoint::create(points_(i,1),points_(i,0),0.));
	}
	points->push_back(OpenMBV::PolygonPoint::create(0.0,points_(n-1,0),0.));

}


} /* namespace MBSim */
