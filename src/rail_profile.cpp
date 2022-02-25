/*
 * RailProfile.cpp
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
 *  Created on: 29 de mar de 2019
 *      Author: leonardo
 *	   Contact: <e-mail>
 */

#include "../include/rail_profile.h"

namespace MBSim {

	MBSIM_OBJECTFACTORY_REGISTERCLASS(MBSIM, RailProfile)

void RailProfile::init(InitStage stage, const InitConfigSet &config) {
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

std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint> > > RailProfile::configurePoints() {


	std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint> > > polygonPoints =
			std::make_shared<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>>();
	int n = points.rows();


	for(int i = 0;i < n;i++){
		polygonPoints->push_back(OpenMBV::PolygonPoint::create(points(i,0),points(i,1),0.0));
	}

	return polygonPoints;
}

} /* namespace MBSim */
