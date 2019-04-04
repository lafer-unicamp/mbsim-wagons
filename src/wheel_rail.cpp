/*
 * wheel_rail.cpp
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
 *      Author: Leonardo Baruffaldi
 *	   Contact: leonardo.baruffaldi@ifsp.edu.br
 */

#include  "wheel_rail.h"

namespace MBSim{

ContactKinematicsWheelRail::~ContactKinematicsWheelRail() {



}

void ContactKinematicsWheelRail::assignContours(
		const std::vector<Contour*>& contour) {
	if(dynamic_cast<RailProfile*>(contour[0])){
		irail = 0;
		iwheel = 1;
		rail = static_cast<RailProfile*>(contour[0]);
		wheel = static_cast<WheelProfile*>(contour[1]);
	} else {
		irail = 1;
		iwheel = 0;
		rail = static_cast<RailProfile*>(contour[1]);
		wheel = static_cast<WheelProfile*>(contour[0]);
	}
}

void ContactKinematicsWheelRail::updateg(double& g,
		std::vector<ContourFrame*>& cFrame, int index) {
	// TODO implement GJK distance algorithm
}

void ContactKinematicsWheelRail::updatewb(fmatvec::Vec& wb, double g,
		std::vector<ContourFrame*>& cFrame) {
}

} /* END namespace MBSim */


