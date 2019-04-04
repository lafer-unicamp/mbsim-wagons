/*
 * RailProfile.h
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

#ifndef INCLUDE_RAIL_PROFILE_H_
#define INCLUDE_RAIL_PROFILE_H_

#include "profilecontour.h"
#include "openmbvcppinterface/extrusion.h"

namespace MBSim {


/**
 * \brief A class that creates a two-dimensional profile for rail from a point map
 */
class RailProfile: public ProfileContour {
public:
	RailProfile(const std::string &name, const std::string &file_, Frame *R=0) :
		ProfileContour(name, file_,R) { readInputFile(); }

	/* INHERITED INTERFACE OF ELEMENT */
	std::string getType() const { return "Rail Profile Contour"; }
	virtual void init(InitStage stage);
	/***************************************************/

	BOOST_PARAMETER_MEMBER_FUNCTION( (void), enableOpenMBV, tag, (optional (diffuseColor,(const fmatvec::Vec3&),"[-1;1;1]")(transparency,(double),0))) {
		std::shared_ptr<OpenMBV::Extrusion> ombv =
				OpenMBV::ObjectFactory::create<OpenMBV::Extrusion>();
		ombv->setDiffuseColor(diffuseColor(0),diffuseColor(1),diffuseColor(2));
		ombv->setHeight(1.0);
		ombv->addContour(configurePoints());
		openMBVRigidBody = ombv;
	}

private:
	/**
	 * \brief PolygonPoint map with the profile points
	 */
	std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>> configurePoints();
};

} /* namespace MBSim */

#endif /* INCLUDE_RAIL_PROFILE_H_ */
