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


struct railData{
	double x,y,z;
	/*
	 * \brief track cant angle
	 */
	double su;
	/*
	 * \brief lateral offset
	 */
	double deltay;
	/*
	 * \brief vertical offset
	 */
	double deltaz;
	/*
	 * \brief rail camber
	 */
	double camberLeft;
};

/**
 * \brief Class that contains information related to global rail topological data
 */
class RailTopology{
public:
	RailTopology() : nominalGauge(1.0) {};
private:
	double nominalGauge;
	std::string filePath;

};


/**
 * \brief A class that creates a two-dimensional profile for rail from a point map
 */
class RailProfile : public ProfileContour {
public:
	RailProfile(const std::string &name="", const std::string &file_="", Frame *R=nullptr) :
		ProfileContour(name, file_,true, true, R) { }

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
		std::shared_ptr<OpenMBV::Extrusion> ombv =
				OpenMBV::ObjectFactory::create<OpenMBV::Extrusion>();
		// TODO: acertar colorização do trilho. Por algum motivo o compilador
		// esté reclamando que o diffuseColor não pode ser usado como função
		//ombv->setDiffuseColor(diffuseColor(0),diffuseColor(1),diffuseColor(2));
		ombv->setHeight(1.0);
		ombv->addContour(configurePoints());
		openMBVRigidBody = ombv;
	}

private:
	/**
	 * \brief PolygonPoint map with the profile points
	 */
	std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>> configurePoints();

	/**
	 * \brief The rail topological data
	 */
	RailTopology* topology;
};




} /* namespace MBSim */

#endif /* INCLUDE_RAIL_PROFILE_H_ */
