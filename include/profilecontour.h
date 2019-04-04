/*
 * profilecontour.h
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
 *  Created on: 26 de mar de 2019
 *      Author: leonardo
 *	   Contact: <e-mail>
 */

#ifndef INCLUDE_PROFILECONTOUR_H_
#define INCLUDE_PROFILECONTOUR_H_

#include <mbsim/contours/rigid_contour.h>
#include "mbsim/utils/boost_parameters.h"
#include <mbsim/utils/openmbv_utils.h>
#include <fstream>
#include <sstream>


namespace MBSim {


/**
 * \brief general class for contours defined by 2d profile points
 * \author Leonardo Baruffaldi
 * \date 2019-03-22 initial commit
 */
class ProfileContour : public MBSim::RigidContour {
public:

	ProfileContour(const std::string& name="", const std::string& file_="", bool solid_=true, Frame* R=0) :
		RigidContour(name,R), file(file_) { readInputFile(); }

	/* INHERITED INTERFACE OF ELEMENT */
	std::string getType() const { return "Profile Contour"; }
	virtual void init(InitStage stage);
	/***************************************************/

	/* INHERITED INTERFACE OF CONTOUR */
	//      virtual fmatvec::Vec3 evalKs(const fmatvec::Vec2 &zeta);
	//      virtual fmatvec::Vec3 evalKt(const fmatvec::Vec2 &zeta) { return Kt; }
	//      virtual fmatvec::Vec3 evalParDer1Kn(const fmatvec::Vec2 &zeta);
	//      virtual fmatvec::Vec3 evalParDer1Ku(const fmatvec::Vec2 &zeta);
	//      virtual fmatvec::Vec2 evalZeta(const fmatvec::Vec3& WrPoint);
	/***************************************************/

	/* GETTER / SETTER */
	void setFile(std::string& filePath_);
	double getNumberOfPoints(){return numberOfPoints;}
	fmatvec::MatVx2 getPoints() {return points;}

	// void setSolid(bool solid_=true) { solid = solid_; }
	// bool getSolid() const { return solid; }
	/***************************************************/

//	BOOST_PARAMETER_MEMBER_FUNCTION( (void), enableOpenMBV, tag, (optional (diffuseColor,(const fmatvec::Vec3&),"[-1;1;1]")(transparency,(double),0))) {
//		OpenMBVRotation ombv(points,diffuseColor,transparency);
//		openMBVRigidBody=ombv.createOpenMBV();
//	}

	// virtual void initializeUsingXML(xercesc::DOMElement *element);


	void readInputFile();

protected:


	/**
	 * \brief number of points read from data file
	 */
	int numberOfPoints;

	/**
	 * \brief coordinates of points
	 */
	fmatvec::MatVx2 points;


private:
/**
	 * \brief name of the input file with the data points
	 */
	std::string file;
	bool solid;




};


} /* namespace MBSim */

#endif /* INCLUDE_PROFILECONTOUR_H_ */
