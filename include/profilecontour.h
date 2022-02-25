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

		ProfileContour(const std::string& name="", const std::string& file_="",
		bool cw_ = false, bool solid_=true, Frame *R=nullptr) :
		RigidContour(name,R), file(file_), cw(cw_) {
			readInputFile();
			setConvexIndexes();
		}

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

		/* GETTER / SETTER */
		void setFile(std::string& filePath_);
		double getNumberOfPoints(){return numberOfPoints;}
		fmatvec::MatVx2 getPoints() {return points;}
		fmatvec::MatVx2 getConvexSetIndexes() {return convexSetIndexes; }
		/***************************************************/


		// BOOST_PARAMETER_MEMBER_FUNCTION( (void), enableOpenMBV, tag, (optional (diffuseColor,(const fmatvec::Vec3&),fmatvec::Vec3(std::vector<double>{-1,1,1}))(transparency,(double),0)(pointSize,(double),0)(lineWidth,(double),0))) {
		// 	OpenMBVLine ombv(1,diffuseColor,transparency,pointSize,lineWidth);
		// 	openMBVRigidBody=ombv.createOpenMBV();
		// }

		// void initializeUsingXML(xercesc::DOMElement *element) override;


		void readInputFile();



	protected:

		/**
		* \brief a matrix that contains the first and last indexes
		* of the convex subsets of the profile.
		**/
		fmatvec::MatVx2 convexSetIndexes;


		/**
		* \brief number of points read from data file
		*/
		int numberOfPoints;

		/**
		* \brief coordinates of points
		*/
		fmatvec::MatVx2 points;

		/**
		* \brief name of the input file with the data points
		*/
		std::string file;

		/**
		* \brief setConvexIndexes sets the variable convexSetIndexes by
		* checking the convexity of consecutive lines inside the profile
		* \param clockwise This flag must be true if the profile inner
		* orientation is supposed to be clockwise. By default the algorithm searches
		* for convex sets counter-clockwisely.
		*/
		public :
		void setConvexIndexes();


	private:

		bool solid;
		/**
		* true if the profile is supposed to be closed clockwise
		*/
		bool cw;

	};


} /* namespace MBSim */

#endif /* INCLUDE_PROFILECONTOUR_H_ */
