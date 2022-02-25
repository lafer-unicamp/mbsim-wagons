/*
 * profilecontour.cpp
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

#include "profilecontour.h"

namespace MBSim {

void ProfileContour::readInputFile() {

	std::ifstream dataPointsFile(this->file);

	// counts the number of lines on the file
	numberOfPoints = 0;
	std::string unused;
	while ( std::getline(dataPointsFile, unused) )
		++numberOfPoints;

	this->points.resize(numberOfPoints,2);

	std::string line;
	// reads line by line
	dataPointsFile.clear();  // sets the eof flag to good
	dataPointsFile.seekg(0, std::ios::beg);

	for(int row = 0; row < numberOfPoints; ++row)
	{
		line.clear();
		std::getline(dataPointsFile, line);
		if ( !dataPointsFile.good() )
			break;

		std::stringstream iss(line);

		for (int col = 0; col < 2; ++col)
		{
			std::string val;
			std::getline(iss, val, ';');

			if ( !iss.good() ){
				//				break;
			}

			std::stringstream convertor(val);
			convertor >> points(row,col);
		}
	}

	dataPointsFile.close();

	//if ( numberOfPoints > 0 ) this->setConvexIndexes();
}

void ProfileContour::setFile(std::string& filePath_) {
	file = filePath_;
	readInputFile();
}

void ProfileContour::init(InitStage stage, const InitConfigSet &config) {
	if(stage==preInit) {
		//sign = solid?1:-1;
		readInputFile();
	}
	else if(stage==plotting) {
		if(plotFeature[openMBV] && openMBVRigidBody) {

		}
	}
	RigidContour::init(stage, config);
}


/**
 * \brief separates the contours into convex sets by listing the first and last indexes
 * of the convex set points. This routine is called once when the profile is set, for it
 * is supposed constant. TODO: continuously deformable wheels would require this to be updated
 * on the fly
 */
void ProfileContour::setConvexIndexes() {

	unsigned int k = 0;
	unsigned int nw = getNumberOfPoints();

	this->convexSetIndexes.resize(nw,2);

	this->convexSetIndexes(0,0) = 0;
	this->convexSetIndexes(0,1) = 0;

	fmatvec::RowVec2 v1, v2;

	for (unsigned int i = 0;i < (nw-2); i++){
		v1 = points.row(i+1) - points.row(i);
		v2 = points.row(i+2) - points.row(i+1);
		v2 = (cw) ? -v2 : v2;
		if (v1(0)*v2(1) - v1(1)*v2(0) > -1e-8)
		{
			this->convexSetIndexes(k,1) = i+2;
		}
		else
		{
			k = k+1;
			this->convexSetIndexes(k,0) = this->convexSetIndexes(k-1,1);
			this->convexSetIndexes(k,1) = this->convexSetIndexes(k-1,1) + 1;
		}
	}

	// TODO redimensionar matriz sem apagar os elementos
	fmatvec::MatVx2 temp;
	temp.resize(k+1,2);
	for (unsigned int i = 0;i<=k;i++){
		temp(i,0) = convexSetIndexes(i,0);
		temp(i,1) = convexSetIndexes(i,1);
	}


	this->convexSetIndexes.resize(k,2);
	convexSetIndexes <<= temp;

}

} /* namespace MBSim */
