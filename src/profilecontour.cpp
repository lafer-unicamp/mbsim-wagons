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

	std::ifstream dataPointsFile(file);

	// counts the number of lines on the file
	this->numberOfPoints = 0;
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
}

void ProfileContour::setFile(std::string& filePath_) {
	file = filePath_;
	readInputFile();
}

void ProfileContour::init(InitStage stage) {
	if(stage==preInit) {
		//sign = solid?1:-1;
		readInputFile();
	}
	else if(stage==plotting) {
		if(plotFeature[openMBV]==enabled && openMBVRigidBody) {

		}
	}
	RigidContour::init(stage);
}

} /* namespace MBSim */


