/*
    Creates a wedged rigid body
    Copyright (C) 2013  Leonardo Baruffaldi leobart@fem.unicamp.br

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#include "wedge.h"
#include <fmatvec/fmatvec.h>

using namespace std;

Wedge::Wedge(const std::string &name) : RigidBody(name), height(1.0), depth(1.0), 
angle1Radians(0.5), angle2Radians(0.0) { 
  pointTable = new fmatvec::Mat();
  gRefFramePosition = new fmatvec::Vec(3,fmatvec::INIT,0.0);
  geometryReferenceFrame = new MBSim::FixedRelativeFrame("RG",*gRefFramePosition,
							 fmatvec::SqrMat(3,fmatvec::EYE),
							 this->getFrameOfReference());
  this->addFrame(geometryReferenceFrame); /// TODO erro na modificação da posição
}

void Wedge::setHeight(double h_)
{
  height = h_;
  updateGrefFramePosition();
}

void Wedge::setAngles(fmatvec::Vec angs_)
{
  angle1Radians = angs_(0);
  angle2Radians = angs_(1);
  updateGrefFramePosition();
}

void Wedge::setDepth(double d_)
{
  depth = d_;
  updateGrefFramePosition();
}

void Wedge::buildContour(void)
{  
  int numberOfPoints = 8;
  
  pointTable->resize(numberOfPoints,3);
  
  double tip = 0.05;
  
  /// First point (upper left)
  (*pointTable)(0,0) = - tip * height * tan(angle2Radians);
  (*pointTable)(0,1) = - tip * height;
  (*pointTable)(0,2) = 0.0;
  
  /// Second point (lower left)
  (*pointTable)(1,0) = - height * tan(angle2Radians);
  (*pointTable)(1,1) = -height;
  (*pointTable)(1,2) = 0.0;
  
  /// Third point (lower right)
  (*pointTable)(2,0) = height * tan(angle1Radians);
  (*pointTable)(2,1) = -height;
  (*pointTable)(2,2) = 0.0;
  
  /// Fourth point (upper right)
  (*pointTable)(3,0) = tip * height * tan(angle1Radians);
  (*pointTable)(3,1) = - tip * height;
  (*pointTable)(3,2) = 0.0;
  

  /// Fifth, sixth, seventh, and eigth points (depth)
  for (int i = 0; i < numberOfPoints/2; i++){
	  (*pointTable)(i+numberOfPoints/2,0) = (*pointTable)(i,0);
	  (*pointTable)(i+numberOfPoints/2,1) = (*pointTable)(i,1);
	  (*pointTable)(i+numberOfPoints/2,2) = depth;
  }


  for(int i=0; i<numberOfPoints; i++) {
	  std::stringstream s;
	  std::stringstream frame_name;
	  s << i+1;
	  frame_name << this->name << "-frame_" << i+1;
	  fmatvec::Vec *position = new fmatvec::Vec(3,fmatvec::INIT,0.0);
	  position->operator()(0) = (*pointTable)(i,0);
	  position->operator()(1) = (*pointTable)(i,1);
	  position->operator()(2) = (*pointTable)(i,2);
	  addFrame(new MBSim::FixedRelativeFrame(frame_name.str(),
					  *position,
					  fmatvec::SqrMat(3, fmatvec::EYE),
					  geometryReferenceFrame));
	  MBSim::Point *point = new MBSim::Point(s.str(),
			  getFrame(frame_name.str()));
	  point->enableOpenMBV();
	  addContour(point);
  }

}

void Wedge::enableOpenMBV(bool enable)
{
  if(enable){
    shared_ptr<vector<shared_ptr<OpenMBV::PolygonPoint>>> vecPoint =
    		make_shared<vector<shared_ptr<OpenMBV::PolygonPoint>>>();
    double x,y,z;
    for (unsigned i=0; i < 4; i++){
	x = this->getPointTable()->operator()(i,0);
	y = this->getPointTable()->operator()(i,1);
	z = this->getPointTable()->operator()(i,2);
	vecPoint->push_back(OpenMBV::PolygonPoint::create(x,y,z));
    }

    std::shared_ptr<OpenMBV::Extrusion> openMBVWedge =
    		OpenMBV::ObjectFactory::create<OpenMBV::Extrusion>();
    openMBVWedge->setHeight(depth);
    openMBVWedge->setInitialTranslation(gRefFramePosition->operator()(0),
					gRefFramePosition->operator()(1),
					gRefFramePosition->operator()(2)
    );
    openMBVWedge->setDiffuseColor(0.25,0.25,0.35);
    openMBVWedge->addContour(vecPoint);
    this->setOpenMBVRigidBody(openMBVWedge);
  }
}

void Wedge::updateGrefFramePosition()
{
      double cSecArea;
      cSecArea = pow(height,2)*(tan(angle1Radians)+tan(angle2Radians))/2;
      gRefFramePosition->operator()(0) = (-pow(height,3)*pow(tan(angle1Radians),2)/6 +
					  pow(height,3)*pow(tan(angle2Radians),2)/6)/cSecArea;
      gRefFramePosition->operator()(1) = 2./3 * height;
      gRefFramePosition->operator()(2) = - depth / 2;
      geometryReferenceFrame->setRelativePosition(*gRefFramePosition);
}
