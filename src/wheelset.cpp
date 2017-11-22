/*
 * wheelset.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: leonardo
 */

#include "wheelset.h"

Wheelset::Wheelset(const std::string& name) : MBSim::RigidBody(name), track(1.0)
{
  // TODO Auto-generated constructor stub

}

void Wheelset::enableOpenMBV()
{
	std::shared_ptr<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>> points =
					std::make_shared<std::vector<std::shared_ptr<OpenMBV::PolygonPoint>>>();
//  std::vector<OpenMBV::PolygonPoint*>* points = new std::vector<OpenMBV::PolygonPoint*>(7);
  points->push_back(OpenMBV::PolygonPoint::create(.1,-0.067,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.4344,-0.067,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.4572,0.0,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.4593,0.043,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.4824,0.043,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.4824,0.077,0.));
  points->push_back(OpenMBV::PolygonPoint::create(.1,0.077,0.));

  std::shared_ptr<OpenMBV::Rotation> openMBVwheel1 =
  				OpenMBV::ObjectFactory::create<OpenMBV::Rotation>();
  openMBVwheel1->setName("Roda 1");;
  openMBVwheel1->setContour(points);
  openMBVwheel1->setInitialRotation(M_PI/2,0.,0.);
  openMBVwheel1->setInitialTranslation(0.,0.,track/2);
//
  std::shared_ptr<OpenMBV::Rotation> openMBVwheel2 =
		  OpenMBV::ObjectFactory::create<OpenMBV::Rotation>();
  openMBVwheel1->setName("Roda 2");
  openMBVwheel2->setContour(points);
  openMBVwheel2->setInitialRotation(3*M_PI/2,0.,0.);
  openMBVwheel2->setInitialTranslation(0.,0.,-track/2);

  std::shared_ptr<OpenMBV::Frustum> openMBVshaft =
		  OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
  openMBVshaft->setName("Eixo");
  openMBVshaft->setBaseRadius(0.125);
  openMBVshaft->setTopRadius(0.125);
  openMBVshaft->setHeight(track);
  openMBVshaft->setInitialTranslation(0,0,track/2.0);
  openMBVshaft->setDiffuseColor(0.25,0.25,0.25);
//
  std::shared_ptr<OpenMBV::CompoundRigidBody> openMBVwheelset =
		  OpenMBV::ObjectFactory::create<OpenMBV::CompoundRigidBody>();
  openMBVwheelset->addRigidBody(openMBVwheel1);
  openMBVwheelset->addRigidBody(openMBVwheel2);
  openMBVwheelset->addRigidBody(openMBVshaft);
//
  this->setOpenMBVRigidBody(openMBVwheelset);
}
