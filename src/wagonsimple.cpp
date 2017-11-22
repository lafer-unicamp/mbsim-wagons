/*
 * wagonsimple.cpp
 *
 *  Created on: Nov 7, 2013
 *      Author: leonardo
 */

#include "wagonsimple.h"

WagonSimple::WagonSimple(const std::string& name) :
    WagonGroup(name)
{
  wagonbox = new MBSim::RigidBody("Wagon box");
  this->addObject(wagonbox);

  geometryReferenceFrame = new MBSim::FixedRelativeFrame("GF",
      fmatvec::Vec3(fmatvec::INIT, 0.0), fmatvec::SqrMat(3, fmatvec::EYE),
      wagonbox->getFrameC());
  wagonbox->addFrame(geometryReferenceFrame);

  frontBolsterConnection = new MBSim::FixedRelativeFrame("FBC",
      fmatvec::Vec3(fmatvec::INIT, 0.0), fmatvec::SqrMat(3, fmatvec::EYE),
      wagonbox->getFrameC());
  wagonbox->addFrame(frontBolsterConnection);

  rearBolsterConnection = new MBSim::FixedRelativeFrame("RBC",
      fmatvec::Vec3(fmatvec::INIT, 0.0), fmatvec::SqrMat(3, fmatvec::EYE),
      wagonbox->getFrameC());
  wagonbox->addFrame(rearBolsterConnection);

}

void
WagonSimple::setTotalMass(double totalMass)
{
  WagonGroup::setTotalMass(totalMass);
  wagonbox->setMass(totalMass);
}

void
WagonSimple::setGeometryReferenceFramePosition(double x_, double y_, double z_)
{
  fmatvec::Vec3 pos_(fmatvec::INIT, 0.0);
  pos_(0) = x_;
  pos_(1) = y_;
  pos_(2) = z_;
  this->setGeometryReferenceFramePosition(pos_);
}

void
WagonSimple::setGeometryReferenceFramePosition(fmatvec::Vec3 pos_)
{
  geometryReferenceFrame->setRelativePosition(pos_);
}

void
WagonSimple::setFrontBolsterConnectionPosition(double x_, double y_, double z_)
{
  fmatvec::Vec3 pos_(fmatvec::INIT, 0.0);
  pos_(0) = x_;
  pos_(1) = y_;
  pos_(2) = z_;
  this->setFrontBolsterConnectionPosition(pos_);
}

void
WagonSimple::setFrontBolsterConnectionPosition(fmatvec::Vec3 pos_)
{
  frontBolsterConnection->setRelativePosition(pos_);
}

void
WagonSimple::setRearBolsterConnectionPosition(double x_, double y_, double z_)
{
  fmatvec::Vec3 pos_(fmatvec::INIT, 0.0);
  pos_(0) = x_;
  pos_(1) = y_;
  pos_(2) = z_;
  this->setRearBolsterConnectionPosition(pos_);
}

void
WagonSimple::setRearBolsterConnectionPosition(fmatvec::Vec3 pos_)
{
	rearBolsterConnection->setRelativePosition(pos_);
}

void WagonSimple::enableOpenMBV(bool enable)
{
	if(enable){
		std::shared_ptr<OpenMBV::Cuboid> wagonboxMBV=
				OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
		wagonboxMBV->setDynamicColor(4);
		wagonboxMBV->setLength(length,height,width);
		wagonbox->setOpenMBVRigidBody(wagonboxMBV);
	};
}
