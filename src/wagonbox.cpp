/*
    Implements wagon freight box object
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


#include "wagonbox.h"

WagonBox::WagonBox(const std::string& name): RigidBody(name), length(1.0), height(1.0), width(1.0)
{
  geometryReferenceFrame = new MBSim::FixedRelativeFrame("GR",fmatvec::Vec(3,fmatvec::INIT,0.0),fmatvec::SqrMat(3,fmatvec::EYE),this->getFrameC());
  this->addFrame(geometryReferenceFrame);
  
  frontBolsterConnectionFrame = new MBSim::FixedRelativeFrame("FBC",fmatvec::Vec(3,fmatvec::INIT,0.0),fmatvec::SqrMat(3,fmatvec::EYE),this->getFrameC());
  this->addFrame(frontBolsterConnectionFrame);
  
  rearBolsterConnectionFrame = new MBSim::FixedRelativeFrame("RBC", fmatvec::Vec(3,fmatvec::INIT,0.0),fmatvec::SqrMat(3,fmatvec::EYE),this->getFrameC());
  this->addFrame(rearBolsterConnectionFrame);
}

//setGeometryReferenceFramePosition
void WagonBox::setGeometryReferenceFramePosition(fmatvec::Vec r_)
{
  geometryReferenceFrame->setRelativePosition(r_);
}

void WagonBox::setGeometryReferenceFramePosition(double x_, double y_, double z_)
{
  fmatvec::Vec r_(3,fmatvec::INIT,0.0);
  r_(0) = x_;
  r_(1) = y_;
  r_(2) = z_;
  
  setGeometryReferenceFramePosition(r_);
}

//setFrontBolsterConnectionPosition
void WagonBox::setFrontBolsterConnectionPosition(fmatvec::Vec r_)
{
  frontBolsterConnectionFrame->setRelativePosition(r_);
}

void WagonBox::setFrontBolsterConnectionPosition(double x_, double y_, double z_)
{
  fmatvec::Vec r_(3,fmatvec::INIT,0.0);
  r_(0) = x_;
  r_(1) = y_;
  r_(2) = z_;
  
  setFrontBolsterConnectionPosition(r_);
}

//setRearBolsterConnectionPosition
void WagonBox::setRearBolsterConnectionPosition(fmatvec::Vec r_)
{
  rearBolsterConnectionFrame->setRelativePosition(r_);
}

void WagonBox::setRearBolsterConnectionPosition(double x_, double y_, double z_)
{
  fmatvec::Vec r_(3,fmatvec::INIT,0.0);
  r_(0) = x_;
  r_(1) = y_;
  r_(2) = z_;
  
  setRearBolsterConnectionPosition(r_);
}

void WagonBox::enableOpenMBV(bool enable)
{
  if(enable){
	std::shared_ptr<OpenMBV::Cuboid> wagonbox =
			  OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
    wagonbox->setDynamicColor(4);
    wagonbox->setLength(length,height,width);
    this->setOpenMBVRigidBody(wagonbox);
  };
}
