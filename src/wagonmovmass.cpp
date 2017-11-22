/*
    Wagon with moving mass
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


#include "wagonmovmass.h"

WagonMovMass::WagonMovMass(const std::string& name, double sFactor_): WagonSimple(name),
splitFactor(sFactor_)
{  
  movingMass = new RigidBody("Moving mass");
  movingMass->setFrameOfReference(wagonbox->getFrame("C"));
  movingMass->setFrameForKinematics( movingMass->getFrame("C") );
  movingMass->setTranslation( new LinearTranslation<fmatvec::VecV> (  ( "[1,0;0,1;0,0]" ) ) );
  fmatvec::Vec initialPosition(3,fmatvec::INIT,0.0);
  initialPosition(0) = 1e-4;
  movingMass->setGeneralizedInitialPosition(initialPosition);
  addObject(movingMass);
#ifdef HAVE_OPENMBVCPPINTERFACE
  movingMass->getFrame("C")->enableOpenMBV();
#endif
  
  /**
   * Definition of the spring connecting the moving mass to the wagon's body
   */
  spring = new MBSim::SpringDamper("Moving mass spring");
  // added small initial deformation to correct initialization error:
  // the initial direction of the force could not be calculated
  spring->connect(wagonbox->getFrame("C"),movingMass->getFrame("C"));
  spring->setForceFunction(new MBSim::LinearSpringDamperForce(10e3,0.01));
  spring->setUnloadedLength(0.00001);
  addLink(spring);
//#ifdef HAVE_OPENMBVCPPINTERFACE
//  OpenMBV::CoilSpring *openMBVSpringM = new OpenMBV::CoilSpring();
//  openMBVSpringM->setCrossSectionRadius ( .005 );
//  openMBVSpringM->setNumberOfCoils ( 5 );
//  openMBVSpringM->setSpringRadius ( 0.01 );
//  spring->setOpenMBVSpring ( openMBVSpringM );
//#endif

  /**
   * Definition of the translational joint with wagon's body
   */
  translational = new MBSim::Joint("Joint: Wagon to moving mass");
  translational->setForceDirection( fmatvec::Mat("[0;1;0]") );
  translational->setForceLaw( new MBSim::BilateralConstraint() );
  translational->connect(wagonbox->getFrame("C"),movingMass->getFrame("C"));
  addLink(translational);
}

void WagonMovMass::setTotalMass(double m_)
{
  this->WagonGroup::setTotalMass(m_);
  wagonbox->setMass(splitFactor * m_);
  movingMass->setMass((1.0-splitFactor)*m_);
}
