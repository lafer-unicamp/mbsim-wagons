/*
    Wagon with embedded sloshing model acc. to Graham-Rodriguez
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

    References:
    Graham, Rodriguez, The Characteristics of Fuel Motion Which Affect Airplane
    Dynamics,J. App. Mechanics, 19(3), 1952.
 */


#include "wagonsloshing.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

WagonSloshing::WagonSloshing(const std::string& name): WagonSimple(name),
fillRatio(1.0), n(2){


}

void WagonSloshing::setTotalMass(double m_)
{
	this->WagonGroup::setTotalMass(m_);
	wagonbox->setMass(m_);
}

void WagonSloshing::Initialize() {
	/*
	 * TODO correct wagon total mass = wagon + solid mass + sloshing
	 * TODO implement solid mass moment of inertia
	 *
	 */

	// output to file
	std::ofstream fileo;
	fileo.open("wagon_properties.txt");


	// Initialize variables
	fluidRho = 780.; // kg/m^3
	fileo << "Fluid density: " << fluidRho << std::endl;
	totalFluidMass = this->getHeight() * this->getWidth() * this->getLength() *
			fillRatio * fluidRho;
	fileo << "Tank length: " << this->getLength() << std::endl;
	fileo << "Tank height: " << this->getHeight() << std::endl;
	fileo << "Tank width: " << this->getWidth() << std::endl;
	fileo << "Filled volume: " << this->getFillRatio() * 100 << "%" << std::endl;
	fileo << "Total mass of fluid: " << totalFluidMass << std::endl;
	aspectRatio = this->getHeight() / this->getLength() * fillRatio;
	fileo << "Aspect ratio: " << aspectRatio << std::endl;

	massContainer.resize(n);
	springs.resize(n);
	joints.resize(n);
	frames.resize(n);

	/*
	 * osstringstream to form element names
	 */
	std::ostringstream convert;


	fmatvec::Vec initialPosition(3,fmatvec::INIT,0.0);

	// SLOSHING MASSES
	for (unsigned i=0;i<n;i++){
		convert.str(std::string());
		convert << "Sloshing mass " << i;
		// new sloshing mass
		massContainer[i] = new MBSim::RigidBody(convert.str());
		massContainer[i]->setMass(8*totalFluidMass*tanh((2*i+1)*M_PI*aspectRatio) /
				(pow(2*i+1,3)*pow(M_PI,3)*aspectRatio));
		fileo << "m" << i << ": " << massContainer[i]->getMass() << std::endl;
		massContainer[i]->setFrameOfReference(wagonbox->getFrameC());
		massContainer[i]->setFrameForKinematics(massContainer[i]->getFrameC());
		massContainer[i]->setTranslation( new MBSim::LinearTranslation<fmatvec::VecV>("[1,0,0;0,1,0;0,0,1]"));

		initialPosition(1) = this->getHeight()*(0.5-(tanh((2*i+1)*M_PI*aspectRatio/2))/
				((2*i+1)*M_PI*aspectRatio/2));
		fileo << "z" << i << ": " << initialPosition(1) << std::endl;
		massContainer[i]->setGeneralizedInitialPosition(initialPosition);

		convert.str(std::string());
		convert << "Frame: Sloshing " << i;
		frames[i] = new MBSim::FixedRelativeFrame(convert.str(),initialPosition,
				fmatvec::SqrMat(3,fmatvec::EYE));
		// adds mass to the wagon
		wagonbox->addFrame(frames[i]);
		addObject(massContainer[i]);
		massContainer[i]->getFrame("C")->enableOpenMBV();

		// SPRINGS
		convert.str(std::string());
		convert << "Spring: sloshing " << i;
		springs[i] = new MBSim::SpringDamper(convert.str());
		// added small initial deformation to correct initialization error:
		// the initial direction of the force could not be calculated
		springs[i]->setForceFunction(new MBSim::LinearSpringDamperForce(
				8*9.81*totalFluidMass*pow(tanh((2*i+1)*M_PI*aspectRatio),2) /
				(pow(2*i+1,2)*pow(M_PI,2)*this->getHeight()*fillRatio)
				,0.01));
		springs[i]->setUnloadedLength(0.001);
		fileo << "k" << i << ": " << 8*9.81*totalFluidMass*pow(tanh((2*i+1)*M_PI*aspectRatio),2) /
				(pow(2*i+1,2)*pow(M_PI,2)*this->getHeight()*fillRatio) << std::endl;

		springs[i]->connect(frames[i],massContainer[i]->getFrameC());

		addLink(springs[i]);


		/**
		 * Definition of the translational joints with wagon's body
		 */
		convert.str(std::string());
		convert << "Joint: translation: " << i;
		joints[i] = new MBSim::Joint(convert.str());
		joints[i]->setForceDirection( fmatvec::Mat("[0,0;1,0;0,1]") );
		joints[i]->setForceLaw( new MBSim::BilateralConstraint() );
		//joints[i]->setImpactForceLaw( new MBSim::BilateralImpact() );
		joints[i]->connect(frames[i],massContainer[i]->getFrameC());
		addLink(joints[i]);

	}// for


	// SOLID MASS

	solidMass = new MBSim::RigidBody("Solidified fluid");
	addObject(solidMass);


	// degrees of freedom
	solidMass->setFrameOfReference(wagonbox->getFrame("C"));
	solidMass->setFrameForKinematics(solidMass->getFrame( "C" ));
	solidMass->setTranslation( new MBSim::LinearTranslation<fmatvec::VecV>("[1,0;0,1;0,0]"));
	solidMass->setRotation( new MBSim::RotationAboutZAxis<fmatvec::VecV>() );

	// total mass
	temp = 0;
	for (unsigned i = 0;i<n;i++){
		temp = temp + massContainer[i]->getMass();
	}
	solidMass->setMass(totalFluidMass - temp);
	fileo << "M: " << solidMass->getMass() << std::endl;

	// loop to calculate the sloshing masses*position products
	temp = 0;
	for (unsigned i = 0;i<n;i++){
		temp = temp + massContainer[i]->getMass()*massContainer[i]->getq0()(1);
	}/* for */

	// initial position
	initialPosition(1) = -(1 / solidMass->getMass()) * temp;
	fileo << "Z: " << initialPosition(1) << std::endl;
	solidMass->setGeneralizedInitialPosition(initialPosition);
	solidMassConnectionFrame = new MBSim::FixedRelativeFrame("Frame: solid fluid",initialPosition,
			fmatvec::SqrMat(3,fmatvec::EYE));
	wagonbox->addFrame(solidMassConnectionFrame);

	// Inertia tensor
	solidMassInertia.resize(3);
	solidMassInertia.init(fmatvec::EYE);
	solidMassInertia(2,2) = totalFluidMass/12 * (pow(this->getLength(),2) +
			pow(this->getHeight()*fillRatio,2)) -
			solidMass->getMass()*pow(solidMass->getq0()(1),2);

	for (unsigned i=0; i<n; i++){
		solidMassInertia(2,2) = solidMassInertia(2,2) -
				massContainer[i]->getMass()*pow(massContainer[i]->getq0()(1),2);
	}/* for */

	solidMass->setInertiaTensor(solidMassInertia);

	temp = solidMass->getInertiaTensor()(2,2);

	fileo << "Solid mass inertia: " << solidMassInertia(2,2) << std::endl;

	/*
	 * Graphics
	 */
	solidMass->getFrame("C")->enableOpenMBV();

	/*
	 * rigid constraint to wagon
	 */
	rigidJoint = new MBSim::Joint("Joint: rigid");

	rigidJoint->setForceDirection( fmatvec::Mat("[1,0;0,1;0,0]") );
	rigidJoint->setForceLaw( new MBSim::BilateralConstraint() );
//	rigidJoint->setImpactForceLaw( new MBSim::BilateralImpact() );

	rigidJoint->setMomentDirection( fmatvec::Mat("[0;0;1]") );
	rigidJoint->setMomentLaw( new MBSim::BilateralConstraint() );
//	rigidJoint->setImpactMomentLaw( new MBSim::BilateralImpact() );

	rigidJoint->connect(solidMassConnectionFrame,solidMass->getFrameC());
	addLink(rigidJoint);

	fileo.close();
}

