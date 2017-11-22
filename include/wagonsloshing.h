/*
 	Wagon with embedded sloshing model acc. to Graham-Rodriguez
    Copyright (C) 2014  Leonardo Baruffaldi leobart@fem.unicamp.br

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

#ifndef WAGONSLOSHING_H
#define WAGONSLOSHING_H

#include "wagonbox.h"
#include "wagonsimple.h"
#include "mbsim/links/spring_damper.h"
#include "mbsim/links/joint.h"
#include "mbsim/constitutive_laws/constitutive_laws.h"
#include "mbsim/frames/frame.h"
#include "mbsim/functions/kinetics/linear_spring_damper_force.h"
#include <math.h> // to use the tanh function
#include <vector>
#include <iostream>

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/coilspring.h>
#endif

class WagonSloshing : public WagonSimple
	{
	public:
		/// Wagon with sloshing
		/// this implementation uses the multi-mass model proposed by Graham and
		/// Rodriguez to simulate fluid sloshing in a  rectangular tank
		/// It uses tank length and height, combined to fluid density and a
		/// tank filling percentage to determine the values of the masses
		/// \brief implements a WagonBox with a sloshing model
		WagonSloshing(const std::string& name);
		/// Overload setter for mass
		/// \param mass
		void setTotalMass(double m_);
		/**
		 * \brief Initilize the wagon object
		 */
		void Initialize();

		double getN() const
		{
			return n;
		}

		void setN(double n_)
		{
			this->n = n_;
		}

		double getFillRatio() const
		{
			return fillRatio;
		}

		void setFillRatio(double fillRatio_)
		{
			this->fillRatio = fillRatio_;
		}

		MBSim::RigidBody* getSolidMass()
		{
			return solidMass;
		}
	protected:
		/**
		 * Value in ]0,1] to determine the tank filling ratio
		 */
		double fillRatio;
		/**
		 * The number of sloshing modes to be considered
		 */
		double n;
		/**
		 * Total mass of fluid
		 */
		double totalFluidMass;
		/**
		 * Fluid density
		 */
		double fluidRho;
		/**
		 * Tank aspect ratio
 		 */
		double aspectRatio;
		/**
		 * Vector of rigid bodies with the sloshing masses
		 */
		std::vector<MBSim::RigidBody*> massContainer;
		/**
		 * Vector with springs
		 */
		std::vector<MBSim::SpringDamper*> springs;
		/**
		 * Vector with translational joints to connect tank to masses
		 */
		std::vector<MBSim::Joint*> joints;
		/**
		 * Vector with frames that are connected to the springs
		 */
		std::vector<MBSim::FixedRelativeFrame*> frames;
		/**
		 * Solidified fluid mass
		 */
		MBSim::RigidBody* solidMass;
		MBSim::FixedRelativeFrame* solidMassConnectionFrame;
		fmatvec::SymMat solidMassInertia;
		/**
		 * Solid mass rigid connection  to tank
		 */
		MBSim::Joint* rigidJoint;
		/**
		 * Temporary variable
		 */
		double temp;
	};
#endif // WAGONMOVMASS_H
