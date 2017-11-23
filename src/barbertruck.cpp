/*
    Barber truck object
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


#include <barbertruck.h>
#include "mbsim/observers/contact_observer.h"
#include "mbsim/observers/mechanical_link_observer.h"
#include "mbsim/functions/kinetics/linear_elastic_function.h"
#include "mbsim/links/generalized_elastic_connection.h"

BarberTruck::BarberTruck( const std::string& projectName ): Truck(projectName) {
	BarberTruck (projectName,false);
};

BarberTruck::BarberTruck ( const std::string& projectName, bool withBushings ) : Truck(projectName), wheelBase(2.5), bolsterWithBushings(withBushings)
{
	setModelType("Barber");
	/// -------------------------- SYSTEM PARAMETERS --------------------------
	/// -------------------------------------------------------------------------

	double angleBolster = ( 36.3 ) *M_PI/180.; // [rad]
	double angleSideframe = ( 1.0 ) *M_PI/180.; // [rad]
	double bolsterWidth = 0.41; // taken at the widest part [m]
	double bolsterHeight = .19; // [m]
	double bolsterLength = 2.654; // [m]
	double bolsterMass = 342. + 6500; // [kg]
	fmatvec::SymMat bolsterInertiaTensor(3,fmatvec::EYE); // [kg.m²]
	bolsterInertiaTensor(0,0) = 97;
	bolsterInertiaTensor(1,1) = 99;
	bolsterInertiaTensor(2,2) = 7;
	double frictionCoefficient = 0.25; // [-]
	double wedgeMass = 2 * 12; // [kg] times two to represent 4 wedges per bolster
	double wedgeHeight = 0.2; // [m]
	double wedgeDepth = 0.10; // [m]
	double wedgeSpringStiffness = 2 * 240000; // [N/m]
	double wheelBase = 1.725;
	// double wheelRadius = 0.45;
	// side frame mass properties are doubled to represent both sides of the truck
	double sideFrameHeight = .495; // [m]
	double sideFrameWidth = .45; // [m]
	double sideFrameMass = 2 * 290; // [kg]
	double sideFrameTrack = 2.197; // [m]
	fmatvec::SymMat sideFrameInertiaTensor(3,fmatvec::EYE); // [kg.m²]
	sideFrameInertiaTensor(0,0) = 2 * 4;
	sideFrameInertiaTensor(1,1) = 2 * 91;
	sideFrameInertiaTensor(2,2) = 2 * 106;
	double truckTrack = 1.575; // [m]
	double bolsterSpringStiffness = 484000; // [N/m]
	double springBedOffset = 0.143; // [m]
	//   double t1,t2,t3 = 0; // temporary storage
	double coefRestitution = 0.05;  // TODO modify this parameter to be setted externally

	//acceleration of gravity
	MBSimEnvironment::getInstance()->setAccelerationOfGravity ( Vec ( "[0.;-9.810;0]" ) );

	/// ---------------------------- DEFINITION OF BODIES -----------------------
	/// -------------------------------------------------------------------------

	wedge2 = new Wedge ( "Wedge 2" );
	wedge1 = new Wedge ( "Wedge 1" );
	wedge3 = new Wedge ( "Wedge 3" );
	wedge4 = new Wedge ( "Wedge 4" );
	sideFrame = new RigidBody ( "Side frame" );
	bolster = new Bolster ( "Bolster" );
	wheelL = new Wheelset ( "Wheel Left ");
	wheelR = new Wheelset ( "Wheel Right" );

	this->addObject ( wedge2 );
	this->addObject ( wedge1 );
	this->addObject ( wedge3 );
	this->addObject ( wedge4 );
	this->addObject ( bolster );
	this->addObject ( sideFrame );
	this->addObject ( wheelL );
	this->addObject ( wheelR );

	/// ------------------------------ FRAMES -----------------------------------
	/// ------------ Frames on environment --------------------------------------
	Vec pos ( 3,INIT,0. );
	this->addFrame ( new FixedRelativeFrame ( "L", Vec ("[0;0;0]"), SqrMat (3,EYE)));
	// "L",Vec ( "[0.;0.;0.]" ),SqrMat ( 3,EYE ) ) );
	/// ------------ Frames for the wedges ----------------------------------------
	pos ( 0 ) = -.17102 + 0.1 * tan(angleSideframe);
	pos ( 1 ) = -0.0219;
	pos ( 2 ) = sideFrameTrack/2;
	bolster->addFrame ( new FixedRelativeFrame ( "W2",pos,SqrMat ( 3,EYE ) ) );
	pos ( 0 ) = -pos ( 0 );
	bolster->addFrame ( new FixedRelativeFrame ( "W1",pos,SqrMat ( 3,EYE ) ) );
	pos ( 2 ) = -sideFrameTrack/2;
	bolster->addFrame ( new FixedRelativeFrame ( "W3",pos,SqrMat ( 3,EYE ) ) );
	pos ( 0 ) = -pos ( 0 );
	bolster->addFrame ( new FixedRelativeFrame ( "W4",pos,SqrMat ( 3,EYE ) ) );
	/// ------------ Frames on the bolster --------------------------------------
	pos ( 0 ) = 0;
	pos ( 1 ) = -0.00785 + 0.1 + 0.00162;
	pos ( 2 ) = 0.0;
	this->addFrame ( new FixedRelativeFrame ( "B",pos,SqrMat ( 3,EYE ) ) );
	sideFrame->addFrame ( new FixedRelativeFrame ( "B",pos,SqrMat ( 3,EYE ) ) );
	sideFrame->addFrame ( new FixedRelativeFrame ( "BS",Vec ( "[0.0;-.18663;0.0]" ),SqrMat ( 3,EYE ),sideFrame->getFrame ( "B" ) ) );
	//
	pos(0) = -wheelBase/2;
	pos(1) = -0.0;
	this->addFrame(new FixedRelativeFrame("RL",pos,SqrMat(3,EYE)));
	sideFrame->addFrame(new FixedRelativeFrame("RL",pos,SqrMat(3,EYE)));
	pos(0) = wheelBase / 2;
	this->addFrame(new FixedRelativeFrame("RR",pos,SqrMat(3,EYE)));
	sideFrame->addFrame(new FixedRelativeFrame("RR",pos,SqrMat(3,EYE)));



	/// ---------------- DEFINITION OF THE WEDGES -------------------------------

	// Wedge 2
	wedge2->setMass ( wedgeMass );
	Vec wedgeAngles ( 2,INIT,0.0 );
	wedgeAngles ( 0 ) = angleBolster + M_PI;
	wedgeAngles ( 1 ) = angleSideframe;
	wedge2->setAngles ( wedgeAngles );
	wedge2->setHeight ( wedgeHeight );
	wedge2->setDepth ( wedgeDepth );
	wedge2->setInertiaTensor ( SymMat ( 3,EYE ) );
	wedge2->buildContour();
	wedge2->setFrameOfReference ( bolster->getFrame ( "W2" ) );
	wedge2->setFrameForKinematics ( wedge2->getFrame ( "C" ) );
	wedge2->getFrame ( "C" )->enableOpenMBV();
	wedge2->setTranslation ( new MBSim::LinearTranslation<VecV> (  ( "[1,0;0,1;0,0]" ) ) );
	wedge2->setRotation ( new RotationAboutZAxis<VecV>() );
	wedge2->enableOpenMBV ( true );

	// Wedge 4
	wedge4->setMass ( wedgeMass );
	wedge4->setAngles ( wedgeAngles );
	wedge4->setHeight ( wedgeHeight );
	wedge4->setDepth ( wedgeDepth );
	wedge4->setInertiaTensor ( SymMat ( 3,EYE ) );
	wedge4->buildContour();
	wedge4->setFrameOfReference ( bolster->getFrame ( "W4" ) );
	wedge4->setFrameForKinematics ( wedge4->getFrame ( "C" ) );
	wedge4->getFrame ( "C" )->enableOpenMBV();
	wedge4->setTranslation ( new MBSim::LinearTranslation<VecV> (  ( "[1,0;0,1;0,0]" ) ) );
	wedge4->setRotation ( new RotationAboutZAxis<VecV>() );
	wedge4->enableOpenMBV ( true );

	// Wedge 1
	wedge1->setMass ( wedgeMass );
	wedgeAngles ( 0 ) = angleSideframe;
	wedgeAngles ( 1 ) = angleBolster + M_PI;
	wedge1->setAngles ( wedgeAngles );
	wedge1->setHeight ( wedgeHeight );
	wedge1->setDepth ( wedgeDepth );
	wedge1->setInertiaTensor ( SymMat ( 3,EYE ) );
	wedge1->buildContour();
	wedge1->setFrameOfReference ( bolster->getFrame ( "W1" ) );
	wedge1->setFrameForKinematics ( wedge1->getFrame ( "C" ) );
	wedge1->getFrame ( "C" )->enableOpenMBV();
	wedge1->setTranslation ( new LinearTranslation<VecV> ( "[1,0;0,1;0,0]"  ) );
	wedge1->setRotation ( new RotationAboutZAxis<VecV>() );
	wedge1->enableOpenMBV ( true );

	// Wedge 3
	wedge3->setMass ( wedgeMass );
	wedge3->setAngles ( wedgeAngles );
	wedge3->setHeight ( wedgeHeight );
	wedge3->setDepth ( wedgeDepth );
	wedge3->setInertiaTensor ( SymMat ( 3,EYE ) );
	wedge3->buildContour();
	wedge3->setFrameOfReference ( bolster->getFrame ( "W3" ) );
	wedge3->setFrameForKinematics ( wedge3->getFrame ( "C" ) );
	wedge3->getFrame ( "C" )->enableOpenMBV();
	wedge3->setTranslation ( new LinearTranslation<VecV> ( "[1,0;0,1;0,0]"  ) );
	wedge3->setRotation ( new RotationAboutZAxis<VecV>() );
	wedge3->enableOpenMBV ( true );


	/// ---------------- DEFINITION OF THE BOLSTER ------------------------------
	// inertia properties
	bolster->setMass(bolsterMass);
	bolster->setInertiaTensor( bolsterInertiaTensor );
	// frames
	bolster->setFrameForKinematics(bolster->getFrameC());
	bolster->setFrameOfReference(getFrame("B"));
	// geometry
	bolster->setWidth(bolsterWidth);
	bolster->setHeight(bolsterHeight);
	bolster->setLength(bolsterLength);
	bolster->setAngle(angleBolster);
	bolster->setWagonConnectionPointPosition(0,bolsterHeight/2,bolsterLength/2);
	// update geometry reference frame
	bolster->setGeometryReferenceFramePosition(0.0,0.0,-bolsterLength/2);
	// display body
	bolster->enableOpenMBV(true);
	bolster->getFrameOfReference()->enableOpenMBV();
	// set the contact planes
	bolster->setContactPlanes();
	// degrees of freedom
	bolster->setTranslation ( new LinearTranslation<VecV> ( "[1,0,0;0,1,0;0,0,1]") );
	bolster->setRotation ( new RotationAboutAxesXZ<VecV>() );

	/// ------------------- DEFINITION OF THE SIDEFRAME -------------------------
	sideFrame->setMass ( sideFrameMass );
	sideFrame->setInertiaTensor( sideFrameInertiaTensor );
	sideFrame->setFrameOfReference(this->getFrameI());
	sideFrame->setFrameForKinematics ( sideFrame->getFrameC() );
	sideFrame->getFrameC()->enableOpenMBV();
	sideFrame->setTranslation( new LinearTranslation<VecV> ("[1,0;0,1;0,0]") );
	sideFrame->setRotation( new RotationAboutZAxis<VecV>() );

	// Contact plane sideframe-left wedge definition
	pos ( 0 ) = -sideFrameWidth/2;
	pos ( 1 ) = -sideFrameHeight/2;
	fmatvec::SqrMat3 rodaroda = BasicRotAIKz(-angleSideframe);
	sideFrame->addFrame(new FixedRelativeFrame("sframe_l_contact_plane",
			pos,
			rodaroda));
	Plane *sideframeL = new Plane ( "Side frame left",
			sideFrame->getFrame("sframe_l_contact_plane") );
	sideframeL->enableOpenMBV(_transparency=0.6);
	sideFrame->addContour ( sideframeL );

	// Contact plane sideframe-right wedge definition
	pos ( 0 ) = -pos ( 0 );
	sideFrame->addFrame(new FixedRelativeFrame("sframe_r_contact_plane",
			pos,
			BasicRotAKIz(-angleSideframe + M_PI)));
	Plane *sideframeR = new Plane ( "Side frame right",
			sideFrame->getFrame("sframe_r_contact_plane") );
	sideframeR->enableOpenMBV(_transparency=0.6);
	sideFrame->addContour ( sideframeR );


	pos ( 0 ) = 0.0;
	sideFrame->addFrame(new FixedRelativeFrame("sframe_central_contact_plane",
			pos,
			BasicRotAKIz(0.5 * M_PI)));
	Plane *ground = new Plane ( "Ground",
			sideFrame->getFrame("sframe_central_contact_plane") );
	ground->enableOpenMBV(_transparency=0.6);
	sideFrame->addContour ( ground );


	/// ---------------- DEFINITION OF THE WHELLS -------------------------------
	wheelL->setFrameOfReference(this->getFrame("RL"));
	wheelL->setFrameForKinematics(wheelL->getFrameC());
	wheelL->setMass(2000);
	wheelL->setTrack(truckTrack);
	wheelL->enableOpenMBV();

	wheelR->setFrameOfReference(this->getFrame("RR"));
	wheelR->setFrameForKinematics(wheelR->getFrameC());
	wheelR->setMass(2000);
	wheelR->setTrack(truckTrack);
	wheelR->enableOpenMBV();

	/// ---------------- DEFINITION OF JOINTS -----------------------------------
	//  Cylindrical joint left wheel-sideframe
	Joint* cylindrical1 = new Joint("CY1");
	cylindrical1->setForceDirection("[1,0;0,1;0,0]");
//	cylindrical1->setMomentDirection("[1;0;0]");
	cylindrical1->setForceLaw(new  BilateralConstraint());
//	cylindrical1->setMomentLaw(new BilateralConstraint());
	//cylindrical1->setImpactForceLaw(new BilateralImpact());
	cylindrical1->connect(sideFrame->getFrame("RL"),wheelL->getFrameC());
	addLink(cylindrical1);

	// Cylindrical joint right wheel-sideframe
	Joint* cylindrical2 = new Joint("CY2");
	cylindrical2->setForceDirection("[1,0;0,1;0,0]");
	cylindrical2->setForceLaw(new  BilateralConstraint());
	//cylindrical2->setImpactForceLaw(new BilateralImpact());
	cylindrical2->connect(sideFrame->getFrame("RR"),wheelR->getFrameC());
	addLink(cylindrical2);

	/// ---------------- DEFINITION OF THE SPRINGS
	///
	// Connection frames of the main spring packs
	setSpringConnectionPointsFrames(bolster,
			springBedOffset, sideFrameTrack/2,bolsterHeight/2,1);
	setSpringConnectionPointsFrames(bolster,
			springBedOffset, -sideFrameTrack/2,bolsterHeight/2,2);
	setSpringConnectionPointsFrames(sideFrame,
			springBedOffset, sideFrameTrack/2,sideFrameHeight/2,1);
	setSpringConnectionPointsFrames(sideFrame,
			springBedOffset, -sideFrameTrack/2,sideFrameHeight/2,2);

	// Connection frames of the springs
	Vec3 relativePosition(INIT,0.);
	relativePosition(2) = sideFrameTrack/2;
	relativePosition(1) = -sideFrameHeight/2;
	relativePosition(0) = springBedOffset;

	std::string wedgeSpringNameRoot("Wedge_spring_");

	sideFrame->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + toStr(1),
			relativePosition,SqrMat(3,EYE)));
	relativePosition(0) = - relativePosition(0);
	sideFrame->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + toStr(2),
				relativePosition,SqrMat(3,EYE)));
	relativePosition(2) = - relativePosition(2);
	sideFrame->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + toStr(4),
				relativePosition,SqrMat(3,EYE)));
	relativePosition(0) = - relativePosition(0);
	sideFrame->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + toStr(3),
					relativePosition,SqrMat(3,EYE)));

	relativePosition = 0. * relativePosition;
	relativePosition(1) = -wedgeHeight/3;
	wedge1->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot,
			relativePosition,SqrMat(3,EYE)));
	wedge2->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot,
				relativePosition,SqrMat(3,EYE)));
	wedge3->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot,
					relativePosition,SqrMat(3,EYE)));
	wedge4->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot,
						relativePosition,SqrMat(3,EYE)));


	// Wedge spring group force law
	LinearSpringDamperForce *wedgeSpringLaw =
			new LinearSpringDamperForce(wedgeSpringStiffness,0.002 * wedgeSpringStiffness);

	// Bolster spring group force law
	LinearSpringDamperForce *bolsterSpringLaw =
			new LinearSpringDamperForce(bolsterSpringStiffness,0.002 * bolsterSpringStiffness);

	// Spring group connecting right wedge to sideframe
	for (unsigned int i = 1; i <= 4; i++){
		SpringDamper *springWedge = new SpringDamper ( "Spring-Wedge-" + toStr(i) );
		springWedge->setUnloadedLength(0.27);
		springWedge->connect ( sideFrame->getFrame ( wedgeSpringNameRoot + toStr(i) ),
				dynamic_cast<RigidBody*>(getObject("Wedge " + toStr(i)))->getFrame ( wedgeSpringNameRoot) );
		springWedge->setForceFunction ( wedgeSpringLaw );
		this->addLink ( springWedge );
		springWedge->enableOpenMBV(_springRadius=0.8 * springBedOffset / 2,
				_crossSectionRadius=0.008,
				_numberOfCoils=5);
	}

	// Spring group connecting bolster to sideframe
	std::stringstream springName;
	std::stringstream frameName;
	for (unsigned k=0; k < 2; k++){
		for (unsigned i=0; i < 3; i++){
			for (unsigned j=0; j < 3; j++){
				springName.str("");
				frameName.str("");
				springName << "Spring-Bolster_" << i+1 << j+1 << "_" << static_cast<char> (k+65);
				frameName << "_SpringBed_" << i+1 << j+1 << "_" << static_cast<char> (65 + k);
				SpringDamper *springBolster = new SpringDamper (springName.str());
				springBolster->setForceFunction(bolsterSpringLaw);
				springBolster->setUnloadedLength(0.254);
				springBolster->enableOpenMBV(_springRadius=0.8 * springBedOffset / 2,
						_crossSectionRadius=0.010,
						_numberOfCoils=5);
				springBolster->connect(
						sideFrame->getFrame(sideFrame->getName() + frameName.str()),
						bolster->getFrame(bolster->getName() + frameName.str())
				);
				// barber trucks doesn't have bolster springs under the wedges
				// therefore, springs on position 02 and 12 doesn't exist
				if ( j == 0 || j == 1) this->addLink(springBolster);
				else if ( i == 2) this->addLink(springBolster);
				springBolster->setPlotFeature("generalizedForce", enabled);
				springBolster->setPlotFeature("deflection",enabled);

				// IF ELASTIC CONNECTIONS

				if (bolsterWithBushings){
					LinearElasticFunction *stiffnessFcn = new LinearElasticFunction();
					SymMat3 stiffnessMatrix(INIT,0.0);
					stiffnessMatrix(0,0) = 40900;  	//stiffnessMatrix(0,2) = 23300;
					stiffnessMatrix(1,0) = -5250;	//stiffnessMatrix(1,2) = -5250;
					stiffnessMatrix(2,0) = 23300;	//stiffnessMatrix(2,2) = 40900;

					stiffnessFcn->setStiffnessMatrix(stiffnessMatrix);
					stiffnessFcn->setDampingMatrix(stiffnessMatrix * 0.004);

					springName << "-Bushing";
					ElasticJoint *bushingBolster =
							new ElasticJoint (springName.str());
					bushingBolster->setGeneralizedForceFunction(stiffnessFcn);
					bushingBolster->setForceDirection("[1,0,0;0,1,0;0,0,1");

					/// The elastic connection element doesn't have an unloaded
					/// length definition. Therefore, because on the spring definition
					/// the frames are offset from each other, there is always a
					/// vertical relative displacement, which causes the force to be
					/// non-zero even if there is no radial deflection. To cope with
					/// this issue, a dummy frame is created in the bolster, but
					/// on the same position that the connection point of the
					/// sideframe is, leading to zero deflection.
					FixedRelativeFrame *dummyFrame = new FixedRelativeFrame(
							bolster->getName() + frameName.str() + "-2",
							Vec3("[0;-0.24627;0]"),SqrMat3(EYE),
							bolster->getFrame(bolster->getName() + frameName.str()));
					bolster->addFrame(dummyFrame);
					bushingBolster->connect(
							sideFrame->getFrame(sideFrame->getName() + frameName.str()),
							dummyFrame);
					if ( j == 0 || j == 1) this->addLink(bushingBolster);
					else if ( i == 2) this->addLink(bushingBolster);
					bushingBolster->setPlotFeatureRecursive("generalizedForce", enabled);
					bushingBolster->setPlotFeatureRecursive("generalizedRelativePosition",enabled);
					bushingBolster->setPlotFeatureRecursive("generalizedRelativeVelocity",enabled);
				}
			}
		}
	}


	/// ------------------------ DEFITION OF THE CONTACTS -----------------------

	// Establish contacts

	// contacts between left wedge and sideframe
	Contact *contact;
	int numberOfContacts = wedge2->getContours().size() / 2;
	std::vector<int> idx(4,0);
	idx[0] = 0;
	idx[1] = 1;
	idx[2] = 4;
	idx[3] = 5;

	for ( int i = 0; i < numberOfContacts; i++ )
	{
		stringstream s;
		string cNumber;
		s << i+1;
		cNumber = s.str();

		contact = new Contact ( std::string ( "Contact_WedgeSideframe_Left-" ) + cNumber );
		contact->connect ( wedge2->getContours() [idx[i]],sideframeL );
		contact->setNormalForceLaw ( new UnilateralConstraint() );
		contact->setNormalImpactLaw ( new UnilateralNewtonImpact ( coefRestitution ) );
		contact->setTangentialImpactLaw ( new PlanarCoulombImpact ( frictionCoefficient ) );
		contact->setTangentialForceLaw ( new PlanarCoulombFriction ( frictionCoefficient ) );
		this->addLink ( contact );
		contact->setPlotFeature("generalizedForce",enabled);
		contact->setPlotFeature("generalizedRelativePosition",enabled);
	}

	// contacts between left wedge and bolster
	idx[0] = 2;
	idx[1] = 3;
	idx[2] = 6;
	idx[3] = 7;

	for ( int i = 0; i < numberOfContacts; i++ )
	{
		stringstream s;
		string cNumber;
		s << i+1;
		cNumber = s.str();

		contact = new Contact ( std::string ( "Contact_WedgeBolster_Left-" ) + cNumber );
		contact->connect ( wedge2->getContours() [idx[i]],bolster->getContour("Contact plane left") );
		contact->setNormalForceLaw ( new UnilateralConstraint() );
		contact->setNormalImpactLaw ( new UnilateralNewtonImpact ( coefRestitution ) );
		contact->setTangentialImpactLaw ( new PlanarCoulombImpact ( frictionCoefficient ) );
		contact->setTangentialForceLaw ( new PlanarCoulombFriction ( frictionCoefficient ) );
		this->addLink ( contact );
		contact->setPlotFeature("generalizedForce",enabled);
		contact->setPlotFeature("generalizedRelativePosition",enabled);

	}

	// contacts between right wedge and bolster
	idx[0] = 0;
	idx[1] = 1;
	idx[2] = 4;
	idx[3] = 5;

	for ( int i = 0; i < numberOfContacts; i++ )
	{
		stringstream s;
		string cNumber;
		s << i+1;
		cNumber = s.str();

		contact= new Contact ( std::string ( "Contact_WedgeBolster_Right-" ) + cNumber );
		contact->connect ( wedge1->getContours() [idx[i]],bolster->getContour("Contact plane right") );
		contact->setNormalForceLaw ( new UnilateralConstraint() );
		contact->setNormalImpactLaw ( new UnilateralNewtonImpact ( coefRestitution ) );
		contact->setTangentialImpactLaw ( new PlanarCoulombImpact ( frictionCoefficient ) );
		contact->setTangentialForceLaw ( new PlanarCoulombFriction ( frictionCoefficient ) );
		this->addLink ( contact);
		contact->setPlotFeature("generalizedForce",enabled);
		contact->setPlotFeature("generalizedRelativePosition",enabled);
	}

	// contacts between right wedge and sideframe
	idx[0] = 2;
	idx[1] = 3;
	idx[2] = 6;
	idx[3] = 7;

	for ( int i = 0; i < numberOfContacts; i++ )
	{
		stringstream s;
		string cNumber;
		s << i+1;
		cNumber = s.str();

		contact = new Contact ( std::string ( "Contact_WedgeSideframe_Right-" ) + cNumber );
		contact->connect ( wedge1->getContours() [idx[i]],sideframeR );
		contact->setNormalForceLaw ( new UnilateralConstraint() );
		contact->setNormalImpactLaw ( new UnilateralNewtonImpact ( coefRestitution ) );
		contact->setTangentialImpactLaw ( new PlanarCoulombImpact ( frictionCoefficient ) );
		contact->setTangentialForceLaw ( new PlanarCoulombFriction ( frictionCoefficient ) );
		this->addLink ( contact );
		contact->setPlotFeature("generalizedForce",enabled);
		contact->setPlotFeature("generalizedRelativePosition",enabled);
	}
}

// getWheelBase
double BarberTruck::getWheelBase()
{
	return wheelBase;
}

// setWheelBase
void BarberTruck::setWheelBase(double wB_)
{
	wheelBase = wB_;
}

void BarberTruck::setSpringConnectionPointsFrames(MBSim::RigidBody *_body,
		double _distanceBetweenSprings,
		double _lateralOffset,
		double _height,
		int lado)
{

	fmatvec::Vec3 offsetTable(fmatvec::INIT,0);
	offsetTable(0) = -1;
	offsetTable(1) = 1;


	fmatvec::Vec3 relativePosition(fmatvec::INIT,0);

	std::stringstream frameName;

	relativePosition(2) = _lateralOffset;
	relativePosition(1) = - _height;

	frameName.str("");
	frameName << "Frame_BolsterSpringCenter_" << static_cast<char> (64 + lado);
	_body->addFrame( new FixedRelativeFrame(frameName.str(),
			relativePosition,
			fmatvec::SqrMat(3,fmatvec::EYE)));
	Frame *refFrame;
	refFrame = _body->getFrame(frameName.str());
	refFrame->enableOpenMBV(_size=0.3);

	for ( unsigned i = 0; i < 3; i++) {
		for ( unsigned j = 0; j < 3; j++){
			relativePosition.init(0);
			frameName.str("");
			frameName << _body->getName() << "_SpringBed_" << i+1 << j+1 << "_" << static_cast<char> (64 + lado);
			relativePosition(0) = _distanceBetweenSprings * offsetTable(i);
			relativePosition(2) = _distanceBetweenSprings * offsetTable(j);

			_body->addFrame(new FixedRelativeFrame(frameName.str(),
					relativePosition,
					fmatvec::SqrMat(3, fmatvec::EYE),
					refFrame));
		}
	}

}



