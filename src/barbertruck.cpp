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


#include "barbertruck.h"

BarberTruck::BarberTruck( const std::string& projectName ): Truck(projectName) {
	BarberTruck (projectName,false, 2.0);
};

BarberTruck::BarberTruck ( const std::string& projectName, bool withBushings, double _wBase ) :
		Truck(projectName), wheelBase(_wBase), bolsterWithBushings(withBushings)
{
	setModelType("Barber");
	/// -------------------------- SYSTEM PARAMETERS --------------------------
	/// -------------------------------------------------------------------------

	/// wedge properties
	double angleBolster = ( 36.3 ) *M_PI/180.; // [rad]
	double angleSideframe = ( 1.0 ) *M_PI/180.; // [rad]
	double frictionCoefficient = 0.01; // [-]
	double wedgeMass = 12; // [kg] times two to represent 4 wedges per bolster
	double wedgeHeight = 0.2; // [m]
	double wedgeDepth = 0.12; // [m]
	double wedgeSpringVerticalStiffness = 393046; // [N/m]
	double wedgeSpringShearStiffness = 118870;  // [N/m]
	double wedgeSpringFreeLength = 0.28; // [m]
	/// bolster properties
	double bolsterWidth = 0.41; // taken at the widest part [m]
	double bolsterHeight = .19; // [m]
	double bolsterLength = 2.654; // [m]
	double bolsterMass = 911.; // [kg]
	fmatvec::SymMat bolsterInertiaTensor(3,fmatvec::EYE); // [kg.m²]
	bolsterInertiaTensor(0,0) = 508.90;
	bolsterInertiaTensor(1,1) = 511.40;
	bolsterInertiaTensor(2,2) = 24.80;
	double bolsterSpringStiffness = 0*572410.86; // rideControl= 4.839868e+05; // [N/m]
	double bolsterSpringDampingFactor = 0.004;
	double bolsterSpringFreeLength = 0.260;
	double springBedOffset = 0.143; // [m]
	double wheelBase = 1.829;
	double wheelsetMass = 1777.55; // [kg]
	fmatvec::SymMat wheelsetInertiaTensor(3,fmatvec::EYE); // [kg.m²]
	wheelsetInertiaTensor(0,0) = 1179.60;
	wheelsetInertiaTensor(1,1) = 1179.60;
	wheelsetInertiaTensor(2,2) = 138.00;
	// double wheelRadius = 0.9652;
	// side frame mass properties are doubled to represent both sides of the truck
	double sideFrameHeight = .495; // [m]
	double sideFrameWidth = .45; // [m]
	double sideFrameMass = 630.00; // [kg]
	double sideFrameTrack = 2 * 1.098; // [m]
	fmatvec::SymMat sideFrameInertiaTensor(3,fmatvec::EYE); // [kg.m²]
	sideFrameInertiaTensor(0,0) = 91.60;
	sideFrameInertiaTensor(1,1) = 323.00;
	sideFrameInertiaTensor(2,2) = 340.30;
	double truckTrack = 1.575; // [m]
	//   double t1,t2,t3 = 0; // temporary storage
	double coefRestitution = 0.00;  // TODO modify this parameter to be setted externally

	//acceleration of gravity
	//MBSimEnvironment::getInstance()->setAccelerationOfGravity ( Vec ( "[0.;-9.810;0]" ) );

	/// ---------------------------- DEFINITION OF BODIES -----------------------
	/// -------------------------------------------------------------------------

	wedge2 = new Wedge ( "Wedge 2" );
	wedge1 = new Wedge ( "Wedge 1" );
	wedge3 = new Wedge ( "Wedge 3" );
	wedge4 = new Wedge ( "Wedge 4" );
	sideFrameLeft = new Sideframe ( "Sideframe Left" );
	sideFrameRight = new Sideframe ( "Sideframe Right" );
	bolster = new Bolster ( "Bolster" );
	wheelRear = new Wheelset ( "Wheelset Rear", truckTrack, sideFrameTrack);
	wheelFront = new Wheelset ( "Wheelset Front", truckTrack, sideFrameTrack);

	this->addObject ( wedge2 );
	this->addObject ( wedge1 );
	this->addObject ( wedge3 );
	this->addObject ( wedge4 );
	this->addObject ( bolster );
	this->addObject ( sideFrameLeft );
	this->addObject( sideFrameRight );
	this->addObject ( wheelRear );
	this->addObject ( wheelFront );

	/// ------------------------------ FRAMES -----------------------------------
	/// ------------ Frames on environment --------------------------------------
	Vec pos ( 3,INIT,0. );
	this->addFrame ( new FixedRelativeFrame ( "L", Vec ("[0;0;0]"), SqrMat (3,EYE)));
	// "L",Vec ( "[0.;0.;0.]" ),SqrMat ( 3,EYE ) ) );
	/// ------------ Frames for the wedges ----------------------------------------
	pos ( 0 ) = -.17102 + 0.1 * tan(angleSideframe);
	pos ( 1 ) = -0.02191;
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
	sideFrameLeft->addFrame ( new FixedRelativeFrame ( "B",pos,SqrMat ( 3,EYE ) ) );
	sideFrameLeft->addFrame ( new FixedRelativeFrame ( "BS",Vec ( "[0.0;-.18663;0.0]" ),SqrMat ( 3,EYE ),sideFrameLeft->getFrame ( "B" ) ) );
	//
	pos(0) = -wheelBase/2;
	pos(1) = -0.0;
	pos(2) =0;
	this->addFrame(new FixedRelativeFrame("RL",pos,SqrMat(3,EYE)));
	sideFrameLeft->addFrame(new FixedRelativeFrame("RL",pos,SqrMat(3,EYE)));
	sideFrameRight->addFrame(new FixedRelativeFrame("RL",pos,SqrMat(3,EYE)));
	pos(0) = wheelBase / 2;
	this->addFrame(new FixedRelativeFrame("RR",pos,SqrMat(3,EYE)));
	sideFrameLeft->addFrame(new FixedRelativeFrame("RR",pos,SqrMat(3,EYE)));
	sideFrameRight->addFrame(new FixedRelativeFrame("RR",pos,SqrMat(3,EYE)));

	sideFrameRight->getFrame("RR")->enableOpenMBV();
	sideFrameLeft->getFrame("RR")->enableOpenMBV();
	/// ------------- Sideframes Frames ----------------------------------------
	pos.init(0);
	pos(2) = -sideFrameTrack / 2;
	this->addFrame(new FixedRelativeFrame("SF1",pos, SqrMat(3,EYE)));
	pos(2) = - pos(2);
	this->addFrame(new FixedRelativeFrame("SF2",pos, SqrMat(3,EYE)));




	/// ---------------- DEFINITION OF THE WEDGES -------------------------------

	VecV translationDirection(3,INIT,0.0);
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
	bolster->setRotation ( new RotationAboutAxesXYZ<VecV>() );

	/// ------------------- DEFINITION OF THE SIDEFRAMES ------------------------
	sideFrameLeft->setMass ( sideFrameMass );
	sideFrameLeft->setInertiaTensor( sideFrameInertiaTensor );
	sideFrameLeft->setFrameOfReference(this->getFrame("SF1"));
	sideFrameLeft->setFrameForKinematics ( sideFrameLeft->getFrameC() );
	sideFrameLeft->getFrameC()->enableOpenMBV();
	sideFrameLeft->setTranslation( new LinearTranslation<VecV> ("[1,0,0;0,1,0;0,0,1]") );
	sideFrameLeft->setRotation( new RotationAboutAxesXYZ<VecV>() );

	sideFrameLeft->setWindowAngle(angleSideframe);
	sideFrameLeft->setWindowHeight(sideFrameHeight);
	sideFrameLeft->setWindowWidth(sideFrameWidth);
	sideFrameLeft->buildContour();

	sideFrameRight->setMass ( sideFrameMass );
	sideFrameRight->setInertiaTensor( sideFrameInertiaTensor );
	sideFrameRight->setFrameOfReference(this->getFrame("SF2"));
	sideFrameRight->setFrameForKinematics ( sideFrameRight->getFrameC() );
	sideFrameRight->getFrameC()->enableOpenMBV();
	sideFrameRight->setTranslation( new LinearTranslation<VecV> ("[1,0,0;0,1,0;0,0,1]") );
	sideFrameRight->setRotation( new RotationAboutAxesXYZ<VecV>() );

	sideFrameRight->setWindowAngle(angleSideframe);
	sideFrameRight->setWindowHeight(sideFrameHeight);
	sideFrameRight->setWindowWidth(sideFrameWidth);
	sideFrameRight->buildContour();


	/// ---------------- DEFINITION OF THE WHEELS -------------------------------
	Vec3 wheelsetOffset(INIT,0.0);
	wheelsetOffset(0) = - wheelBase / 2;
	this->addFrame(
			new FixedRelativeFrame("RWS",wheelsetOffset,SqrMat(3,EYE))
	);
	this->addFrame(
			new FixedRelativeFrame("FWS",-wheelsetOffset,SqrMat(3,EYE))
	);


	wheelRear->setMass(wheelsetMass);
	wheelRear->setInertiaTensor(wheelsetInertiaTensor);
	wheelRear->setFrameOfReference(this->getFrame("RWS"));
	wheelRear->setFrameForKinematics(wheelRear->getFrameC());

	wheelFront->setMass(wheelsetMass);
	wheelFront->setInertiaTensor(wheelsetInertiaTensor);
	wheelFront->setFrameOfReference(this->getFrame("FWS"));
	wheelFront->setFrameForKinematics(wheelFront->getFrameC());

	/// ---------------- DEFINITION OF JOINTS -----------------------------------
	///
	///
	SymMatV wheelBoxStiffMatrix(6,INIT,0.0);
	wheelBoxStiffMatrix(0,0) = 13134500; // longitudinal
	wheelBoxStiffMatrix(2,2) = 8756300; // lateral
	wheelBoxStiffMatrix(1,1) = 175127000; // vertical
	wheelBoxStiffMatrix(4,4) = 55475; // yaw
	wheelBoxStiffMatrix(3,3) = 55475; // roll

	//  Wheelbox joint rear wheelset-sideframe left/
	addLink(new WheelBox("Wheelbox: left sideframe to rear wheelset",
			sideFrameLeft->getFrame("RL"),
			wheelRear->getFrame("SFL"), wheelBoxStiffMatrix));
	// Wheelbox joint front wheelset-sideframe left
	addLink(new WheelBox("Wheelbox: left sideframe to front wheelset",
			sideFrameLeft->getFrame("RR"),wheelFront->getFrame("SFL"),
			wheelBoxStiffMatrix));

	//  Wheelbox joint rear wheelset-sideframe right
	addLink(new WheelBox("Wheelbox: right sideframe to rear wheelset",
			sideFrameRight->getFrame("RL"),wheelRear->getFrame("SFR"),
			wheelBoxStiffMatrix));

	// Wheelbox joint front wheelset-sideframe right
	addLink(new WheelBox("Wheelbox: right sideframe to front wheelset",
			sideFrameRight->getFrame("RR"),wheelFront->getFrame("SFR"),
			wheelBoxStiffMatrix));

	/// ---------------- DEFINITION OF THE SPRINGS
	///
	// Connection frames of the main spring packs
	setSpringConnectionPointsFrames(bolster,
			springBedOffset, sideFrameTrack/2,bolsterHeight/2,1);
	setSpringConnectionPointsFrames(bolster,
			springBedOffset, -sideFrameTrack/2,bolsterHeight/2,2);
	setSpringConnectionPointsFrames(sideFrameRight,
			springBedOffset, 0*sideFrameTrack,sideFrameHeight/2,1);
	setSpringConnectionPointsFrames(sideFrameLeft,
			springBedOffset, 0,sideFrameHeight/2,2);

	// Connection frames of the springs
	Vec3 relativePosition(INIT,0.);
	relativePosition(2) = sideFrameTrack;
	relativePosition(1) = -sideFrameHeight/2;
	relativePosition(0) = springBedOffset;

	std::string wedgeSpringNameRoot("Wedge_spring_");

	sideFrameLeft->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + to_string(1),
			relativePosition,SqrMat(3,EYE)));
	relativePosition(0) = - relativePosition(0);
	sideFrameLeft->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + to_string(2),
			relativePosition,SqrMat(3,EYE)));
	relativePosition(2) = 0;
	sideFrameLeft->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + to_string(4),
			relativePosition,SqrMat(3,EYE)));
	relativePosition(0) = - relativePosition(0);
	sideFrameLeft->addFrame(new FixedRelativeFrame(wedgeSpringNameRoot + to_string(3),
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
	SymMat3 wedgeSpringStiffness(INIT,0.0);
	wedgeSpringStiffness(0,0) = wedgeSpringShearStiffness;
	wedgeSpringStiffness(1,1) = wedgeSpringVerticalStiffness;
	wedgeSpringStiffness(2,2) = wedgeSpringShearStiffness;
	LinearElasticFunctionWithClearances *wedgeSpringLaw = new LinearElasticFunctionWithClearances();
	wedgeSpringLaw->setStiffnessMatrix(wedgeSpringStiffness);
	wedgeSpringLaw->setDampingMatrix(wedgeSpringStiffness * bolsterSpringDampingFactor);
	VecV wedgeSpringClearances(3,INIT,0);
	wedgeSpringClearances(1) = wedgeSpringFreeLength;
	wedgeSpringLaw->setClearances(wedgeSpringClearances);

	// Bolster spring group force law
	LinearSpringDamperForce *bolsterSpringLaw =
			new LinearSpringDamperForce(bolsterSpringStiffness,0.002 * bolsterSpringStiffness);

	// Spring group connecting right wedge to sideframe
	for (unsigned int i = 1; i <= 4; i++){
		ElasticJoint *springWedge = new ElasticJoint ( "Spring-Wedge-" + to_string(i) );
		springWedge->connect ( sideFrameLeft->getFrame ( wedgeSpringNameRoot + to_string(i) ),
				dynamic_cast<RigidBody*>(getObject("Wedge " + to_string(i)))->getFrame ( wedgeSpringNameRoot) );
		springWedge->setForceDirection("[1,0,0;0,1,0;0,0,1]");
		springWedge->setGeneralizedForceFunction(wedgeSpringLaw);
		this->addLink ( springWedge );
//		springWedge->enableOpenMBV(_springRadius=0.8 * springBedOffset / 2,
//				_crossSectionRadius=0.008,
//				_numberOfCoils=5);
	}

	// Spring group connecting bolster to sideframe
	std::stringstream springName;
	std::stringstream frameName;
	for (unsigned k=0; k < 2; k++){ // this loops alternates between the sideframes
		for (unsigned i=0; i < 3; i++){
			for (unsigned j=0; j < 3; j++){
				springName.str("");
				frameName.str("");
				springName << "Spring-Bolster_" << i+1 << j+1 << "_" << static_cast<char> (k+65);
				frameName << "_SpringBed_" << i+1 << j+1 << "_" << static_cast<char> (65 + k);
				SpringDamper *springBolster = new SpringDamper (springName.str());
				springBolster->setForceFunction(bolsterSpringLaw);
				springBolster->setUnloadedLength(bolsterSpringFreeLength);
				springBolster->enableOpenMBV(_springRadius=0.8 * springBedOffset / 2,
						_crossSectionRadius=0.010,
						_numberOfCoils=5);
				if ( k == 1 ) springBolster->connect(
						sideFrameLeft->getFrame(sideFrameLeft->getName() + frameName.str()),
						bolster->getFrame(bolster->getName() + frameName.str()));
				else if ( k == 0 ) springBolster->connect(
						sideFrameRight->getFrame(sideFrameRight->getName() + frameName.str()),
						bolster->getFrame(bolster->getName() + frameName.str()));
				// therefore, springs on position 02 and 12 don't exist
				if ( j == 0 || j == 1) this->addLink(springBolster);
				else if ( i == 2) this->addLink(springBolster);
				springBolster->setPlotFeature(generalizedForce, true);
				springBolster->setPlotFeature(deflection,true);

				// IF ELASTIC CONNECTIONS

				if (bolsterWithBushings){
					/* Molas RideControl
4.085513e+04	2.535569e+04	1.818989e-12	1.607263e+04	7.427344e+05	2.014681e+07
2.535569e+04	4.839868e+05	-5.261543e+03	-1.011572e+05	5.402705e+05	3.508579e+06
1.818989e-12	-5.261543e+03	4.085513e+04	-2.014681e+07	2.423125e+05	1.607263e+04
1.607263e+04	-1.011572e+05	-2.014681e+07	3.278376e+09	1.070903e+08	-2.980232e-08
7.427344e+05	5.402705e+05	2.423125e+05	1.070903e+08	1.692246e+09	4.999814e+07
2.014681e+07	3.508579e+06	1.607263e+04	-2.980232e-08	4.999814e+07	3.278376e+09

					 * */

					/* Molas Motion Control Ride Master x e z são radiais e y é vertical
					 *
					 *208548.3305736644	-2766.1381645451406	-3521.3778699214	560.3736702170961	-856.761867682231	1123.0811471920993
					 *-2766.1381645451406	572410.85895094	-2766.1381645451406	-1733.25931022422	-1392.11431652101	-1733.25931022422
					 *-3521.3778699214	-2766.1381645451406	208548.3305736644	1123.0811471920993	-856.761867682231	560.3736702170961
					 *560.3736702170961	-1733.25931022422	1123.0811471920993	4717.90900078866	13.63578809905135	42.544610441601094
					 *-856.761867682231	-1392.11431652101	-856.761867682231	13.63578809905135	1994.4260229744702	13.63578809905135
					 *1123.0811471920993	-1733.25931022422	560.3736702170961	42.544610441601094	13.63578809905135	4717.90900078866
					 *
					 *
					 */
					LinearElasticFunctionWithClearances *stiffnessFcn = new LinearElasticFunctionWithClearances();
					SymMat stiffnessMatrix(6,INIT,0.0);
					stiffnessMatrix(0,0) = 208548.3305;
					stiffnessMatrix(1,0) = -2766.1381;	stiffnessMatrix(1,1) = 572410.8589;
					stiffnessMatrix(2,0) = -3521.3778;	stiffnessMatrix(2,1) = -2766.1381;	stiffnessMatrix(2,2) = 208548.3305;
					stiffnessMatrix(3,0) = 560.3736; 	stiffnessMatrix(3,1) = -1733.2593;	stiffnessMatrix(3,2) = 1123.0811;
						stiffnessMatrix(3,3) = 4717.9090;
					stiffnessMatrix(4,0) = -856.7618;	stiffnessMatrix(4,1) = -1392.1143;	stiffnessMatrix(4,2) = -856.7618;
						stiffnessMatrix(4,3) = 13.6357; stiffnessMatrix(4,4) = 1994.4260;
					stiffnessMatrix(5,0) = 1123.0811;	stiffnessMatrix(5,1) = -1733.2593;	stiffnessMatrix(5,2) = 560.3736;
						stiffnessMatrix(5,3) = 42.5446;			stiffnessMatrix(5,4) = 13.6357;	stiffnessMatrix(5,5) = 4717.9090;

					// rotate the sitffness matrix phi around the axial axis
					// each side frame has a different phi, so that forces tend to cancel each other on the abscence
					// of external excitation.
					double phi;
					if (k == 0 ) phi = 0;
					if (k == 1 ) phi = M_PI;
					SqrMat3 rotateY = BasicRotAIKy(phi);
					SqrMat randomRotationMatrix (6,EYE);
					for (unsigned ri=0;ri<=2;ri++){
						for (unsigned rj=0;rj<=2;rj++){
							randomRotationMatrix(ri,rj) = rotateY(ri,rj);
							randomRotationMatrix(ri+3,rj+3) = rotateY(ri,rj);
						}
					}


					SymMat rotatedStiffness(6,INIT,0.0);
					rotatedStiffness = fmatvec::JTMJ(stiffnessMatrix,randomRotationMatrix);

					// sets the bushings with the rotated stiffness matrix
					stiffnessFcn->setStiffnessMatrix(rotatedStiffness);
					stiffnessFcn->setDampingMatrix(rotatedStiffness * bolsterSpringDampingFactor);
					VecV freeLengths(6,INIT,0);
					freeLengths(1) = -bolsterSpringFreeLength;
					stiffnessFcn->setClearances(freeLengths);



					springName << "-Bushing";
					ElasticJoint *bushingBolster =
							new ElasticJoint (springName.str());
					bushingBolster->setGeneralizedForceFunction(stiffnessFcn);
					bushingBolster->setForceDirection("[1,0,0;0,1,0;0,0,1]");
					bushingBolster->setMomentDirection("[1,0,0;0,1,0;0,0,1]");


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
					if ( k == 1 ) bushingBolster->connect(bolster->getFrame(bolster->getName() + frameName.str()),
							sideFrameLeft->getFrame(sideFrameLeft->getName() + frameName.str()));
//							dummyFrame);
					else if ( k == 0 ) bushingBolster->connect(bolster->getFrame(bolster->getName() + frameName.str()),
							sideFrameRight->getFrame(sideFrameRight->getName() + frameName.str()));
//							dummyFrame);
					if ( j == 0 || j == 1) this->addLink(bushingBolster);
					else if ( i == 2) this->addLink(bushingBolster);
					bushingBolster->setPlotFeatureRecursive(generalizedForce, true);
					bushingBolster->setPlotFeatureRecursive(generalizedRelativePosition,true);
					bushingBolster->setPlotFeatureRecursive(generalizedRelativeVelocity,true);
				}
			}
		}
	}
// END OF SPRING CONNECTIONS

	/// ------------------------ DEFITION OF THE CONTACTS -----------------------

	// Establish contacts
	bool contactsSideframe=false;
	bool contactsBolster=false;


	if(contactsSideframe){
		setWedgeContacts(wedge1->getContour("Right face"),sideFrameRight->getContour("Side frame right"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge2->getContour("Left face"),sideFrameRight->getContour("Side frame left"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge3->getContour("Right face"),sideFrameLeft->getContour("Side frame right"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge4->getContour("Left face"),sideFrameLeft->getContour("Side frame left"),
				frictionCoefficient,coefRestitution);
	}
	else{
		/// If contact is not active, planar joints are set
		Joint *wedgeJointS = new Joint("Wedge 1 to sideframe");
		wedgeJointS->connect(wedge1->getFrame("Right face central frame"),sideFrameRight->getFrame("sframe_r_contact_plane"));
		wedgeJointS->setForceLaw(new MBSim::BilateralConstraint);
		// wedgeJointS->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointS->setForceDirection("[1;0;0]");
		// wedgeJointS->setMomentDirection("[0,0;1,0;0,1]");
		wedgeJointS->setPlotFeature(force,true);
		wedgeJointS->setPlotFeature(generalizedRelativePosition,true);
		addLink(wedgeJointS);

		wedgeJointS = new Joint("Wedge 2 to sideframe");
		wedgeJointS->connect(wedge2->getFrame("Left face central frame"),sideFrameRight->getFrame("sframe_l_contact_plane"));
		wedgeJointS->setForceLaw(new MBSim::BilateralConstraint);
		// wedgeJointS->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointS->setForceDirection("[1;0;0]");
		// wedgeJointS->setMomentDirection("[0,0;1,0;0,1]");
		addLink(wedgeJointS);

		wedgeJointS = new Joint("Wedge 3 to sideframe");
		wedgeJointS->connect(wedge3->getFrame("Right face central frame"),sideFrameLeft->getFrame("sframe_r_contact_plane"));
		wedgeJointS->setForceLaw(new MBSim::BilateralConstraint);
		// wedgeJointS->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointS->setForceDirection("[1;0;0]");
		// wedgeJointS->setMomentDirection("[0,0;1,0;0,1]");
		addLink(wedgeJointS);

		wedgeJointS = new Joint("Wedge 4 to sideframe");
		wedgeJointS->connect(wedge4->getFrame("Left face central frame"),sideFrameLeft->getFrame("sframe_l_contact_plane"));
		wedgeJointS->setForceLaw(new MBSim::BilateralConstraint);
		// wedgeJointS->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointS->setForceDirection("[1;0;0]");
		// wedgeJointS->setMomentDirection("[0,0;1,0;0,1]");
		addLink(wedgeJointS);
	}

	if(contactsBolster){
		setWedgeContacts(wedge1->getContour("Left face"),bolster->getContour("Contact plane right"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge2->getContour("Right face"),bolster->getContour("Contact plane left"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge3->getContour("Left face"),bolster->getContour("Contact plane right"),
				frictionCoefficient,coefRestitution);
		setWedgeContacts(wedge4->getContour("Right face"),bolster->getContour("Contact plane left"),
				frictionCoefficient,coefRestitution);
	}
	else {
		Joint *wedgeJointB = new Joint("Wedge 1 to bolster");
		wedgeJointB->connect(wedge1->getFrame("Left face central frame"),bolster->getFrame("right_wedge_plane_frame"));
		wedgeJointB->setForceLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setForceDirection("[1;0;0]");
		wedgeJointB->setMomentDirection("[0;0;1]");
		wedgeJointB->setPlotFeature(force,true);
		addLink(wedgeJointB);

		wedgeJointB = new Joint("Wedge 2 to bolster");
		wedgeJointB->connect(wedge2->getFrame("Right face central frame"),bolster->getFrame("left_wedge_plane_frame"));
		wedgeJointB->setForceLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setForceDirection("[1;0;0]");
		wedgeJointB->setMomentDirection("[0;0;1]");
		addLink(wedgeJointB);

		wedgeJointB = new Joint("Wedge 3 to bolster");
		wedgeJointB->connect(wedge3->getFrame("Left face central frame"),bolster->getFrame("right_wedge_plane_frame"));
		wedgeJointB->setForceLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setForceDirection("[1;0;0]");
		wedgeJointB->setMomentDirection("[0;0;1]");
		addLink(wedgeJointB);

		wedgeJointB = new Joint("Wedge 4 to bolster");
		wedgeJointB->connect(wedge4->getFrame("Right face central frame"),bolster->getFrame("left_wedge_plane_frame"));
		wedgeJointB->setForceLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setMomentLaw(new MBSim::BilateralConstraint);
		wedgeJointB->setForceDirection("[1;0;0]");
		wedgeJointB->setMomentDirection("[0;0;1]");
		addLink(wedgeJointB);
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

void BarberTruck::setWedgeContacts(Contour *wedgeFace,
		Contour *otherFace,
		double frictionCoefficient,
		double coefRestitution,
		bool observerActive){

	Contact *contact;
	std::stringstream contactName;

	for (unsigned i=0;
			i < dynamic_cast<MBSim::CompoundContour*>(wedgeFace)->getNumberOfElements();
			i++){
		contactName.str("");
		contactName << "Contact_" << wedgeFace->getParent()->getName() <<
				wedgeFace->getName() <<
				"-" << otherFace->getName() << "-" << i;
		contact = new Contact ( contactName.str() );
		contact->connect ( dynamic_cast<MBSim::CompoundContour*>(wedgeFace)->getContour(i),
				otherFace );
		contact->setNormalForceLaw ( new UnilateralConstraint );
		contact->setNormalImpactLaw ( new UnilateralNewtonImpact ( coefRestitution ) );
		contact->setTangentialImpactLaw ( new SpatialCoulombImpact ( frictionCoefficient ) );
		contact->setTangentialForceLaw ( new SpatialCoulombFriction ( frictionCoefficient ) );
		contact->setPlotFeature(generalizedForce,true);
		contact->setPlotFeature(generalizedRelativePosition,true);
		addLink(contact);

		if (observerActive)
		{
			ContactObserver *observer = new ContactObserver(std::string("Contact_" + contactName.str()));
			observer->setContact(contact);
			observer->enableOpenMBVContactPoints(_size=0.05);
			observer->enableOpenMBVNormalForce(_size=0.05);
			observer->enableOpenMBVTangentialForce(_size=0.05);
			addObserver(observer);
		}
	}
}
