/*
 * sideframe.cpp
 *
 *  Created on: 26 de jan de 2018
 *      Author: leonardo
 */

#include "sideframe.h"


using namespace fmatvec;
using namespace MBSim;

Sideframe::Sideframe (const std::string &name) : RigidBody(name), sideFrameWidth(1.0), sideFrameHeight(1.0), angleSideFrame(0.0), wheelBase(2.0)
{

	std::shared_ptr<OpenMBV::IvBody> cad=OpenMBV::ObjectFactory::create<OpenMBV::IvBody>();
	cad->setIvFileName("wrl/H-Lateral.wrl");
	cad->setScaleFactor(1.0/1000);
	cad->setInitialRotation(3*M_PI/2,0,M_PI/2);
	cad->setInitialTranslation(0,0.32,0);
	this->setOpenMBVRigidBody(cad);

}

void Sideframe::buildContour(void)
{
	Vec3 pos(INIT,0.0);

	// Contact plane sideframe-left wedge definition
		pos ( 0 ) = -sideFrameWidth/2;
		pos ( 1 ) = -sideFrameHeight/2;
		fmatvec::SqrMat3 rodaroda = BasicRotAIKz(-angleSideFrame);
		this->addFrame(new FixedRelativeFrame("sframe_l_contact_plane",
				pos,
				rodaroda));
		sideframeL = new Plane ( "Side frame left",
				this->getFrame("sframe_l_contact_plane") );
		sideframeL->enableOpenMBV(_transparency=0.6);
		this->addContour ( sideframeL );

		// Contact plane sideframe-right wedge definition
		pos ( 0 ) = -pos ( 0 );
		this->addFrame(new FixedRelativeFrame("sframe_r_contact_plane",
				pos,
				BasicRotAKIz(-angleSideFrame + M_PI)));
		sideframeR = new Plane ( "Side frame right",
				this->getFrame("sframe_r_contact_plane") );
		sideframeR->enableOpenMBV(_transparency=0.6);
		this->addContour ( sideframeR );


		pos ( 0 ) = 0.0;
		this->addFrame(new FixedRelativeFrame("sframe_central_contact_plane",
				pos,
				BasicRotAKIz(0.5 * M_PI)));
		ground = new Plane ( "Ground",
				this->getFrame("sframe_central_contact_plane") );
		ground->enableOpenMBV(_transparency=0.6);
		this->addContour ( ground );
}

Contour* Sideframe::getSideframeContour(int i)
{
	switch (i){
	case 0:
		return dynamic_cast<Contour*>(this->sideframeL);
	case 1:
		return dynamic_cast<Contour*>(sideframeR);
	case 2:
		return dynamic_cast<Contour*>(ground);
	default:
		std::cout << "No contour found with specified id" << std::endl;
		return 0;
	}

}


