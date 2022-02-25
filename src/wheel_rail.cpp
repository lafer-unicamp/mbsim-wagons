/*
 * wheel_rail.cpp
 * Copyright (C) 2019 MBSim Development Team
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 *  Created on: 29 de mar de 2019
 *      Author: Leonardo Baruffaldi
 *	   Contact: leonardo.baruffaldi@ifsp.edu.br
 */

#include  "wheel_rail.h"

namespace MBSim{

ContactKinematicsWheelRail::~ContactKinematicsWheelRail() {



}

void ContactKinematicsWheelRail::assignContours(
		const std::vector<Contour*>& contour) {
	if(dynamic_cast<RailProfile*>(contour[0])){
		irail = 0;
		iwheel = 1;
		rail = static_cast<RailProfile*>(contour[0]);
		wheel = static_cast<WheelProfile*>(contour[1]);
	} else {
		irail = 1;
		iwheel = 0;
		rail = static_cast<RailProfile*>(contour[1]);
		wheel = static_cast<WheelProfile*>(contour[0]);
	}
}

void ContactKinematicsWheelRail::updateg(double& g,
		std::vector<ContourFrame*>& cFrame, int index) {
	/* TODO implement GJK distance algorithm
	 * 1. Calcular ponto de contato
	 * 1.1. Determinar posição global do perfil da roda
	 * 1.2. Determinar posição global do perfil do trilho
	 * 1.3. Rodar algoritmo de localização do ponto de contato
	 * 2. Calcular normal
	 * 3. Atualizar os cFrames dos contornos para serem compatíveis com os pontos de contato
	 */

	// Parte 1: determinação do ponto de contato roda-trilho
	// Parte 1.1.: posição global do perfil da roda
	// wheelProfilePoints é uma relação de pontos no plano
	// de contato roda-trilho.
	Vec3 rWheel = wheel->getFrame()->evalPosition();
	MatVx2 wheelProfilePoints = wheel->getPoints();
	// Parte 1.2.: posição global do perfil do trilho
	/*
	 * estabelecer a posição de um frame em relação à
	 * linha média do trilho no local mais próximo ao
	 * centro de massa do rodeiro
	 *
	 */
	Vec3 rRail = rail->getFrame()->evalPosition();
	MatVx2 railProfilePoints = rail->getPoints();

	cFrame[iwheel]->setPosition(rWheel);
	cFrame[irail]->setPosition(rRail);
	cFrame[iwheel]->enableOpenMBV();
	cFrame[irail]->enableOpenMBV();

	findContactPoints(wheelProfilePoints,railProfilePoints);

	g = 0.003;

}

void ContactKinematicsWheelRail::updatewb(Vec& wb, double g,
		std::vector<ContourFrame*>& cFrame) {
}

/**
 * Implementation of the improved search method from Gino van den Bergen's
 * paper on the fast implementation of the GJK distance algorithm
 * G. V. den Bergen, “A Fast and Robust GJK Implementation for Collision Detection
 * of Convex Objects”, Journal of Graphics Tools, vol. 4, nº 2, p. 7–25, jan. 1999.
 */
void ContactKinematicsWheelRail::findContactPoints(MatVx2 wheelProfilePoints,
		MatVx2 railProfilePoints) {
	/* CHECK THE SMALLEST DISTANCE BETWEEN WHEEL AND RAIL
	 * as a preliminary adjustment, we separate the convex sets on the wheel and
	 * rail to accelerate the contact calculation
	 */


	/**
	 * \brief a matrix containig the sets of initial and final
	 * points of the convex segments of the profile
	 */

	// TODO: move this code as a function under ProfileContour
	MatVx2 wheelConvex;
	MatVx2 railConvex;

	wheelConvex = wheel->getConvexSetIndexes();
	railConvex = rail->getConvexSetIndexes();

	/* we start with a coarse check by numbering the segments of the wheel and rail pro-
	 * files and checking their proximity
	 */

	//		double mindist = 1e6;
	//		int noderail, nodewheel;
	//		bool madeContact = false;
	//
	//		noderail = 0;
	//		nodewheel = 0;

	double g;

	for (int i=0;i<railConvex.rows();i++){
		for (int j=0;j<wheelConvex.rows();j++){
			g = gjkDistance(wheelProfilePoints,
					railProfilePoints,1e-3,false);
			g = g*1;
		}
	}


}

/**
 * Calculation of the Minkowski sum of two polyhedra A and B
 * @param idxmap is a pass-by-reference matrix of integers with the mapping of the indexes os A and B
 * on the resulting sum <A-B>
 */
MatVx2 ContactKinematicsWheelRail::minkowskiSum(MatVx2 A,
		MatVx2 B, MatVx2I& idxmap){


	unsigned k = A.rows();
	unsigned p = B.rows();
	int c = 0;
	MatVx2 mSum(k*p,0);
	idxmap.resize(k*p);

	for (unsigned i=0;i<p;i++){
		for (unsigned j=0;j<k;j++){
			mSum(c,0) = -A(j,0) + B(i,0);
			mSum(c,1) = -A(j,1) + B(i,1);
			idxmap(c,0) = j;
			idxmap(c,1) = i;
			c++;
		}
	}

	return mSum;

}

/**
 * based on Gino van den Bergen's paper on GJK algorithm
 * G. V. den Bergen, “A Fast and Robust GJK Implementation for Collision Detection
 * of Convex Objects”, Journal of Graphics Tools, vol. 4, nº 2, p. 7–25, jan. 1999.
 * @param A : coordinates of the vertices of some polyhedron
 * @param B : coordinates of the vertices of some other polyhedron
 * @param maxPenetration : maximum allowed penetration between polyhedra
 * @param verbose : should verbose mode be active?
 * @return d : the distance between to polyhedra
 */

double ContactKinematicsWheelRail::gjkDistance(MatVx2 A,
		MatVx2 B, double maxPenetration, bool verbose) {

	bool contact = true;
	bool closeEnough = false;
	//	unsigned int dim = A.cols(); // the dimensionality of the problem
	double epsilon = 1e-8;
	VecV u, v, w, vtemp;
	MatVx2 Y, W, Ymod(3), Wtemp, Yred(2);
	VecV idY, idW, lambda(3,INIT,0), lout(3), deltaX;
	RowVecV deltaY;
	MatV combs, keep, D;
	MatVx2I idxmab;
	unsigned int wIdx = 0;
	unsigned int n = 0;
	double mu = 0;
	unsigned int iter = 0; // iterations
	double normv; // the euclidian norm of v
	double delta; // normalized projection of w in v


	// the Minskowski sum of the two polyhedra
	MatVx2 minkowski = minkowskiSum(A,B,idxmab);
	unsigned int nmab = minkowski.rows();

	// select the first test points
	unsigned int ptAid = 0;
	unsigned int ptBid = 0;
	MatVx2 ptA = A(Range<Var, Var>(ptAid,0),Range<Var,Var>(ptAid,1));
	MatVx2 ptB = B(Range<Var, Var>(ptBid,0),Range<Var,Var>(ptBid,1));

	// algorithm start
	v = (VecV)minkowski.T().col(0);

	while (!closeEnough & contact & (iter < nmab)){

		polytopeMap(-v,minkowski,w,wIdx);

		normv = nrm2(v);

		if (normv != 0){
			delta =  v.T()*w / normv;
		}
		else {
			delta = 0;
		}

		// delta >= 0 means the polyhedra are not in contact
		if (delta >= 0){
			contact = false;
		}

		mu = std::max(mu,delta);

		closeEnough = ((normv-mu<=epsilon * normv) | (normv < maxPenetration));

		if (closeEnough) std::cout << "Contato" << std:: endl;

		if (!closeEnough){

			Y = W;
			idY = idW;
			n = Y.rows();
			Y.resize(n+1,Y.cols());
			idY.resize(n+1);
			Y.set(n,w.T());
			idY(n) = wIdx;



			// assemble combinations lists
			// first column: number of elements in each combination
			// i.e.: 2 1 2 0 means combination of 1 and 2 between 1, 2, and 3.
			// 2 1 3 0 means combination of 1 and 3 between 1, 2 ,and 3

			switch (n+1)
			{
			case 1:
				combs = MatV("[1,0]");
				break;
			case 2:
				combs = MatV("[2 0 1]");
				break;
			case 3:
				combs = MatV("[2,1,3,0;2,2,3,0]");
				break;
			}

			Ymod.init(0);
			Ymod.set(Range<Var,Var>(0,0),Range<Var,Var>(n,1),Y);

			vtemp = Ymod.row(0).T();
			Wtemp.resize(1,1);
			lout = lambda;

			for (int k = 0; k < combs.rows(); k++){

				lambda.init(0);

				// which vectors to test inside the Y set
				keep.resize(1,1);
				if (combs.cols() == 2) {
					keep.resize(1,1);
					keep(0,0) = combs(0,1);

				}
				else {
					keep = combs(Range<Var,Var>(k,0),Range<Var,Var>(k,combs(k,1)));
				}


				// we need to check if the two compared vectors are not the same
				if ((keep.cols() == 2) && (Ymod.row(keep(0,0)) == Ymod.row(keep(1,0)))) {
					keep.resize(1,1);
				}



				// if keep contains a single point (a 1-simplex), then the value of
				// u is that single point.
				if (keep.cols() == 1){
					u = Ymod.row(keep(0,0)).T();
					lambda(keep(0,0)) = 1;
					Wtemp.resize(1,2);
					Wtemp.set(0,u.T());
					idW.resize(1);
					idW(0) = idY(keep(0,0));
				}
				else{
					// otherwise, we have to calculate the vector that is
					// perpendicular to the line joining the two vertices
					// and that also belongs to the 2-simplex that is formed
					// by those two vertices
					Yred.init(0);
					D.resize(2,2);
					D.init(1);
					for (int k = 0; k < keep.rows();k++){
						Yred.set(k,Ymod.row(keep(k,0)));
					}


					deltaY = (Ymod.row(keep(1,0))-Y.row(keep(0,0)));
					D.set(0,deltaY*Yred.T());
					D.set(1,RowVec(2,INIT,1));
					//TODO finish implementation
				}

				if ( (nrm2(u) < nrm2(vtemp)) && (min(lambda) >= 0)) {
					vtemp = u;
					for (int l=0; l < keep.cols(); l++){
						Wtemp.set(l,Ymod.row(keep(0,l)));
						lout(l) = lambda(keep(0,l));
					}
					//					ptAid.resize
					//					ptAid = idxmab(idW,0);
					//					ptBid = idxmab(idW,1);
					//					ptA = lout * A.row(ptAid).T();
					//					ptB = lout * B.row(ptBid).T();
				}

			}

			v = vtemp;
			W = Wtemp;

		}// end of !closeEnough

		iter++;
	}

	if (!contact && verbose) {
		std::cout << "Routine gjkdistance terminated due to no contact found." << std::endl;
	} else if (iter >= nmab) {
		std::cout << "Routine gjkdistance terminated due to max iterations" << std::endl;
	} else if (normv < maxPenetration) {
		std::cout << "Routine gjkdistance terminated because penetration reached the maximum value allowed."
				<< std::endl;
	} else if (normv-mu<=epsilon * normv) {
		std::cout << "Routine gjkdistance ended successfully by reaching the tolerance gap" << std::endl;
	}

	if (contact) {
		return -normv;
	} else {
		return normv;
	}

}


/** polytope map (van den Bergen)
 * G. V. den Bergen, “A Fast and Robust GJK Implementation for Collision Detection
 * of Convex Objects”, Journal of Graphics Tools, vol. 4, nº 2, p. 7–25, jan. 1999.
 * @param vector - a vertex
 * @param polyhedron - an array of vertices from a polytope
 * @param sa - proximal point of the map
 * @param index - index of sa on the polyhedron list of vertices
 */
void ContactKinematicsWheelRail::polytopeMap(VecV vector, MatVx2 polyhedron, VecV& sa, unsigned int& index) {


	VecV X = (VecV)(polyhedron * vector);

	index = fmatvec::maxIndex(X);
	sa = polyhedron.T().col(index);

}

} /* END namespace MBSim */
