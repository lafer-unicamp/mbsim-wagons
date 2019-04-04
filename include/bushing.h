/*
 * bushing.h
 *
 *  Created on: 19 de set de 2018
 *      Author: leonardo
 */

#ifndef BUSHING_H_
#define BUSHING_H_

#include <mbsim/links/elastic_joint.h>
#include <mbsim/functions/kinetics/linear_elastic_function.h>
#include <fmatvec/fmatvec.h>

using namespace MBSim;


class LinearElasticFunctionWithClearances: public MBSim::LinearElasticFunction {
public:
	LinearElasticFunctionWithClearances();

	void setClearances(fmatvec::VecV _clear){ this->clearances = _clear;}
	virtual fmatvec::VecV operator() (const fmatvec::VecV& q, const fmatvec::VecV& u) {return K * (q - clearances) + D * u;};

protected:
	fmatvec::VecV clearances;

};



class Bushing: public MBSim::ElasticJoint {
public:
	Bushing(const std::string &name);
	void setClearances(fmatvec::VecV _clear){ this->func->setClearances(_clear); }
	void setStiffnessMatrix(const fmatvec::SymMatV &K_){ this->func->setStiffnessMatrix(K_); }
	void setDampingMatrix(const fmatvec::SymMatV &D_){ this->func->setStiffnessMatrix(D_); }

protected:
	fmatvec::VecV clearances;
	LinearElasticFunctionWithClearances *func{nullptr};
};

#endif /* BUSHING_H_ */
