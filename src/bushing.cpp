/*
 * bushing.cpp
 *
 *  Created on: 19 de set de 2018
 *      Author: leonardo
 */

#include "bushing.h"

Bushing::Bushing(const std::string &name) : ElasticJoint(name){

	this->setGeneralizedForceFunction(func);

}


LinearElasticFunctionWithClearances::LinearElasticFunctionWithClearances() : LinearElasticFunction(){

}

