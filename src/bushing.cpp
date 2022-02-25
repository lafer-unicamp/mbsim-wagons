/*
 * bushing.cpp
 *
 *  Created on: 19 de set de 2018
 *      Author: leonardo
 */

#include "bushing.h"

/*!
   \brief Default constructor
   \param name name of the element
   \pre "Pre-conditions"
   \post "Post-conditions"
   \return The created bushing
*/

Bushing::Bushing(const std::string &name) : ElasticJoint(name){

	this->setGeneralizedForceFunction(func);

}


LinearElasticFunctionWithClearances::LinearElasticFunctionWithClearances() : LinearElasticFunction(){

}
