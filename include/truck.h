/*
 * truck.h
 *
 *  Created on: Oct 30, 2013
 *      Author: leonardo
 */

#ifndef TRUCK_H
#define TRUCK_H

#include "mbsim/group.h"
#include "wedge.h"
#include "wheelset.h"
#include "sideframe.h"

class Truck : public MBSim::Group
{
public:
  Truck(const std::string& name);

  /**
   * Get type of truck
   * @return name of the type of truck
   */
  std::string getModelType(void){ return typeOfTruck; };
private:
  std::string typeOfTruck;
protected:
  /**
   * Set the type of truck (to be use in constructors)
   * @param typeOfTruck_ string with the name of the type
   */
  void setModelType (const std::string& typeOfTruck_){ typeOfTruck = typeOfTruck_;};
};

#endif /* TRUCK_H */
