/*
 * wheel.h
 *
 *  Created on: Nov 2, 2013
 *      Author: leonardo
 */

#ifndef WHEEL_H
#define WHEEL_H

#include <string.h>
#include <mbsim/objects/rigid_body.h>

#include <openmbvcppinterface/rotation.h>

class Wheel : public MBSim::RigidBody
{
public:
  Wheel(const std::string& name);
  virtual
  ~Wheel();
private:
  /**
   * Name of the profile (standard) of the wheel
   */
  std::string profileName;
protected:
  void setProfileName(std::string profName_){profileName = profName_;};
  std::string getProfileName(void){return profileName;};

  /**
   * Vector of PolygonPoints with the wheel profile
   */
  std::vector<OpenMBV::PolygonPoint*>* profile;
  void setProfileContour( std::vector<OpenMBV::PolygonPoint*>* prof_ ){ profile = prof_; };
  std::vector<OpenMBV::PolygonPoint*>* getProfileContour(void){ return profile;};
};

#endif /* WHEEL_H */
