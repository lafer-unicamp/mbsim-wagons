/* Copyright (C) 2019 Leonardo Baruffaldi
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
 * Contact: leonardo.baruffaldi@ifsp.edu.br
 */

#ifndef _CONTACT_KINEMATICS_WHEEL_RAIL_H_
#define _CONTACT_KINEMATICS_WHEEL_RAIL_H_

#include "mbsim/frames/contour_frame.h"
#include "mbsim/contours/point.h"
#include "mbsim/contours/circle.h"
#include "mbsim/contact_kinematics/contact_kinematics.h"
#include "wheel_profile.h"
#include "rail_profile.h"

namespace MBSim {

  class WheelProfile;
  class RailProfile;

  /**
   * \brief pairing 2D wheel profile to 2D rail profile
   * \author Leonardo Baruffaldi
   * * \date 2019-03-22 initial commit
   */
  class ContactKinematicsWheelRail : public ContactKinematics {
    public:
	  /**
	   * \brief constructor
	   */
	  ContactKinematicsWheelRail() {}

	  /**
	   * \brief destructor
	   */
	  ~ContactKinematicsWheelRail() override;


      /* INHERITED INTERFACE */
      virtual void assignContours(const std::vector<Contour*> &contour);
      virtual void updateg(double &g, std::vector<ContourFrame*> &cFrame, int index = 0);
      virtual void updatewb(fmatvec::Vec &wb, double g, std::vector<ContourFrame*> &cFrame);
      /***************************************************/

    private:
      /**
       * \brief contour index
       */
      int iwheel, irail;

      /**
       * \brief contour classes
       */
      WheelProfile *wheel;
      RailProfile *rail;

  };

}

#endif
