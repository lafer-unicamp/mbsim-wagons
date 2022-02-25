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

using namespace fmatvec;

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

      void findContactPoints(MatVx2 wheelProfilePoints, MatVx2 railProfilePoints);

    private:
      /**
       * \brief contour index
       */
      int iwheel, irail;

      /**
       * \brief contour class
       */
      WheelProfile *wheel;
      RailProfile *rail;

      /**
       * \brief Implementation of the Minkowski sum
       */
      fmatvec::MatVx2 minkowskiSum(fmatvec::MatVx2 A, fmatvec::MatVx2 B, MatVx2I& idxmap);

      /**
       * \brief Implementation of GJK distance algorithm according to
       * G. V. den Bergen, “A Fast and Robust GJK Implementation for Collision Detection
       * of Convex Objects”, Journal of Graphics Tools, vol. 4, nº 2, p. 7–25, jan. 1999.
       */
      double gjkDistance(fmatvec::MatVx2 A,
    		  fmatvec::MatVx2 B,
			  double maxPenetration,
			  bool verbose);

      void polytopeMap(VecV vector, MatVx2 polyhedron, VecV& sa, unsigned int& index);

  };

}

#endif
