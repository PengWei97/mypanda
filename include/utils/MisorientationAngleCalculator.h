#pragma once

#include "Moose.h"
#include "EulerAngles.h"
#include "EulerAngleProvider.h"

typedef Eigen::Quaternion<Real> QuatReal; // Type alias for quaternion using Real type

// Enum for twin types
enum class TwinType {Sigma3_FCC, Sigma9_FCC, TT1_HCP, CT1_HCP, NONE};

// Enum for key orientation types
enum class QuatType {getTwinning, getSSymm, getCSymm};
// enum class CrystalType {FCC, BCC, HCP};

struct MisorientationAngleData{Real _misor = -1.0; bool _is_twin = false; TwinType _twin_type = TwinType::NONE;};

/**
 * Namespace for misorientation angle calculations.
 * This namespace contains functions to calculate misorientation angles 
 * and to perform quaternion operations.
 */
namespace MisorientationAngleCalculator
{
  // Enum for crystallographic structure types
  enum CrystalType {FCC, BCC, HCP};

  // function 1: Calculate the misorientation angle between two Euler angles.
  MisorientationAngleData calculateMisorientaion(EulerAngles & Euler1, EulerAngles & Euler2, MisorientationAngleData & s, CrystalType crystal_type = HCP);

  // function 2: Obtain key quaternions based on the type of key orientation.
  std::vector<QuatReal> getKeyQuat(const QuatType & quat_type,  CrystalType crystal_type = HCP);

  // function 3.1: Compute the scalar dot product using quaternions.
  Real dotQuaternion(const QuatReal & o1, const QuatReal & o2, 
                     const std::vector<QuatReal> & qcs, 
                     const std::vector<QuatReal> & qss);

  // function 3.2: computes inv(o1) .* o2 usig quaternion
  QuatReal itimesQuaternion(const QuatReal & q1, const QuatReal & q2);

  // function 3.3: Compute the product of the inverse of the first quaternion and the second quaternion.
  Real dotOuterQuaternion(const QuatReal & rot1, const std::vector<QuatReal> & rot2);

  // function 3.4: X*Y is the matrix product of X and Y. ~twice~
  Real mtimes2Quaternion(const QuatReal & q1, const std::vector<QuatReal> & q2, const QuatReal & qTwin);  
}