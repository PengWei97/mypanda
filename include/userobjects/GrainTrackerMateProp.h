#pragma once

#include "GrainTrackerElasticity.h"
#include "EulerAngleProvider.h"

/**
 * Manages material properties based on grain_id, including elasticity, 
 * Euler angles & grain rotation matrices
 */
class GrainTrackerMateProp : public GrainTrackerElasticity
{
public:
  static InputParameters validParams();

  GrainTrackerMateProp(const InputParameters & parameters);

  // Returns the Euler angles for the given grain_id
  const EulerAngles & updateEulerAngles(const unsigned int & grain_id) const;
  
protected:
  EulerAngles _angles_tmp;
};