//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GBAnisotropyMisoriBase.h"
#include "MooseEnum.h"
#include "EulerAngleProvider.h"
#include "GrainTracker.h"
#include "MisorientationAngleCalculator.h"

// Forward Declarations

/**
 * GBAnisotropyMisori is created based on GrainTracker for real-time acquisition of
 * sigma_ij and mob_ij based on misorientation angle.
 * - Function 1: GB anisotropy based on classic theories
 * - Function 2: Low-energy and low-mobility properties of twin interfaces
 */
class GBAnisotropyMisori : public GBAnisotropyMisoriBase
{
public:
  static InputParameters validParams();

  GBAnisotropyMisori(const InputParameters & parameters);

protected:
  // Compute sigma_ij and mob_ij for specific grain boundaries
  virtual void computeGBProperties() override;

  // Calculate GB energy based on the Read-Shockley model
  virtual Real calculatedGBEnergy(const MisorientationAngleData & misori_s);

  // Calculate GB mobility based on the sigmoidal law
  virtual Real calculatedGBMobility(const MisorientationAngleData & misori_s);

  // Determine the twinning type based on misorientation data
  Real determineTwinningType(const MisorientationAngleData & misori_s);

  // Fill symmetric properties with averaged values
  void fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max);

  // Store misorientation angle data, including twinning type
  MisorientationAngleData _misori_s;
  MisorientationAngleCalculator::CrystalType _crystal_structure;
  
  // Twin boundary properties for HCP_Ti
  const Real _TT1_sigma;
  const Real _CT1_sigma;
  const Real _TT1_mob;
  const Real _CT1_mob;

  // Twin boundary properties for FCC_Ni
  const Real _Sigma3_sigma;
  const Real _Sigma9_sigma;
  const Real _Sigma3_mob;
  const Real _Sigma9_mob;

  // GrainTracker user object for acquiring grain IDs
  const GrainTracker & _grain_tracker;
  
  // Euler angle provider user object for acquiring Euler angles
  const EulerAngleProvider & _euler; 

  // Flags to consider GB energy and mobility anisotropy
  const bool _gb_energy_anisotropy;
  const bool _gb_mobility_anisotropy;    

  // Material properties to store misorientation angle and twinning type
  MaterialProperty<Real> & _misori_angle;
  MaterialProperty<Real> & _twinning_type;
};

