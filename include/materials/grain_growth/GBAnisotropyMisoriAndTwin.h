#pragma once

#include "GBAnisotropyMisoriAng.h"

/**
 * GBAnisotropyMisoriAng:
 * Extends GBAnisotropyMisoriAng to account for twin boundary properties.
 */
class GBAnisotropyMisoriAndTwin : public GBAnisotropyMisoriAng
{
public:
  static InputParameters validParams();

  GBAnisotropyMisoriAndTwin(const InputParameters & parameters);

protected:
  // initialize _twin_boundary_type
  virtual void initOthersMaterialProperties() override;
  
  // Override methods for grain boundary energy and mobility calculations
  virtual Real calculateGBEnergy(const MisorientationAngleData & misori_s) override;
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s) override;

  // Twin boundary properties for HCP_Ti
  const Real _TT1_sigma;
  const Real _CT1_sigma;
  const Real _TT1_mob;
  const Real _CT1_mob;

  // Twin boundary properties for FCC_Ni
  const Real _Sigma9_sigma;
  const Real _Sigma3_sigma;
  const Real _Sigma9_mob;
  const Real _Sigma3_mob;

  // Material property to indicate the twin boundary type
  MaterialProperty<Real> & _twin_boundary_type;
};
