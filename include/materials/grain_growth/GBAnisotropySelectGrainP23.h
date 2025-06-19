#pragma once

#include "GBAnisotropyMisoriAndTwin.h"

/**
 * GBAnisotropySelectGrainP23:
 * xxxx
 */
class GBAnisotropySelectGrainP23 : public GBAnisotropyMisoriAndTwin
{
public:
  static InputParameters validParams();

  GBAnisotropySelectGrainP23(const InputParameters & parameters);

protected:
  // Override methods for grain boundary energy and mobility calculations
  virtual Real calculateGBEnergy(const MisorientationAngleData & misori_s) override;
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s) override;

  bool isInList(const std::vector<Real> & list, const Real val) const;

  const bool _is_select_grain_type; // Flag to select grain type 1
  const std::vector<Real> _select_grain_type1; // Selected grain type 1
  const std::vector<Real> _select_grain_type2; // Selected grain type 2
  const std::vector<Real> _scaling_factors;
};
