#pragma once

#include "GBAnisotropyWtDeltaRho.h"

/**
 * GBAnisotropy5SelectGrain:
 * xxx
 */
class GBAnisotropy5SelectGrain : public GBAnisotropyWtDeltaRho
{
public:
  static InputParameters validParams();

  GBAnisotropy5SelectGrain(const InputParameters & parameters);

protected:
  // Calculate grain boundary mobility with selective grain mobility adjustment
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s) override;

  const bool _is_select_grains_for_high_mobility; // Flag to indicate if high mobility grains are selected
  std::unordered_map<unsigned int, Real> _selected_grain_mobility_factors; // Map to store mobility for selected grains
  const Real _mobility_decay_rate; // High mobility degradation factor
};