#pragma once

#include "DeformedGrainEBSDRho.h"

/**
 * Computes deformation energy based on GNDs from EBSD data.
 */
class DeformedGrainEBSDRho_ExtendR2 : public DeformedGrainEBSDRho
{
public:
  static InputParameters validParams();

  DeformedGrainEBSDRho_ExtendR2(const InputParameters & parameters);

protected:
  virtual Real getRhoWtTime(const unsigned int & grain_id) const override;
  Real computeRhoWithRecovery(const Real & rho_init) const;

  std::vector<unsigned int> _feature_ids; // Grain IDs for selected grains
  std::vector<Real> _set_rhos; // Set rhos for the selected grains
  const Real _threshold_DeltaRho; // Threshold for delta_rho to apply the set_rhos

  /// Optional concurrent recovery
  const bool _enable_concurrent_recovery;
  const Real _rho_default;
  const Real _rho_end;
  const Real _a_rho;

  /// Grain selector by feature ID (fast lookup)
  std::unordered_map<unsigned int, Real> _rho_map;

  const MaterialProperty<Real> & _delta_rho; // Delta rho for grain boundary anisotropy
};
