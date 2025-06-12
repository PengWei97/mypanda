#pragma once

#include "GBAnisotropyMisoriAndTwin.h"
#include "EBSDReaderMaterialProperty.h"

/**
 * GBAnisotropyWtDeltaRho3:
 * A class to calculate grain boundary (GB) mobility considering 
 * the difference in dislocation density (Delta Rho) and twin boundary effects.
 */
class GBAnisotropyWtDeltaRho3 : public GBAnisotropyMisoriAndTwin
{
public:
  static InputParameters validParams();

  GBAnisotropyWtDeltaRho3(const InputParameters & parameters);

protected:
  // Initialize additional material properties (e.g., Delta Rho)
  virtual void initOthersMaterialProperties() override;
  
  // Compute material properties depending on the given grain pair
  virtual void calculateOthersMaterialProperties(const unsigned int & grain_i, const unsigned int & grain_j) override;

  // Compute grain boundary mobility, optionally modified by Delta Rho
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s) override;

  // Apply the Delta Rho effect to mobility if criteria are met
  virtual void calculatedGBMobilityWtRho(const Real & delta_rho, Real & mob_ij);

  // Reference to the EBSD data reader
  const EBSDReaderMaterialProperty & _ebsd_reader;

  // Flag indicating if GB mobility calculation considers Delta Rho
  const bool _is_gb_mob_with_delta_rho;

  /// Set of grain IDs for which Delta Rho influence is considered
  const std::vector<unsigned int> _spec_grains;

  /// Amplification factor applied to GB mobility for specific grains
  const Real _amplifier_factor;

  /// Delta Rho threshold above which the mobility is modified
  const Real _threshold_for_sGrains;

  /// Material property to store the computed Delta Rho value at each quadrature point
  MaterialProperty<Real> & _delta_rho;

  /// IDs of the current grain pair being processed
  unsigned int _grain_i;
  unsigned int _grain_j;
};