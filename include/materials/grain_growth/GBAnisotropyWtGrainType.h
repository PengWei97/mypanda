#pragma once

#include "GBAnisotropyMisoriBase.h"
#include "GrainTracker.h"
#include "EBSDReader.h"

/**
 * @class GBAnisotropyWtGrainType
 * @brief Extends GBAnisotropyMisoriBase to account for grain-type-specific properties
 *        during the computation of grain boundary (GB) energy and mobility.
 *
 * This class introduces the concept of grain type during the computation
 * of GB properties, leveraging information from the EBSDReader to adjust
 * the GB characteristics accordingly.
 */

class GBAnisotropyWtGrainType : public GBAnisotropyMisoriBase
{
public:
  static InputParameters validParams();

  GBAnisotropyWtGrainType(const InputParameters & parameters);

protected:
  virtual void computeGBProperties() override;

  void computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                               const std::vector<unsigned int> & grain_ids,
                               Real & sigma_min, Real & sigma_max,
                               Real & mob_min, Real & mob_max);

  Real calculateGBEnergy(const Real & gi, const Real & gj);
  
  Real calculateGBMobility(const Real & gi, const Real & gj);

  void fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max);

  /// References to user objects for tracking grain information.
  const EBSDReader & _ebsd_reader;
  const GrainTracker & _grain_tracker;
  const Real _select_grain_type;
  const bool _gb_energy_anisotropy;
  const bool _gb_mobility_anisotropy;
  const Real _execution_time;

  /// Material property that stores the type of grain at each quadrature point.
  MaterialProperty<Real> & _grain_type;
};