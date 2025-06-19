#pragma once

#include "GBAnisotropyMisoriBase.h"
#include "MooseEnum.h"
#include "EulerAngleProvider.h"
#include "GrainTracker.h"
#include "MisorientationAngleCalculator.h"

/**
 * GBAnisotropyMisoriAng:
 * This class computes the anisotropic properties of grain boundaries (GB) 
 * based on misorientation angles. It supports both energy and mobility 
 * anisotropy calculations according to the crystal structure and misorientation state.
 */
class GBAnisotropyMisoriAng : public GBAnisotropyMisoriBase
{
public:
  static InputParameters validParams();

  GBAnisotropyMisoriAng(const InputParameters & parameters);

protected:
  /// Main method to compute sigma (GB energy) and mobility for specific grain boundaries.
  virtual void computeGBProperties() override;

  /// Compute sigma and mobility for each pair of grains based on their IDs.
  void computeSigmaAndMobility(const std::vector<unsigned int> & var_index,
                                      const std::vector<unsigned int> & grain_id_index,
                                      Real &sigma_min, Real & sigma_max,
                                      Real &mob_min, Real & mob_max);

  /// Calculate GB energy based on the Read-Shockley model.
  virtual Real calculateGBEnergy(const MisorientationAngleData & misori_s);

  /// Calculate GB mobility based on the sigmoidal law.
  virtual Real calculateGBMobility(const MisorientationAngleData & misori_s);

  /// Fill symmetric properties (sigma and mobility) with averaged values for consistency.
  void fillSymmetricProperties(Real sigma_min, Real sigma_max, Real mob_min, Real mob_max);

  /// initialize other material properties.
  virtual void initOthersMaterialProperties() {};

  /// calculate other material properties based on grain id.
  virtual void calculateOthersMaterialProperties(const unsigned int & grain_i, const unsigned int & grain_j) {};

  /// Misorientation data for grain boundaries.
  MisorientationAngleData _misori_s;

  /// Crystal structure type (FCC, BCC, HCP).
  MisorientationAngleCalculator::CrystalType _crystal_structure;

  /// Execution time threshold for triggering GB property updates.
  const Real _execution_time;

  /// References to user objects for tracking grain information and orientation.
  const GrainTracker & _grain_tracker;
  const EulerAngleProvider & _euler;

  /// Flags to determine if energy and mobility anisotropy should be considered.
  const bool _gb_energy_anisotropy;
  const bool _gb_mobility_anisotropy;

  /// Constants for the sigmoidal law used in mobility calculations.
  const Real _B;
  const Real _n;

  /// Material properties for storing misorientation angle and grain type.
  MaterialProperty<Real> & _misori_angle;

  /// current grain id
  mutable unsigned int _current_grain_i = 0;
  mutable unsigned int _current_grain_j = 0;
};
