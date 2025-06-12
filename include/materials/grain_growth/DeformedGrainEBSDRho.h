#pragma once

#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "EBSDReaderMaterialProperty.h"

// Forward Declarations
class GrainTrackerInterface;

/**
 * Computes deformation energy based on GNDs from EBSD data.
 */
class DeformedGrainEBSDRho : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  DeformedGrainEBSDRho(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual Real getRhoWtTime(const unsigned int & grain_id) const;
  
  const unsigned int _op_num; // total number of grains
  const std::vector<const VariableValue *> _vals; // order parameter values

  // Simulation parameters
  const Real _length_scale;
  const Real _time_scale;
  const Real _Elas_Mod; // the elastic modulus
  const Real _Burg_vec; // the Length of Burger's Vector
  const Real _stored_factor;
  const Real _execution_time; // time to perform grain boundary anisotropy
  const Real _JtoeV; // Joule to eV conversion
  
  const GrainTrackerInterface & _grain_tracker; // Grain tracker object
  const EBSDReaderMaterialProperty & _GNDs_provider;

  // Material properties
  MaterialProperty<Real> & _rho_eff; // the average effective dislocation density
  std::vector<MaterialProperty<Real> *> _D_stored_energy;
  MaterialProperty<Real> & _feature_id;
};
