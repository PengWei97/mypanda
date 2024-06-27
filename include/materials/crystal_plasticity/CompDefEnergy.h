#pragma once

#include "Material.h"
#include "GrainTrackerMateProp.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "DerivativeMaterialPropertyNameInterface.h"
/**
 * This material is used to compute the weigthed material properties 
 * for crystal plasticity model.
 */
class CompDefEnergy : public Material,
                         public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  CompDefEnergy(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const std::string _base_name;
  std::string _elasticity_tensor_name;

  const GrainTrackerMateProp & _grain_tracker;
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;  

  std::vector<const MaterialProperty<Real> *> _elastic_energy_gr;
  std::vector<const MaterialProperty<Real> *> _plastic_energy_gr;
  MaterialProperty<Real> & _elastic_energy;
  MaterialProperty<Real> & _plastic_energy;
  std::vector<MaterialProperty<Real> *> _D_elastic_energy_gr;
  std::vector<MaterialProperty<Real> *> _D_plastic_energy_gr;

  // Conversion factor from Joules to eV
  const Real _JtoeV;
  Real _length_scale;
  Real _pressure_scale;
};