//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeElasticityTensor.h"
#include "GrainTrackerMateProp.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"

/**
 * CompPolyElastTensorCPCoupled defines an elasticity tensor material
 * object for crystal plasticity and phase field coupled models. 
 * This class accepts either Bunge Euler angles or a 3x3 rotation matrix, 
 * from which the 'passive' rotation matrix is generated. 
 * This rotation matrix is used by the crystal plasticity models
 * to rotate the crystal slip system direction and plane normals into the
 * user-specified orientation.
 */
class CompPolyElastTensorCPCoupled : public ComputeElasticityTensor
{
public:
  static InputParameters validParams();

  CompPolyElastTensorCPCoupled(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  /**
   * Defines the constant rotation matrix from the user specified
   * Bunge Euler Angles or user-supplied rotation matrix.
   */
  virtual void computeQpElasticityTensor() override;

  /**
   * Calculates the contributions of each grain to the elasticity tensor and updates stateful properties.
   * @param op_to_grains Vector of active order parameters corresponding to grain IDs.
   * @return Sum of interpolation factors for grains.
   */
  Real computeGrainContributions(const std::vector<unsigned int> & op_to_grains);

  /**
   * Updates the derivatives of the elasticity tensor with respect to the order parameters.
   * @param op_to_grains Vector of active order parameters corresponding to grain IDs.
   * @param sum_h Normalization factor for the derivatives.
   */
  void updateElasticityTensorDerivatives(const std::vector<unsigned int> & op_to_grains, Real sum_h);


  const GrainTrackerMateProp & _grain_tracker;
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;

  Real _length_scale;
  Real _pressure_scale;
  const Real _JtoeV;

  /// Material property that stores the values of the Euler Angles for postprocessing
  std::vector<MaterialProperty<RealVectorValue> *> _Euler_angles_mat_prop_gr;

  /// Crystal Rotation Matrix used to rotate the slip system direction and normal
  std::vector<MaterialProperty<RankTwoTensor> *> _crysrot_gr;

  std::vector<MaterialProperty<RankFourTensor> *> _D_elasticity_tensor_gr;
};
