//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeFiniteStrainElasticStress.h"

#include "CPStressFullCoupledUpdateBase.h"
#include "ComputeCrystalPlasticityEigenstrainBase.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "GrainTrackerMateProp.h"

/**
 * CompMultiCPStressFullCoupled (used together with CPStressFullCoupledUpdateBase)
 * uses the multiplicative decomposition of the deformation gradient and solves the PK2 stress
 * residual equation at the intermediate configuration to evolve the material state. The internal
 * variables are updated using an iterative predictor-corrector algorithm. Backward Euler
 * integration rule is used for the rate equations.
 *
 * This material that is not called by MOOSE because of the compute=false flag
 * set in the parameter list.
 */
class CompMultiCPStressFullCoupled : public ComputeFiniteStrainElasticStress
{
public:
  static InputParameters validParams();

  CompMultiCPStressFullCoupled(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /**
   * Updates the stress (PK2) at a quadrature point by calling constiutive
   * relationship as defined in a child class
   * Solves stress residual equation using Newton - Raphson: Updates slip
   * system resistances iteratively
   */
  virtual void updateStress(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult);

  /**
   * initializes the stateful properties such as PK2 stress, resolved shear
   * stress, plastic deformation gradient, slip system resistances, etc.
   * This class is often overwritten by inherting classes.
   */
  virtual void initQpStatefulProperties() override;

  /**
   * Calls the residual and jacobian functions used in the stress update
   * algorithm.
   */
  void calculateResidualAndJacobian();

  /**
   * Reset the PK2 stress and the inverse deformation gradient to old values and
   * provide an interface for inheriting classes to reset material properties
   */
  void preSolveQp();

  /**
   * Solve the stress and internal state variables (e.g. slip increment,
   * slip system resistance) at each qp points
   */
  void solveQp();

  /**
   * Save the final stress and internal variable values after the iterative solve.
   */
  void postSolveQp(RankTwoTensor & stress_new, RankFourTensor & jacobian_mult);

  /**
   * Solves the internal variables stress as a function of the slip specified
   * by the constitutive model defined in the inheriting class
   */
  void solveStateVariables();

  /**
   * solves for stress, updates plastic deformation gradient.
   */
  void solveStress();

  /**
   * Calculate stress residual as the difference between the stored material
   * property PK2 stress and the elastic PK2 stress calculated from the
   * constitutively defined equivalent_slip_increment.
   * The equivalent_slip_increment is passed in as an input arguement.
   */
  void calculateResidual();

  /**
   * Calculates the jacobian as
   * $\mathbf{J} = \mathbf{I} - \mathbf{C} \frac{d\mathbf{E}^e}{d\mathbf{F}^e}
   * \frac{d\mathbf{F}^e}{d\mathbf{F}^P^{-1}} \frac{d\mathbf{F}^P^{-1}}{d\mathbf{PK2}}$
   */
  void calculateJacobian();

  ///@{Calculates the tangent moduli for use as a preconditioner, using the elastic or elastic-plastic option as specified by the user
  void calcTangentModuli(RankFourTensor & jacobian_mult);
  void elasticTangentModuli(RankFourTensor & jacobian_mult);
  void elastoPlasticTangentModuli(RankFourTensor & jacobian_mult);
  ///@}

  /// performs the line search update
  bool lineSearchUpdate(const Real & rnorm_prev, const RankTwoTensor & dpk2);

  /**
   * Calculates the deformation gradient due to eigenstrain
   */
  void calculateEigenstrainDeformationGrad();

  /// number of plastic models
  const unsigned _num_models;

  /// The user supplied cyrstal plasticity consititutive models
  std::vector<CPStressFullCoupledUpdateBase *> _models;

  /// number of eigenstrains
  const unsigned _num_eigenstrains;

  /// The user supplied cyrstal plasticity eigenstrains
  std::vector<ComputeCrystalPlasticityEigenstrainBase *> _eigenstrains;

  /// optional parameter to define several mechanical systems on the same block, e.g. multiple phases
  const std::string _base_name;

  /// add some paramerts for the coupled phase field model
  const GrainTrackerMateProp & _grain_tracker;
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;
  unsigned int _active_op_index;
  unsigned int _grain_id;

  /// Elasticity tensor as defined by a separate class
  std::vector<const MaterialProperty<RankFourTensor> *> _elasticity_tensor_gr;

  /**
   * Crystal rotation in the original, or reference, configuration as defined by
   * Euler angle arguments in the ComputeElasticityTensor classes
   */
  std::vector<const MaterialProperty<RankTwoTensor> *> _crysrot_gr;

  /// Stress residual equation relative tolerance
  Real _rtol;
  /// Stress residual equation absolute tolerance
  Real _abs_tol;

  /// Residual tensor
  RankTwoTensor _residual_tensor;
  /// Jacobian tensor
  RankFourTensor _jacobian;

  /// Maximum number of iterations for stress update
  unsigned int _maxiter;
  /// Maximum number of iterations for internal variable update
  unsigned int _maxiterg;

  /// Type of tangent moduli calculation
  const enum class TangentModuliType { EXACT, NONE } _tan_mod_type;

  /// Maximum number of substep iterations
  unsigned int _max_substep_iter;

  /// time step size during substepping
  Real _substep_dt;

  /// Flag to activate line serach
  bool _use_line_search;

  /// Minimum line search step size
  Real _min_line_search_step_size;

  /// Line search bisection method tolerance
  Real _line_search_tolerance;

  /// Line search bisection method maximum iteration number
  unsigned int _line_search_max_iterations;

  /// strain formulation
  const enum class LineSearchMethod { CutHalf, Bisection } _line_search_method;

  ///@{Plastic deformation gradient RankTwoTensor for the crystal
  std::vector<MaterialProperty<RankTwoTensor> *> _plastic_deformation_gradient_gr;
  std::vector<const MaterialProperty<RankTwoTensor> *> _plastic_deformation_gradient_gr_old;
  ///@}

  ///@{ Generalized eigenstrain deformation gradient RankTwoTensor for the crystal
  MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient_old;
  ///@}

  ///@{Total deformation gradient RankTwoTensor for the crystal
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  ///@}

  ///@{Second Piola-Kirchoff stress measure
  std::vector<MaterialProperty<RankTwoTensor> *> _pk2_gr;
  std::vector<const MaterialProperty<RankTwoTensor> *> _pk2_gr_old;
  ///@}

  /// Lagrangian total strain measure for the entire crystal
  MaterialProperty<RankTwoTensor> & _total_lagrangian_strain;
  MaterialProperty<RankTwoTensor> & _elastic_lagrangian_strain;

  /// Cauchy stress for grain id
  std::vector<MaterialProperty<RankTwoTensor> *> _stress_gr;
  std::vector<MaterialProperty<RankFourTensor> *> _Jacobian_mult_gr;

  /**
   * Tracks the rotation of the crystal during deformation
   * Note: this rotation tensor is not applied to the crystal lattice
   */
  std::vector<MaterialProperty<RankTwoTensor> *> _updated_rotation_gr;

  ///@{Helper deformation gradient tensor variables used in iterative solve
  RankTwoTensor _temporary_deformation_gradient;
  RankTwoTensor _elastic_deformation_gradient;
  RankTwoTensor _inverse_plastic_deformation_grad;
  RankTwoTensor _inverse_plastic_deformation_grad_old;
  RankTwoTensor _inverse_eigenstrain_deformation_grad;
  ///@}

  /// Flag to print to console warning messages on stress, constitutive model convergence
  const bool _print_convergence_message;

  /// Flag to check whether convergence is achieved or if substepping is needed
  bool _convergence_failed;

  /// Flag to chech whether the order parameter is actived
  bool _is_active_op;

  ///@{ Used for substepping; Uniformly divides the increment in deformation gradient
  RankTwoTensor _delta_deformation_gradient;
  RankTwoTensor _temporary_deformation_gradient_old;
  ///@}

  /// Scales the substepping increment to obtain deformation gradient at a substep iteration
  Real _dfgrd_scale_factor;

  /// define elastic energy and plastic energy
  std::vector<MaterialProperty<Real> *> _elastic_energy_gr;
  std::vector<MaterialProperty<Real> *> _plastic_energy_gr;
};
