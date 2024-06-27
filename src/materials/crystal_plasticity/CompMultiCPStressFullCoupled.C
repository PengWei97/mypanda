//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CompMultiCPStressFullCoupled.h"

#include "CPStressFullCoupledUpdateBase.h"
#include "libmesh/utility.h"
#include "Conversion.h"
#include "MooseException.h"

registerMooseObject("mypandaApp", CompMultiCPStressFullCoupled);

InputParameters
CompMultiCPStressFullCoupled::validParams()
{
  InputParameters params = ComputeFiniteStrainElasticStress::validParams();

  params.addClassDescription(
      "Crystal Plasticity base class: handles the Newton iteration over the stress residual and "
      "calculates the Jacobian based on constitutive laws from multiple material classes "
      "that are inherited from CPStressFullCoupledUpdateBase");
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");

  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors and euler angles");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");

  params.addRequiredParam<std::vector<MaterialName>>(
      "crystal_plasticity_models",
      "The material objects to use to calculate crystal plasticity stress and strains.");
  params.addParam<std::vector<MaterialName>>(
      "eigenstrain_names", {}, "The material objects to calculate eigenstrains.");
  params.addParam<MooseEnum>("tan_mod_type",
                             MooseEnum("exact none", "none"),
                             "Type of tangent moduli for preconditioner: default elastic");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residual relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residual absolute tolerance");
  params.addParam<unsigned int>("maxiter", 100, "Maximum number of iterations for stress update");
  params.addParam<unsigned int>(
      "maxiter_state_variable", 100, "Maximum number of iterations for state variable update");
  params.addParam<unsigned int>(
      "maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
  params.addParam<Real>("line_search_tol", 0.5, "Line search bisection method tolerance");
  params.addParam<unsigned int>(
      "line_search_maxiter", 20, "Line search bisection method maximum number of iteration");
  params.addParam<MooseEnum>("line_search_method",
                             MooseEnum("CUT_HALF BISECTION", "CUT_HALF"),
                             "The method used in line search");
  params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on the stress measure and general constitutive model quantinties.");
  return params;
}

CompMultiCPStressFullCoupled::CompMultiCPStressFullCoupled(
    const InputParameters & parameters)
  : ComputeFiniteStrainElasticStress(parameters),
    _num_models(getParam<std::vector<MaterialName>>("crystal_plasticity_models").size()),
    _num_eigenstrains(getParam<std::vector<MaterialName>>("eigenstrain_names").size()),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),

    _grain_tracker(getUserObject<GrainTrackerMateProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _active_op_index(0), // initialize to 0
    _grain_id(0), // initialize to 0

    _elasticity_tensor_gr(_op_num),
    _crysrot_gr(_op_num), // defined in the elasticity tensor classes for crystal plasticity

    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxiter_state_variable")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type").getEnum<TangentModuliType>()),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_line_search_step_size(getParam<Real>("min_line_search_step_size")),
    _line_search_tolerance(getParam<Real>("line_search_tol")),
    _line_search_max_iterations(getParam<unsigned int>("line_search_maxiter")),
    _line_search_method(getParam<MooseEnum>("line_search_method").getEnum<LineSearchMethod>()),
    _plastic_deformation_gradient_gr(_op_num),
    _plastic_deformation_gradient_gr_old(_op_num),
    _eigenstrain_deformation_gradient(
        _num_eigenstrains ? &declareProperty<RankTwoTensor>("eigenstrain_deformation_gradient")
                          : nullptr),
    _eigenstrain_deformation_gradient_old(
        _num_eigenstrains
            ? &getMaterialPropertyOld<RankTwoTensor>("eigenstrain_deformation_gradient")
            : nullptr),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _pk2_gr(_op_num),
    _pk2_gr_old(_op_num),
    _total_lagrangian_strain(
        declareProperty<RankTwoTensor>("total_lagrangian_strain")), // Lagrangian strain
    _elastic_lagrangian_strain(
        declareProperty<RankTwoTensor>("elastic_lagrangian_strain")), // Elastic lagrangian strain
    _stress_gr(_op_num),
    _Jacobian_mult_gr(_op_num),
    _updated_rotation_gr(_op_num),
    _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages")),
    _elastic_energy_gr(_op_num),
    _plastic_energy_gr(_op_num)
{
  _convergence_failed = false;

  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    _elasticity_tensor_gr[op_index] = &getMaterialProperty<RankFourTensor>(
        _base_name + "elasticity_tensor_gr" + std::to_string(op_index));
    _crysrot_gr[op_index] = &getMaterialProperty<RankTwoTensor>(_base_name + "crysrot_gr" + std::to_string(op_index));

    _plastic_deformation_gradient_gr[op_index] = &declareProperty<RankTwoTensor>(
        _base_name + "plastic_deformation_gradient_gr" + std::to_string(op_index));
    _plastic_deformation_gradient_gr_old[op_index] = &getMaterialPropertyOld<RankTwoTensor>(
        _base_name + "plastic_deformation_gradient_gr" + std::to_string(op_index));
    _pk2_gr[op_index] = &declareProperty<RankTwoTensor>(_base_name + "pk2_gr" + std::to_string(op_index));
    _pk2_gr_old[op_index] = &getMaterialPropertyOld<RankTwoTensor>(_base_name + "pk2_gr" + std::to_string(op_index));

    _stress_gr[op_index] = &declareProperty<RankTwoTensor>(_base_name + "stress_" + coupledName("v", op_index));
    _Jacobian_mult_gr[op_index] = &declareProperty<RankFourTensor>(_base_name + "Jacobian_mult_" + coupledName("v", op_index));

    _updated_rotation_gr[op_index] = &declareProperty<RankTwoTensor>(_base_name + "updated_rotation_gr" + std::to_string(op_index));

    // Elastic energy & plastic energy
    _elastic_energy_gr[op_index] = &declareProperty<Real>(_base_name + "elastic_energy_" + coupledName("v", op_index));
    _plastic_energy_gr[op_index] = &declareProperty<Real>(_base_name + "plastic_energy_" + coupledName("v", op_index));
  }
}

void
CompMultiCPStressFullCoupled::initQpStatefulProperties()
{
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    (*_plastic_deformation_gradient_gr[op_index])[_qp].zero();
    (*_plastic_deformation_gradient_gr[op_index])[_qp].addIa(1.0);

    (*_pk2_gr[op_index])[_qp].zero();

    (*_updated_rotation_gr[op_index])[_qp].zero();
    (*_updated_rotation_gr[op_index])[_qp].addIa(1.0);

    // Elastic energy & plastic energy
    (*_elastic_energy_gr[op_index])[_qp] = 0.0;
    (*_plastic_energy_gr[op_index])[_qp] = 0.0;
  }

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setQp(_qp);
    _models[i]->initQpStatefulProperties();
  }

  _total_lagrangian_strain[_qp].zero();
  _elastic_lagrangian_strain[_qp].zero();

  if (_num_eigenstrains)
  {
    (*_eigenstrain_deformation_gradient)[_qp].zero();
    (*_eigenstrain_deformation_gradient)[_qp].addIa(1.0);
  }
  else
  {
    // set to identity if no eigenstrain is added
    _inverse_eigenstrain_deformation_grad.zero();
    _inverse_eigenstrain_deformation_grad.addIa(1.0);
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setQp(_qp);
    _eigenstrains[i]->initQpStatefulProperties();
  }
}

void
CompMultiCPStressFullCoupled::initialSetup()
{
  // get crystal plasticity models
  std::vector<MaterialName> model_names =
      getParam<std::vector<MaterialName>>("crystal_plasticity_models");

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    CPStressFullCoupledUpdateBase * model =
        dynamic_cast<CPStressFullCoupledUpdateBase *>(&getMaterialByName(model_names[i]));

    if (model)
    {
      _models.push_back(model);
      // TODO: check to make sure that the material model is compatible with this class
    }
    else
      mooseError("Model " + model_names[i] +
                 " is not compatible with CompMultiCPStressFullCoupled");
  }

  // get crystal plasticity eigenstrains
  std::vector<MaterialName> eigenstrain_names =
      getParam<std::vector<MaterialName>>("eigenstrain_names");

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    ComputeCrystalPlasticityEigenstrainBase * eigenstrain =
        dynamic_cast<ComputeCrystalPlasticityEigenstrainBase *>(
            &getMaterialByName(eigenstrain_names[i]));

    if (eigenstrain)
      _eigenstrains.push_back(eigenstrain);
    else
      mooseError("Eigenstrain" + eigenstrain_names[i] +
                 " is not compatible with CompMultiCPStressFullCoupled");
  }
}

void
CompMultiCPStressFullCoupled::computeQpStress()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setQp(_qp);
    _models[i]->setMaterialVectorSize();
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
    _eigenstrains[i]->setQp(_qp);

  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real sum_h = 0.0;
  _stress[_qp].zero();
  _Jacobian_mult[_qp].zero();
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    _grain_id = op_to_grains[op_index];

    _active_op_index = op_index;
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->setActiveOpIndex(_active_op_index);

    // if (_grain_id == FeatureFloodCount::invalid_id)
    // {
    //   for (unsigned int i = 0; i < _num_models; ++i)
    //     _models[i]->resetQpStatefulProperties(op_index);

    //   // continue;
    // }

    _elastic_lagrangian_strain[_qp].zero();
    updateStress((*_stress_gr[_active_op_index])[_qp], (*_Jacobian_mult_gr[_active_op_index])[_qp]); // This is NOT the exact jacobian

    if (_grain_id == FeatureFloodCount::invalid_id)
    {
      for (unsigned int i = 0; i < _num_models; ++i)
        _models[i]->resetQpStatefulProperties(op_index);

      // continue;
      (*_elastic_energy_gr[op_index])[_qp] = 0.0;
      (*_plastic_energy_gr[op_index])[_qp] = 0.0;
    }
    else
    {
      (*_elastic_energy_gr[op_index])[_qp] = 0.5 * ((*_pk2_gr[_active_op_index])[_qp]).doubleContraction(_elastic_lagrangian_strain[_qp]);
      for (unsigned int i = 0; i < _num_models; ++i)
        _models[i]->calculatePlasticityEnergy((*_plastic_energy_gr[op_index])[_qp], op_index);
    }

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all cauchy stress and Jocabian matrix
    _stress[_qp] += (*_stress_gr[_active_op_index])[_qp] * h;
    _Jacobian_mult[_qp] += (*_Jacobian_mult_gr[_active_op_index])[_qp] * h;
    sum_h += h;
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _stress[_qp] /= sum_h;
  _Jacobian_mult[_qp] /= sum_h;
}

void
CompMultiCPStressFullCoupled::updateStress(RankTwoTensor & cauchy_stress,
                                                     RankFourTensor & jacobian_mult)
{
  // Does not support face/boundary material property calculation
  if (isBoundaryMaterial())
    return;

  // Initialize substepping variables
  unsigned int substep_iter = 1;
  unsigned int num_substep = 1;

  _temporary_deformation_gradient_old = _deformation_gradient_old[_qp];
  if (_temporary_deformation_gradient_old.det() == 0)
    _temporary_deformation_gradient_old.addIa(1.0);

  _delta_deformation_gradient = _deformation_gradient[_qp] - _temporary_deformation_gradient_old;

  // Loop through all models and calculate the schmid tensor for the current state of the crystal
  // lattice
  // Not sure if we should pass in the updated or the original rotation here
  // If not, then we should not need to compute the flow direction every iteration here
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->calculateFlowDirection((*_crysrot_gr[_active_op_index])[_qp]);

  do
  {
    _convergence_failed = false;
    preSolveQp();

    _substep_dt = _dt / num_substep;
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->setSubstepDt(_substep_dt);

    // calculate F^{eigen} only when we have eigenstrain
    _inverse_eigenstrain_deformation_grad.zero();
    _inverse_eigenstrain_deformation_grad.addIa(1.0);
    if (_num_eigenstrains)
      calculateEigenstrainDeformationGrad();

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _temporary_deformation_gradient =
          (static_cast<Real>(istep) + 1) / num_substep * _delta_deformation_gradient;
      _temporary_deformation_gradient += _temporary_deformation_gradient_old;

      solveQp();

      if (_convergence_failed)
      {
        if (_print_convergence_message)
          mooseWarning(
              "The crystal plasticity constitutive model has failed to converge. Increasing "
              "the number of substeps.");

        substep_iter++;
        num_substep *= 2;
        break;
      }
    }

    if (substep_iter > _max_substep_iter && _convergence_failed)
      mooseException("CompMultiCPStressFullCoupled: Constitutive failure");
  } while (_convergence_failed);

  postSolveQp(cauchy_stress, jacobian_mult);
}

void
CompMultiCPStressFullCoupled::preSolveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->setInitialConstitutiveVariableValues();

  (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr_old[_active_op_index])[_qp];
  _inverse_plastic_deformation_grad_old = (*_plastic_deformation_gradient_gr_old[_active_op_index])[_qp].inverse();

  if (_grain_id == FeatureFloodCount::invalid_id)
  {
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->resetInitialConstitutiveVariableValues();
  }
}

void
CompMultiCPStressFullCoupled::solveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setSubstepConstitutiveVariableValues();
    _models[i]->calculateSlipResistance();
  }

  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  solveStateVariables();
  if (_convergence_failed)
    return; // pop back up and take a smaller substep

  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->updateSubstepConstitutiveVariableValues();

  // save off the old F^{p} inverse now that have converged on the stress and state variables
  _inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;
}

void
CompMultiCPStressFullCoupled::postSolveQp(RankTwoTensor & cauchy_stress,
                                                    RankFourTensor & jacobian_mult)
{
  cauchy_stress = _elastic_deformation_gradient * (*_pk2_gr[_active_op_index])[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;

  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  (*_updated_rotation_gr[_active_op_index])[_qp] = rot * (*_crysrot_gr[_active_op_index])[_qp];
}

void
CompMultiCPStressFullCoupled::solveStateVariables()
{
  unsigned int iteration;
  bool iter_flag = true;

  iteration = 0;
  // Check for slip system resistance update tolerance
  do
  {
    solveStress();
    if (_convergence_failed)
      return;

    (*_plastic_deformation_gradient_gr[_active_op_index])[_qp] =
        _inverse_plastic_deformation_grad.inverse(); // the postSoveStress

    // Update slip system resistance and state variable after the stress has been finalized
    // We loop through all the models for each calculation
    // in order to make sure that when coupling appears, the state variables are updated based on
    // the same components
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->cacheStateVariablesBeforeUpdate();

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->calculateStateVariableEvolutionRateComponent();

    for (unsigned int i = 0; i < _num_models; ++i)
      if (!_models[i]->updateStateVariables())
        _convergence_failed = true;

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->calculateSlipResistance();

    if (_convergence_failed)
      return;

    for (unsigned int i = 0; i < _num_models; ++i)
    {
      // iter_flag = false, stop iteration if all values are converged and good to go
      // iter_flag = true, continue iteration if any value is not converged
      if (!_models[i]->areConstitutiveStateVariablesConverged())
      {
        iter_flag = true;
        break;
      }
      else
        iter_flag = false; // iter_flag = false, stop iteration only when all models returns true
    }

    if (iter_flag)
    {
      if (_print_convergence_message)
        mooseWarning("CompMultiCPStressFullCoupled: State variables (or the system "
                     "resistance) did not converge at element ",
                     _current_elem->id(),
                     " and qp ",
                     _qp,
                     "\n");
    }
    iteration++;
  } while (iter_flag && iteration < _maxiterg);

  if (iteration == _maxiterg)
  {
    if (_print_convergence_message)
      mooseWarning(
          "CompMultiCPStressFullCoupled: Hardness Integration error. Reached the "
          "maximum number of iterations to solve for the state variables at element ",
          _current_elem->id(),
          " and qp ",
          _qp,
          "\n");

    _convergence_failed = true;
  }
}

void
CompMultiCPStressFullCoupled::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidualAndJacobian();
  if (_convergence_failed)
  {
    if (_print_convergence_message)
      mooseWarning("CompMultiCPStressFullCoupled: the slip increment exceeds tolerance "
                   "at element ",
                   _current_elem->id(),
                   " and Gauss point ",
                   _qp);

    return;
  }

  rnorm = _residual_tensor.L2norm();
  rnorm0 = rnorm;

  // Check for stress residual tolerance; different from user object version which
  // compares the absolute tolerance of only the original rnorm value
  while (rnorm > _rtol * rnorm0 && rnorm > _abs_tol && iteration < _maxiter)
  {
    // Calculate stress increment
    dpk2 = -_jacobian.invSymm() * _residual_tensor;
    (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] + dpk2;

    calculateResidualAndJacobian();

    if (_convergence_failed)
    {
      if (_print_convergence_message)
        mooseWarning("CompMultiCPStressFullCoupled: the slip increment exceeds tolerance "
                     "at element ",
                     _current_elem->id(),
                     " and Gauss point ",
                     _qp);

      return;
    }

    rnorm_prev = rnorm;
    rnorm = _residual_tensor.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))
    {
      if (_print_convergence_message)
        mooseWarning("CompMultiCPStressFullCoupled: Failed with line search");

      _convergence_failed = true;
      return;
    }

    if (_use_line_search)
      rnorm = _residual_tensor.L2norm();

    iteration++;
  }

  if (iteration >= _maxiter)
  {
    if (_print_convergence_message)
      mooseWarning("CompMultiCPStressFullCoupled: Stress Integration error rmax = ",
                   rnorm,
                   " and the tolerance is ",
                   _rtol * rnorm0,
                   " when the rnorm0 value is ",
                   rnorm0,
                   " for element ",
                   _current_elem->id(),
                   " and qp ",
                   _qp);

    _convergence_failed = true;
  }
}

// Calculates stress residual equation and jacobian
void
CompMultiCPStressFullCoupled::calculateResidualAndJacobian()
{
  calculateResidual();
  if (_convergence_failed)
    return;
  calculateJacobian();
}

void
CompMultiCPStressFullCoupled::calculateResidual()
{
  RankTwoTensor ce, elastic_strain, ce_pk2, equivalent_slip_increment_per_model,
      equivalent_slip_increment, pk2_new;

  equivalent_slip_increment.zero();

  // calculate slip rate in order to compute F^{p-1}  
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    equivalent_slip_increment_per_model.zero();

    // calculate shear stress with consideration of contribution from other physics
    _models[i]->calculateShearStress(
        (*_pk2_gr[_active_op_index])[_qp], _inverse_eigenstrain_deformation_grad, _num_eigenstrains);

    _convergence_failed = !_models[i]->calculateSlipRate();

    if (_convergence_failed)
      return;

    _models[i]->calculateEquivalentSlipIncrement(equivalent_slip_increment_per_model);
    equivalent_slip_increment += equivalent_slip_increment_per_model;
  }
  
  RankTwoTensor residual_equivalent_slip_increment =
      RankTwoTensor::Identity() - equivalent_slip_increment;
  _inverse_plastic_deformation_grad =
      _inverse_plastic_deformation_grad_old * residual_equivalent_slip_increment;

  _elastic_deformation_gradient = _temporary_deformation_gradient *
                                  _inverse_eigenstrain_deformation_grad *
                                  _inverse_plastic_deformation_grad;

  ce = _elastic_deformation_gradient.transpose() * _elastic_deformation_gradient;

  elastic_strain = ce - RankTwoTensor::Identity();
  elastic_strain *= 0.5;
  _elastic_lagrangian_strain[_qp] = elastic_strain;

  pk2_new = (*_elasticity_tensor_gr[_active_op_index])[_qp] * elastic_strain;
  _residual_tensor = (*_pk2_gr[_active_op_index])[_qp] - pk2_new;
}

void
CompMultiCPStressFullCoupled::calculateJacobian()
{
  // may not need to cache the dfpinvdpk2 here. need to double check
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2, dfpinvdpk2_per_model;

  RankTwoTensor ffeiginv = _temporary_deformation_gradient * _inverse_eigenstrain_deformation_grad;

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
        dfedfpinv(i, j, k, j) = ffeiginv(i, k);

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->calculateTotalPlasticDeformationGradientDerivative(
        dfpinvdpk2_per_model,
        _inverse_plastic_deformation_grad_old,
        _inverse_eigenstrain_deformation_grad,
        _num_eigenstrains);
    dfpinvdpk2 += dfpinvdpk2_per_model;
  }

  _jacobian =
      RankFourTensor::IdentityFour() - ((*_elasticity_tensor_gr[_active_op_index])[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

void
CompMultiCPStressFullCoupled::calcTangentModuli(RankFourTensor & jacobian_mult)
{
  switch (_tan_mod_type)
  {
    case TangentModuliType::EXACT:
      elastoPlasticTangentModuli(jacobian_mult);
      break;
    default:
      elasticTangentModuli(jacobian_mult);
  }
}

void
CompMultiCPStressFullCoupled::elastoPlasticTangentModuli(RankFourTensor & jacobian_mult)
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2, feiginvfpinv;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  usingTensorIndices(i_, j_, k_, l_);
  dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                (*_elasticity_tensor_gr[_active_op_index])[_qp] * deedfe;

  pk2fet = (*_pk2_gr[_active_op_index])[_qp] * _elastic_deformation_gradient.transpose();
  fepk2 = _elastic_deformation_gradient * (*_pk2_gr[_active_op_index])[_qp];

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
      {
        tan_mod(i, j, i, l) += pk2fet(l, j);
        tan_mod(i, j, j, l) += fepk2(i, l);
      }

  tan_mod += dsigdpk2dfe;

  const auto je = _elastic_deformation_gradient.det();
  if (je > 0.0)
    tan_mod /= je;

  feiginvfpinv = _inverse_eigenstrain_deformation_grad * _inverse_plastic_deformation_grad;
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = feiginvfpinv(l, j);

  jacobian_mult = tan_mod * dfedf;
}

void
CompMultiCPStressFullCoupled::elasticTangentModuli(RankFourTensor & jacobian_mult)
{
  // update jacobian_mult
  jacobian_mult = (*_elasticity_tensor_gr[_active_op_index])[_qp];
}

bool
CompMultiCPStressFullCoupled::lineSearchUpdate(const Real & rnorm_prev,
                                                         const RankTwoTensor & dpk2)
{
  if (_line_search_method == LineSearchMethod::CutHalf)
  {
    Real rnorm;
    Real step = 1.0;

    do
    {
      (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] - step * dpk2;
      step /= 2.0;
      (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] + step * dpk2;

      calculateResidual();
      rnorm = _residual_tensor.L2norm();
    } while (rnorm > rnorm_prev && step > _min_line_search_step_size);

    // has norm improved or is the step still above minumum search step size?
    return (rnorm <= rnorm_prev || step > _min_line_search_step_size);
  }
  else if (_line_search_method == LineSearchMethod::Bisection)
  {
    unsigned int count = 0;
    Real step_a = 0.0;
    Real step_b = 1.0;
    Real step = 1.0;
    Real s_m = 1000.0;
    Real rnorm = 1000.0;

    calculateResidual();
    auto s_b = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm1 = _residual_tensor.L2norm();
    (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] - dpk2;
    calculateResidual();
    auto s_a = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm0 = _residual_tensor.L2norm();
    (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] + dpk2;

    if ((rnorm1 / rnorm0) < _line_search_tolerance || s_a * s_b > 0)
    {
      calculateResidual();
      return true;
    }

    while ((rnorm / rnorm0) > _line_search_tolerance && count < _line_search_max_iterations)
    {
      (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] - step * dpk2;
      step = 0.5 * (step_b + step_a);
      (*_pk2_gr[_active_op_index])[_qp] = (*_pk2_gr[_active_op_index])[_qp] + step * dpk2;
      calculateResidual();
      s_m = _residual_tensor.doubleContraction(dpk2);
      rnorm = _residual_tensor.L2norm();

      if (s_m * s_a < 0.0)
      {
        step_b = step;
        s_b = s_m;
      }
      if (s_m * s_b < 0.0)
      {
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    // below tolerance and max iterations?
    return ((rnorm / rnorm0) < _line_search_tolerance && count < _line_search_max_iterations);
  }
  else
    mooseError("Line search method is not provided.");
}

void
CompMultiCPStressFullCoupled::calculateEigenstrainDeformationGrad()
{
  _inverse_eigenstrain_deformation_grad.zero();
  _inverse_eigenstrain_deformation_grad.addIa(1.0);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setSubstepDt(_substep_dt);
    _eigenstrains[i]->computeQpProperties();
    _inverse_eigenstrain_deformation_grad *= _eigenstrains[i]->getDeformationGradientInverse();
  }
  (*_eigenstrain_deformation_gradient)[_qp] = _inverse_eigenstrain_deformation_grad.inverse();
}
