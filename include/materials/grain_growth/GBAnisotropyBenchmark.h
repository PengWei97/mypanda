//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations

/**
 * Function[kappa, gamma, m, L] = parameters (sigma, mob, w_GB, sigma0)
 * Parameter determination method is elaborated in Phys. Rev. B, 78(2), 024113, 2008, by N. Moelans
 * Thanks to Prof. Moelans for the explanation of her paper.
 */
class GBAnisotropyBenchmark : public Material
{
public:
  static InputParameters validParams();

  GBAnisotropyBenchmark(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  // Reads data related to Grain Boundary (GB) anisotropy from an input file and stores the values in the appropriate member variables.
  virtual void readGBAnisotropyData();

  // Computes and updates the parameters related to GB anisotropy, such as interface energy (sigma), mobility, and other physical constants.
  virtual void computeGBAnisotropyParameters();

  const unsigned int _mesh_dimension;

  const Real _wGB;
  const Real _length_scale;
  const Real _time_scale;
  const Real _delta_sigma;
  const Real _delta_mob;
  const Real _kb;
  const Real _JtoeV;

  const FileName _Anisotropic_GB_file_name;

  const bool _inclination_anisotropy;

  const VariableValue & _T;

  std::vector<std::vector<Real>> _sigma;
  std::vector<std::vector<Real>> _mob;
  std::vector<std::vector<Real>> _Q;
  std::vector<std::vector<Real>> _kappa_gamma;
  std::vector<std::vector<Real>> _a_g2;

  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _gamma;
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _mu;

  Real _mu_qp;

  const unsigned int _op_num;

  const std::vector<const VariableValue *> _vals;
  const std::vector<const VariableGradient *> _grad_vals;
};
