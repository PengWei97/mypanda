//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Action that sets up ACSEDGPolyEBSD Kernels (based ACSEDGPoly) that adds the stored energy contribution
 * to grain growth models. This allows such models to simulate recrystallization as well.
 * The input GNDs are based on EBSD data.
 */
class PolyStoredEnergyEBSDAction : public Action,
                                   public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  PolyStoredEnergyEBSDAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of grains to create
  const unsigned int _op_num;

  /// base name for the order parameter variables
  const std::string _var_name_base;

  const std::string _stored_energy_name;
};
