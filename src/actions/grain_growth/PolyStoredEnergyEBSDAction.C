//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PolyStoredEnergyEBSDAction.h"
#include "Factory.h"
#include "Conversion.h"
#include "FEProblem.h"

registerMooseAction("mypandaApp", PolyStoredEnergyEBSDAction, "add_kernel");

InputParameters
PolyStoredEnergyEBSDAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Action that adds the contribution of stored energy associated with dislocations to grain growth models");
  params.addRequiredParam<unsigned int>("op_num",
                                        "specifies the total number of OPs representing");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");

  params.addParam<bool>(
      "use_displaced_mesh", false, "Whether to use displaced mesh in the kernels");
  return params;
}

PolyStoredEnergyEBSDAction::PolyStoredEnergyEBSDAction(const InputParameters & params)
  : Action(params),
    _op_num(getParam<unsigned int>("op_num")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _stored_energy_name("stored_energy")
{
}

void
PolyStoredEnergyEBSDAction::act()
{
  for (unsigned int op = 0; op < _op_num; ++op)
  {
    // Create variable names
    std::string var_name = _var_name_base + Moose::stringify(op);

    // Create stored energy derivative name
    MaterialPropertyName D_stored_energy_name = 
        derivativePropertyNameFirst(_stored_energy_name, var_name);

    // Step up ACSEDGPolyEBSD -- Allen-Cahn Stored Energy in Deformed Grains kernels using EBSD
    InputParameters params = _factory.getValidParams("ACSEDGPolyEBSD");

    // Set the actual parameters for the kernel
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<MaterialPropertyName>("D_stored_energy_name") = D_stored_energy_name;

    std::string kernel_name = "ACStoredEnergy_" + var_name;
    _problem->addKernel("ACSEDGPolyEBSD", kernel_name, params);
  }
}
