//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACSEDGPolyEBSD.h"

registerMooseObject("mypandaApp", ACSEDGPolyEBSD);

InputParameters
ACSEDGPolyEBSD::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Stored Energy contribution to grain growth");
  params.addRequiredParam<MaterialPropertyName>(
      "D_stored_energy_name", "The stored energy derivative for the specific order parameter");
  return params;
}

ACSEDGPolyEBSD::ACSEDGPolyEBSD(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_stored_energy(getMaterialProperty<Real>("D_stored_energy_name"))
{
}

Real
ACSEDGPolyEBSD::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
    {
      if (_fe_problem.time() > 10.0)
        return _u[_qp] * _D_stored_energy[_qp];
      else
        return 0.0;
    }

    case Jacobian:
      return 0.0;
  }
  mooseError("Invalid type passed in");
}