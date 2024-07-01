//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGrGrDeformationEnergyDrivingFC.h"

#include "Material.h"

registerMooseObject("mypandaApp", ACGrGrDeformationEnergyDrivingFC);

InputParameters
ACGrGrDeformationEnergyDrivingFC::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds elastic energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_deformation_energy_name", "The elastic/plastic energy derivative for the specific order parameter");
  return params;
}

ACGrGrDeformationEnergyDrivingFC::ACGrGrDeformationEnergyDrivingFC(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_deformation_energy(getMaterialProperty<Real>("D_deformation_energy_name"))
{
}

Real
ACGrGrDeformationEnergyDrivingFC::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return _D_deformation_energy[_qp]; // Compute the deformation energy driving force

    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}
