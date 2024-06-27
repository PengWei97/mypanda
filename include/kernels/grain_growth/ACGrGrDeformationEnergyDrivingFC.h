//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ACBulk.h"

/**
 * Calculates the porton of the Allen-Cahn equation that results from the elastic deformation energy.
 */
class ACGrGrDeformationEnergyDrivingFC : public ACBulk<Real>
{
public:
  static InputParameters validParams();

  ACGrGrDeformationEnergyDrivingFC(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);

private:
  const MaterialProperty<Real> & _D_deformation_energy;
};
