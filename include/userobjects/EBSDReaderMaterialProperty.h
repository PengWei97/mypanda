#pragma once

#include "EBSDReader.h"

// Extension of EBSDReader to include dislocation density calculations.
class EBSDReaderMaterialProperty : public EBSDReader
{
public:
  // Defines the valid input parameters.
  static InputParameters validParams();

  // Constructor with input parameters.
  EBSDReaderMaterialProperty(const InputParameters & parameters);

  // Gets the initial dislocation density for a specified grain.
  Real getRhoInit(unsigned int grain_id) const;

  // Gets the time-evolved dislocation density for a specified grain.
  Real getRhoWtTime(unsigned int grain_id) const;

protected:
  // Flag indicating if concurrent recovery is enabled.
  const bool _is_concurrent_recovery;

  const Real _rho_end;
  const Real _a_rho;
  const Real _rho_default;
};
