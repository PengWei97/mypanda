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
  const Real getRhoInit(unsigned int grain_id) const;

  // Gets the time-evolved dislocation density for a specified grain.
  const Real getRhoWtTime(const unsigned int & grain_id, const Real & y_coord) const;
  const Real getRhoWtTime(const unsigned int & grain_id) const;
  const bool isType2Grain(const unsigned int & grain_id, const Real & y_coord) const;

protected:
  // Flag indicating if concurrent recovery is enabled.
  const bool _is_concurrent_recovery;
  const Real _rho_end1, _rho_end2;
  const Real _a_rho1, _a_rho2;
  const Real _rho_default;
  const bool _is_select_grains;

  std::unordered_map<unsigned int, Real> _grain_id_to_rho_init_map; // Example grain IDs for selection
};
