// GeneralDislocationRecovery.h
#pragma once

#include "GeneralUserObject.h"
#include "EBSDReaderMaterialProperty.h"

class GeneralDislocationRecovery : public GeneralUserObject
{
public:
  static InputParameters validParams();

  GeneralDislocationRecovery(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {};
  virtual void finalize() override {}

  Real getRhoWtTime(const unsigned int & grain_id, const Real & y_coord) const;
  Real getRhoWtTime(const unsigned int & grain_id) const;

protected:
  const bool _is_concurrent_recovery;
  const Real _rho_end1, _rho_end2;
  const Real _a_rho1, _a_rho2;
  const Real _rho_default;

  const EBSDReaderMaterialProperty & _GNDs_provider; // GNDs provider for EBSD reader
};
