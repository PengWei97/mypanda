#pragma once

#include "GBAnisotropyMisoriBase.h"
#include "GrainTracker.h"
#include "EBSDReader.h"

/**
 * @class GBAnisotropyWtGrainType
 * @brief Extends GBAnisotropyMisoriBase to account for grain-type-specific properties
 *        during the computation of grain boundary (GB) energy and mobility.
 *
 * This class introduces the concept of grain type during the computation
 * of GB properties, leveraging information from the EBSDReader to adjust
 * the GB characteristics accordingly.
 */

class GBAnisotropyWtGrainType : public GBAnisotropyMisoriBase
{
public:
  static InputParameters validParams();

  GBAnisotropyWtGrainType(const InputParameters & parameters);

protected:
  virtual void computeGBProperties() override;

  /// References to user objects for tracking grain information.
  const GrainTracker & _grain_tracker;
  const EBSDReader & _ebsd_reader;

  /// Material property that stores the type of grain at each quadrature point.
  MaterialProperty<Real> & _grain_type;
};