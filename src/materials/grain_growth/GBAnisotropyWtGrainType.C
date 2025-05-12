#include "GBAnisotropyWtGrainType.h"

registerMooseObject("mypandaApp", GBAnisotropyWtGrainType);

InputParameters
GBAnisotropyWtGrainType::validParams()
{
  InputParameters params = GBAnisotropyMisoriBase::validParams();
  params.addClassDescription(
    "Extends GBAnisotropyMisoriBase to account for grain-type-specific properties "
    "using data provided by an EBSD reader object.");
  params.addRequiredParam<UserObjectName>("ebsd_reader", 
    "Name of the EBSD reader user object that provides grain-specific data.");  
  params.addRequiredParam<UserObjectName>("grain_tracker", "GrainTracker user object");
  return params;
}

GBAnisotropyWtGrainType::GBAnisotropyWtGrainType(const InputParameters & parameters)
: GBAnisotropyMisoriBase(parameters),
  _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
  _grain_tracker(getUserObject<GrainTracker>("grain_tracker")),
  _grain_type(declareProperty<Real>("grain_type"))
{
}

void
GBAnisotropyWtGrainType::computeGBProperties()
{
  // Call the base class method to compute properties
  GBAnisotropyMisoriBase::computeGBProperties();

  _grain_type[_qp] = 0.0;
  Real sum_h = 0.0;

  // get the GB location based on the GrainTracker in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id()); 
  // Step 2: Compute effective dislocation density (rho_eff)
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Get order parameter value for current grain
    _grain_type[_qp] = _ebsd_reader.getAvgData(grain_id)._custom[0] * h;
    sum_h += h;
  }
  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _grain_type[_qp] /= sum_h;

  // If there is isotropy, use the HAGB values
  std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, _GBsigma_HAGB));
  std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, _GBmob_HAGB));
}