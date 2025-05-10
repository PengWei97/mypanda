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

  _grain_type[_qp] = -1.0;

  // get the GB location based on the GrainTracker in the quadrature point
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id()); 

  // Extract grain IDs and corresponding indices
  for (MooseIndex(op_to_grains) i = 0; i < op_to_grains.size(); ++i)
  {
    if (op_to_grains[i] != FeatureFloodCount::invalid_id)
    {
      _grain_type[_qp] = _ebsd_reader.getAvgData(op_to_grains[i])._custom[0];
      break;
    }
  }

  // If there is isotropy, use the HAGB values
  std::fill(_sigma.begin(), _sigma.end(), std::vector<Real>(_op_num, _GBsigma_HAGB));
  std::fill(_mob.begin(), _mob.end(), std::vector<Real>(_op_num, _GBmob_HAGB));
}