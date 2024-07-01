#include "GrainTrackerMateProp.h"

registerMooseObject("mypandaApp", GrainTrackerMateProp);

InputParameters
GrainTrackerMateProp::validParams()
{
  InputParameters params = GrainTrackerElasticity::validParams();
  params.addClassDescription("Manages material properties based on grain_id");
  return params;
}

GrainTrackerMateProp::GrainTrackerMateProp(const InputParameters & parameters)
  : GrainTrackerElasticity(parameters),
  _angles_tmp(EulerAngles())
{
  _angles_tmp.random();
}

const EulerAngles &
GrainTrackerMateProp::updateEulerAngles(const unsigned int & grain_id) const
{
  if (grain_id < _euler.getGrainNum())
    return _euler.getEulerAngles(grain_id);

  return _angles_tmp;
}