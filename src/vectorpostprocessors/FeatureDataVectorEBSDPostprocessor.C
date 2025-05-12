#include "FeatureDataVectorEBSDPostprocessor.h"

registerMooseObject("mypandaApp", FeatureDataVectorEBSDPostprocessor);

InputParameters
FeatureDataVectorEBSDPostprocessor::validParams()
{
  InputParameters params = FeatureDataVectorPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("ebsd_reader",
                                          "The EBSDReader UserObject to get values from.");
  params.addClassDescription("This object is designed to pull information from the data structures "
                             "of a \"FeatureDataVectorPostprocessor\" or derived object (e.g. individual "
                             "feature volumes) and the EBSDReader object");
  return params;
}

FeatureDataVectorEBSDPostprocessor::FeatureDataVectorEBSDPostprocessor(
    const InputParameters & parameters)
  : FeatureDataVectorPostprocessor(parameters),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _grain_type(declareVector("grain_type"))
{
}

void
FeatureDataVectorEBSDPostprocessor::execute()
{
  FeatureDataVectorPostprocessor::execute();

  // Get the number of features
  const auto num_features = _feature_counter.getTotalFeatureCount();

  // Resize the grain type vector
  _grain_type.assign(num_features, -1);

  // Loop over each feature and assign the grain type
  for (MooseIndex(num_features) feature_num = 0; feature_num < num_features; ++feature_num)
  {
    _grain_type[feature_num] = _ebsd_reader.getAvgData(feature_num)._custom[0];
  }
}
