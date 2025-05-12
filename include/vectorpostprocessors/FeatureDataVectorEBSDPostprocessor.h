#pragma once

#include "FeatureDataVectorPostprocessor.h"
#include "EBSDReader.h"

class FeatureDataVectorEBSDPostprocessor : public FeatureDataVectorPostprocessor
{
public:
  static InputParameters validParams();

  FeatureDataVectorEBSDPostprocessor(const InputParameters & parameters);

  virtual void execute() override;
protected:
  const EBSDReader & _ebsd_reader;

  VectorPostprocessorValue & _grain_type;
};

