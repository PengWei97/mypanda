#include "EBSDReaderMaterialProperty.h"

registerMooseObject("mypandaApp", EBSDReaderMaterialProperty);

InputParameters
EBSDReaderMaterialProperty::validParams()
{
  InputParameters params = EBSDReader::validParams();
  params.addClassDescription("EBSD reader for material property");
  params.addParam<bool>("is_concurrent_recovery", false, "Enable concurrent recovery model");
  params.addParam<Real>("a_rho1", 3.3e-3, "Evolution coefficient during medium time recovery");
  params.addParam<Real>("rho_end1", 2.10e12, "Dislocation density after long-term concurrent recovery");
  params.addParam<Real>("a_rho2", 3.3e-3, "Evolution coefficient during medium time recovery");
  params.addParam<Real>("rho_end2", 2.10e12, "Dislocation density after long-term concurrent recovery");
  params.addParam<Real>("rho_default", 2.0e15, "Dislocation density at the beginning of simulation");
  params.addParam<bool>("is_select_grains", false,  "Enable selection of specific grains for initial rho");
  params.addParam<std::vector<unsigned int>>("select_grain_ids", {99, 255, 347}, "Example grain IDs for selection");
  params.addParam<std::vector<Real>>("set_rho_init_vectors", {1.0e12, 1.0e12, 1.0e12}, "Initial dislocation density for selected grains");
  return params;
}

EBSDReaderMaterialProperty::EBSDReaderMaterialProperty(const InputParameters & parameters)
  : EBSDReader(parameters),
    _is_concurrent_recovery(getParam<bool>("is_concurrent_recovery")),
    _rho_end1(getParam<Real>("rho_end1")),
    _rho_end2(getParam<Real>("rho_end2")),
    _a_rho1(getParam<Real>("a_rho1")),
    _a_rho2(getParam<Real>("a_rho2")),
    _rho_default(getParam<Real>("rho_default")),
    _is_select_grains(getParam<bool>("is_select_grains"))
{
  const auto & ids = getParam<std::vector<unsigned int>>("select_grain_ids");
  const auto & rhos = getParam<std::vector<Real>>("set_rho_init_vectors");

  if (ids.size() != rhos.size())
    mooseError("select_grain_ids and set_rho_init_vectors must have the same length");
    
    for (std::size_t i = 0; i < ids.size(); ++i)
      _grain_id_to_rho_init_map[ids[i]] = rhos[i];
}

const Real
EBSDReaderMaterialProperty::getRhoInit(unsigned int grain_id) const
{
  // 优先使用手动指定的映射
  if (_is_select_grains)
  {
    auto it = _grain_id_to_rho_init_map.find(grain_id);
    if (it != _grain_id_to_rho_init_map.end())
      return it->second;
  }

  if (_custom_columns > 0 && grain_id < getGrainNum())
    return getAvgData(grain_id)._custom[0];

  return _rho_default;
}

const Real
EBSDReaderMaterialProperty::getRhoWtTime(const unsigned int & grain_id, const Real & y_coord) const
{
  const Real time = _fe_problem.time();
  const Real rho_init = getRhoInit(grain_id);

  if (!_is_concurrent_recovery)
    return rho_init;

  // 判断晶粒类型（根据y坐标）
  const bool is_type2 = (y_coord > 100.0 && rho_init > 1.0e13);
  const Real rho_end = is_type2 ? _rho_end2 : _rho_end1;
  const Real a_rho = is_type2 ? _a_rho2 : _a_rho1;

  Real rho_now = (rho_init > rho_end)
                   ? (rho_init - rho_end) * std::exp(-a_rho * time) + rho_end
                   : rho_end;

  return std::clamp(rho_now, rho_end, _rho_default);
}

const Real
EBSDReaderMaterialProperty::getRhoWtTime(const unsigned int & grain_id) const
{
  const Real time = _fe_problem.time();
  const Real rho_init = getRhoInit(grain_id);

  if (!_is_concurrent_recovery)
    return rho_init;

  Real rho_now = (rho_init > _rho_end1)
                   ? (rho_init - _rho_end1) * std::exp(-_a_rho1 * time) + _rho_end1
                   : _rho_end1;

  return std::clamp(rho_now, _rho_end1, _rho_default);
}

const bool 
EBSDReaderMaterialProperty::isType2Grain(const unsigned int & grain_id, const Real & y_coord) const
{ 
  const Real rho_init = getRhoInit(grain_id);
  return (y_coord > 100.0 && rho_init > 1.0e13);
}