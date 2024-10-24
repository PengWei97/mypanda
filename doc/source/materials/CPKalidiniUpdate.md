# 继承关系
CrystalPlasticityKalidindiUpdate > CrystalPlasticityStressUpdateBase > Material
1. CrystalPlasticityStressUpdateBase: 
   1. 函数 `setQp` ~ 将积分点索引赋予给 `CrystalPlasticityStressUpdateBase`;
   2. 函数 `initQpStatefulProperties` ~ 调用 `setMaterialVectorSize`
      1. 函数 `setMaterialVectorSize`: 初始化 `_tau, _flow_direction, _slip_resistance, _slip_increment`


# 计算流程
> 主要被 `ComputeMultiCPStress` 所调用来执行

1. `calculateFlowDirection(const RankTwoTensor & crysrot)`: 调用 `calculateSchmidTensor(_number_slip_systems, _slip_plane_normal, _slip_direction, _flow_direction[_qp], crysrot);`
2. `calculateSchmidTensor`: 输入 `number_slip_systems, plane_normal_vector, direction_vector, crysrot` 来计算 `schmid_tensor[_qp]` ~ 通过旋转矩阵旋转之后的施密特张量
3. 