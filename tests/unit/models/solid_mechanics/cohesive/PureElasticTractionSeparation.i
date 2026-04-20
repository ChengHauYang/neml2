[Tensors]
  [delta]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [T_expected]
    type = Vec
    values = '10.0 10.0 15.0'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 100
    tangent_stiffness = 50
  []
[]
