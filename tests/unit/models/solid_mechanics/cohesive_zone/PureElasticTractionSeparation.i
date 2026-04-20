[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [traction]
    type = Vec
    values = '10.0 10.0 15.0'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = '100.0'
    tangent_stiffness = '50.0'
  []
[]
