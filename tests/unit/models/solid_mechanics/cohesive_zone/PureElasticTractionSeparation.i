[Tensors]
  [delta]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [traction]
    type = Vec
    values = '10 10 15'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    output_Vec_names = 'forces/traction'
    output_Vec_values = 'traction'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 100
    tangential_stiffness = 50
  []
[]
