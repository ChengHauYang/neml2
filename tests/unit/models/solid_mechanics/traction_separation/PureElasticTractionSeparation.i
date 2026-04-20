[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'forces/interface_traction'
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
    values = '100 40 60'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 1000
    tangent_stiffness = 200
  []
[]
