[Tensors]
  [jump]
    type = Vec
    values = '0.1 -0.2 0.05'
  []
  [traction]
    type = Vec
    values = '1.0 -0.8 0.2'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/interface_traction'
    output_Vec_values = 'traction'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 10.0
    tangent_stiffness = 4.0
  []
[]
