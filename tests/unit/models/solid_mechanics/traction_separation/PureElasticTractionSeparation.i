[Tensors]
  [jump]
    type = Vec
    values = '0.1 0.2 -0.3'
  []
  [traction]
    type = Vec
    values = '100 100 -150'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/interface/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/interface/traction'
    output_Vec_values = 'traction'
    check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 1000
    tangent_stiffness = 500
  []
[]
