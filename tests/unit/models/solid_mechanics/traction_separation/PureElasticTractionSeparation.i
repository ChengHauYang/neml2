[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'

    output_Vec_names = 'state/traction'
    output_Vec_values = 'T'

    check_values = true
    check_derivatives = true
    check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = 1000.0
    tangent_stiffness = 500.0
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.1 0.05 0.02'
  []
  [T]
    type = Vec
    values = '100.0 25.0 10.0'
  []
[]
