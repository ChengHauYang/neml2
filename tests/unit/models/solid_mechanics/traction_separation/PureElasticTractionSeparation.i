[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
    check_values = true
    check_derivatives = true
    check_second_derivatives = true
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.01 0.005 0.003'
  []
  # T_n = 100*0.01 = 1.0, T_s1 = 50*0.005 = 0.25, T_s2 = 50*0.003 = 0.15
  [T_expected]
    type = Vec
    values = '1.0 0.25 0.15'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = '100.0'
    tangent_stiffness = '50.0'
  []
[]
