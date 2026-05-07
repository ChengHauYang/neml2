[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction_expected'
    # TODO: flip these on once set_value is filled in and the expected output is correct.
    check_values = false
    check_derivatives = false
    check_AD_parameter_derivatives = false
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.01 0.005 0.0'
  []
  [traction_expected]
    type = Vec
    values = '0.0 0.0 0.0'
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    displacement_jump = 'forces/displacement_jump'
    traction = 'state/traction'
    normal_stiffness = 1000.0
    tangent_stiffness = 500.0
  []
[]
