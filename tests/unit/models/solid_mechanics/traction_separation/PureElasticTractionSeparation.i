[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'u'

    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'T_expected'

    check_values = true
    check_derivatives = true
    check_second_derivatives = false
  []
[]

[Models]
  [model]
    type = PureElasticTractionSeparation
    normal_stiffness = '10.0'
    tangent_stiffness = '5.0'
  []
[]

[Tensors]
  [u]
    type = Vec
    values = '1.0 2.0 -3.0'
  []
  [T_expected]
    type = Vec
    values = '10.0 10.0 -15.0'
  []
[]
