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

    value_rel_tol = 1e-4
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = '0.5'
    tangential_gap_at_maximum_shear_traction = '0.5'
    maximum_normal_traction = '10.0'
    maximum_shear_traction = '5.0'
  []
[]

[Tensors]
  [u]
    type = Vec
    values = '0.1 0.1 0.1'
  []
  [T_expected]
    type = Vec
    values = '4.276641 1.296938 1.296938'
  []
[]
