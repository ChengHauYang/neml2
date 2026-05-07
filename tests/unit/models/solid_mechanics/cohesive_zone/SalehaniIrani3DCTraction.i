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
    values = '0.001 0.0005 0.0'
  []
  [traction_expected]
    type = Vec
    values = '0.0 0.0 0.0'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    displacement_jump = 'forces/displacement_jump'
    traction = 'state/traction'
    normal_gap_at_maximum_normal_traction = 0.002
    tangential_gap_at_maximum_shear_traction = 0.002
    maximum_normal_traction = 50.0
    maximum_shear_traction = 30.0
  []
[]
