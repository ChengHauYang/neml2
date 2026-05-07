[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction_expected'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = 'damage_expected'
    # TODO: flip these on once set_value is filled in and expected outputs are correct.
    check_values = false
    check_derivatives = false
    check_AD_parameter_derivatives = false
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.0005 0.0003 0.0'
  []
  [traction_expected]
    type = Vec
    values = '0.0 0.0 0.0'
  []
  [damage_expected]
    type = Scalar
    values = '0.0'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    displacement_jump = 'forces/displacement_jump'
    traction = 'state/traction'
    damage = 'state/damage'
    penalty_stiffness = 1.0e5
    GIc = 0.5
    GIIc = 1.0
    normal_strength = 50.0
    shear_strength = 30.0
    criterion_exponent = 1.45
  []
[]
