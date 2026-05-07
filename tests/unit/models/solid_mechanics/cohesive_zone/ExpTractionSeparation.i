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
    values = '0.001 0.0005 0.0'
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
    type = ExpTractionSeparation
    displacement_jump = 'forces/displacement_jump'
    traction = 'state/traction'
    damage = 'state/damage'
    fracture_energy = 100.0
    characteristic_length = 0.001
    tangential_weight = 1.0
  []
[]
