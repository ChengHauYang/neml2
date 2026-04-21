[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'u'

    input_Scalar_names = 'old_state/internal/delta_eff_max'
    input_Scalar_values = '0.0'

    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'T_expected'

    output_Scalar_names = 'state/internal/delta_eff_max'
    output_Scalar_values = '0.27386127'

    check_values = true
    check_derivatives = true
    check_second_derivatives = false
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = '1.0'
    softening_length = '0.5'
    weighting_factor = '0.5'
    eps = '0.0'
    irreversible_damage = false
  []
[]

[Tensors]
  [u]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [T_expected]
    type = Vec
    values = '0.2313054 0.4626108 0.6939162'
  []
[]
