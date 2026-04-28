[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'displacement_jump'
    input_Vec_values = 'jump'

    input_Scalar_names = 'effective_displacement_jump_max~1'
    input_Scalar_values = '0.005'

    output_Vec_names = 'traction'
    output_Vec_values = 'T'

    output_Scalar_names = 'effective_displacement_jump_max'
    output_Scalar_values = '0.022912878474781382'

    # FD truncation error on the sqrt-based effective jump is non-trivial; relax to match
    # the looser tolerance used by other smooth-but-coupled models.
    derivative_abs_tol = 1e-6
    derivative_rel_tol = 1e-4
    parameter_derivative_abs_tol = 1e-6
    parameter_derivative_rel_tol = 1e-4

    check_values = true
    check_derivatives = true
    # Second derivatives are not declared by the model
    check_second_derivatives = false
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 100.0
    softening_length_scale = 0.05
    tangential_weighting_factor = 1.0
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.02 0.01 0.005'
  []
  [T]
    type = Vec
    values = '505.907657802721 252.9538289013605 126.47691445068025'
  []
[]
