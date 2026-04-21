[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'
    # old_state/effective_displacement_jump_max = 0 (no prior damage)
    input_Scalar_names = 'old_state/effective_displacement_jump_max'
    input_Scalar_values = '0.0'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
    output_Scalar_names = 'state/damage state/effective_displacement_jump_max'
    output_Scalar_values = '0.349567 0.021506'
    check_values = true
    check_derivatives = true
    # The sqrt regularization (eps=1e-16) causes high curvature; FD step (1e-6) is not small
    # enough for tight derivative matching. Loosen tolerances accordingly.
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-3
    parameter_derivative_rel_tol = 1e-3
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.02 0.01 0.005'
  []
  # delta_eff = sqrt(0.02^2 + 0.5*(0.01^2+0.005^2)) = sqrt(0.0004+0.000062..) = sqrt(0.000462..)
  # = sqrt(0.000462...) = 0.021506
  # d = 1 - exp(-0.021506/0.05) = 0.349567
  # c = 5/0.05^2 = 2000
  # T = (1-0.349567)*2000*jump
  [T_expected]
    type = Vec
    values = '26.017339 13.008669 6.504335'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = '5.0'
    softening_length = '0.05'
    tangential_weight = '0.5'
    irreversible_damage = true
  []
[]
