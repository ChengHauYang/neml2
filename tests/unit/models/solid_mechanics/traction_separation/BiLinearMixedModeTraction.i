[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump old_forces/displacement_jump'
    input_Vec_values = 'jump_new jump_old'
    input_Scalar_names = 'old_state/damage forces/t old_forces/t'
    input_Scalar_values = '0.0 0.1 0.0'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = '0.724561'
    check_values = true
    check_derivatives = true
    # The smooth Heaviside regularization (alpha=1e-10) is far below the FD step (1e-6),
    # so the FD sees a near-discontinuity at delta_n=0. Loosen tolerances accordingly.
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-3
  []
[]

[Tensors]
  [jump_new]
    type = Vec
    values = '0.001 0.02 0.01'
  []
  [jump_old]
    type = Vec
    values = '-0.001 0.02 0.01'
  []
  # Lagged mixity: dn_old=-0.001 < 0 => has_opening=false => beta=0 (pure shear)
  # Separate jump_new/jump_old avoids dn_old=0 FD discontinuity in has_opening.
  # delta_init = S/K = 80/1e4 = 0.008
  # delta_final = sqrt(2)*2*GIIc/S = sqrt(2)*2*2/80 = 0.070711
  # delta_m = sqrt(0.001^2 + 0.02^2 + 0.01^2) = 0.022383 (in softening)
  # d_linear = 0.070711*(0.022383-0.008)/(0.022383*(0.070711-0.008)) = 0.724561
  # T_n = (1-0.724561)*1e4*0.001 = 2.75439
  # T_s1 = (1-0.724561)*1e4*0.02 = 55.0878, T_s2 = (1-0.724561)*1e4*0.01 = 27.5439
  [T_expected]
    type = Vec
    values = '2.75439 55.0878 27.5439'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    stiffness = '10000.0'
    mode_I_fracture_energy = '1.0'
    mode_II_fracture_energy = '2.0'
    normal_strength = '100.0'
    shear_strength = '80.0'
    criterion_exponent = '1.5'
    lag_mode_mixity = true
    lag_displacement_jump = false
  []
[]
