[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump old_forces/interface_displacement_jump'
    input_Vec_values = 'jump jump_old'
    input_Scalar_names = 'old_state/internal/damage forces/t old_forces/t'
    input_Scalar_values = '0 1 0'
    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/damage state/internal/beta state/internal/delta_init state/internal/delta_final state/internal/delta_m'
    output_Scalar_values = '0.8571628019 0.4 0.01 0.2 0.0538516481'
    derivative_abs_tol = 5e-3
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.05 0.02 0'
  []
  [jump_old]
    type = Vec
    values = '0.04 0.015 0'
  []
  [traction]
    type = Vec
    values = '7.1418599041 2.8567439616 0'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTractionSeparation
    criterion = BK
    penalty_stiffness = 1000.0
    mode_I_fracture_energy = 1.0
    mode_II_fracture_energy = 1.0
    normal_strength = 10.0
    shear_strength = 10.0
    power_law_exponent = 1.0
    viscosity = 0.0
    mode_mixity = 'state/internal/beta'
    critical_displacement_jump = 'state/internal/delta_init'
    final_displacement_jump = 'state/internal/delta_final'
    effective_displacement_jump = 'state/internal/delta_m'
    lag_mode_mixity = false
  []
[]
