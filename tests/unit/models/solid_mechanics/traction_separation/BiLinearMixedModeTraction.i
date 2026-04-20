[Tensors]
  [jump]
    type = Vec
    values = '0.11 0.03 0.015'
  []
  [jump_old]
    type = Vec
    values = '0.09 0.02 0.01'
  []
  [traction]
    type = Vec
    values = '6.969977288219108 1.9009028967870294 0.9504514483935147'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    derivative_abs_tol = 1e-3
    input_Scalar_names = 'forces/t old_forces/t old_state/internal/damage'
    input_Scalar_values = '1.0 0.5 0.2'
    input_Vec_names = 'forces/interface_displacement_jump old_forces/interface_displacement_jump'
    input_Vec_values = 'jump jump_old'
    output_Scalar_names = 'state/internal/damage state/internal/mode_mixity state/internal/effective_jump_init state/internal/effective_jump_final state/internal/effective_jump'
    output_Scalar_values = '0.36636570107099015 0.30491836056815314 0.09769008779248471 0.16581502023350478 0.115'
    output_Vec_names = 'state/interface_traction'
    output_Vec_values = 'traction'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    penalty_stiffness = 100.0
    critical_energy_release_rate_I = 0.8
    critical_energy_release_rate_II = 1.2
    normal_strength = 10.0
    shear_strength = 8.0
    eta = 1.5
    viscosity = 0.0
    lag_mode_mixity = false
    lag_displacement_jump = false
  []
[]
