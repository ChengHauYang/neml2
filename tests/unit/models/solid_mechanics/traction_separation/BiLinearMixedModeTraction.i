[Tensors]
  [jump]
    type = Vec
    values = '0.02 0.01 0.005'
  []
  [traction]
    type = Vec
    values = '8.509568224393188 4.254784112196594 2.127392056098297'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/interface/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/interface/traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/interface_damage state/internal/mode_mixity state/internal/damage_initiation_jump state/internal/full_degradation_jump state/internal/effective_displacement_jump'
    output_Scalar_values = '0.5745215887803405 0.5590169943754439 0.010384889714355793 0.21496204790122736 0.022912878474784043'
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-3
    parameter_derivative_rel_tol = 1e-3
    parameter_derivative_abs_tol = 1e-3
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    # AD-backed Jacobian tracing captures differentiable temporaries as constants for this law.
    jit = false
    penalty_stiffness = 1000
    mode_I_critical_energy_release_rate = 1
    mode_II_critical_energy_release_rate = 2
    normal_strength = 10
    shear_strength = 12
    power_law_exponent = 1.5
    lag_mode_mixity = false
    regularization = 1e-3
  []
[]
