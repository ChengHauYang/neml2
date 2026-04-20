[Tensors]
  [delta]
    type = Vec
    values = '0.04 0.0 0.0'
  []
  [jump_zero]
    type = Vec
    values = '0.0 0.0 0.0'
  []
  [T_expected]
    type = Vec
    values = '40.0 0.0 0.0'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump old_state/internal/displacement_jump_store'
    input_Vec_values = 'delta jump_zero'
    input_Scalar_names = 'forces/time_step old_state/damage'
    input_Scalar_values = '1.0 0.0'
    output_Vec_names = 'state/traction state/internal/displacement_jump_store'
    output_Vec_values = 'T_expected delta'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = '0.0'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    penalty_stiffness = 1000
    normal_fracture_energy = 1.0
    shear_fracture_energy = 1.0
    tensile_strength = 100
    shear_strength = 50
    power_law_exponent = 2.0
    viscosity = 0.0
    criterion = 'BK'
    lag_mode_mixity = true
    lag_displacement_jump = false
  []
[]
