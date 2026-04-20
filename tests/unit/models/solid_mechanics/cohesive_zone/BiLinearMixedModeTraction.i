[Tensors]
  [delta]
    type = Vec
    values = '0.001 0.001 0.001'
  []
  [d_old]
    type = Scalar
    values = '0.0'
  []
  [traction]
    type = Vec
    values = '1.0 1.0 1.0'
  []
  [damage]
    type = Scalar
    values = '0.0'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    input_Scalar_names = 'old_state/internal/d'
    input_Scalar_values = 'd_old'
    output_Vec_names = 'forces/traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/d'
    output_Scalar_values = 'damage'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    penalty_stiffness = 1000
    mode_I_critical_energy = 1.0
    mode_II_critical_energy = 1.0
    tensile_strength = 100
    shear_strength = 100
    power_law_exponent = 1.0
  []
[]
