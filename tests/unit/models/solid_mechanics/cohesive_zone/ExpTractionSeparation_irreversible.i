[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/displacement_jump'
    input_Vec_values = 'jump'
    input_Scalar_names = 'old_state/damage'
    input_Scalar_values = 'damage_old'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = 'damage'
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-8
    parameter_derivative_rel_tol = 1e-3
    parameter_derivative_abs_tol = 1e-8
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.005 0.008 0.006'
  []
  [damage_old]
    type = Scalar
    # Old damage > trial damage (~0.579) -> frozen at d_old=0.8
    values = '0.8'
  []
  [traction]
    type = Vec
    # frozen at d_old=0.8: T = (1-0.8)*5000*g = 1000*(0.005, 0.008, 0.006)
    values = '5.0 8.0 6.0'
  []
  [damage]
    type = Scalar
    values = '0.8'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = '0.5'
    softening_length_scale = '0.01'
    tangential_opening_weight = '0.5'
    irreversible = true
  []
[]
