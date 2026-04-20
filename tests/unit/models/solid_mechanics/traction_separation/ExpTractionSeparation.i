[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    input_Scalar_names = 'old_state/internal/delta_eff_max'
    input_Scalar_values = 'zero'
    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/d state/internal/delta_eff_max'
    output_Scalar_values = 'd_val delta_eff_val'
    check_derivatives = true
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-3
    parameter_derivative_rel_tol = 1e-3
    parameter_derivative_abs_tol = 1e-3
  []
[]

[Tensors]
  [zero]
    type = Scalar
    values = '0'
  []
  [jump]
    type = Vec
    values = '0.1 0 0'
  []
  [traction]
    type = Vec
    values = '9.048374180359595 0 0'
  []
  [d_val]
    type = Scalar
    values = '0.09516258196404048'
  []
  [delta_eff_val]
    type = Scalar
    values = '0.1'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 100
    softening_length = 1
    tangential_weight = 1
    epsilon = 0
    jit = false
  []
[]
