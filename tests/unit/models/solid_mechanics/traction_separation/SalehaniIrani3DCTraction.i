[Tensors]
  [jump]
    type = Vec
    values = '0.05 0.03 -0.02'
  []
  [traction]
    type = Vec
    values = '4.9594295888429 1.8048276601915427 -1.2032184401276953'
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
    derivative_rel_tol = 1e-4
    derivative_abs_tol = 1e-4
    parameter_derivative_rel_tol = 1e-4
    parameter_derivative_abs_tol = 1e-4
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    # AD-backed Jacobian tracing captures differentiable temporaries as constants for this law.
    jit = false
    normal_gap_at_maximum_normal_traction = 0.2
    tangential_gap_at_maximum_shear_traction = 0.1
    maximum_normal_traction = 10
    maximum_shear_traction = 5
  []
[]
