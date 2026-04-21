[Tensors]
  [jump]
    type = Vec
    values = '0.1 0.05 -0.02'
  []
  [traction]
    type = Vec
    values = '0.5962577014123308 0.2981288507061654 -0.11925154028246615'
  []
  [effective_jump]
    type = Vec
    values = '0.1 0.1 -0.04'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/interface/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/interface/traction state/interface/effective_displacement_jump'
    output_Vec_values = 'traction effective_jump'
    output_Scalar_names = 'state/internal/effective_displacement_jump_max state/internal/interface_damage'
    output_Scalar_values = '0.14696938456699105 0.25467787323458657'
    derivative_rel_tol = 1e-4
    derivative_abs_tol = 1e-4
    parameter_derivative_rel_tol = 1e-4
    parameter_derivative_abs_tol = 1e-4
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    # AD-backed Jacobian tracing captures differentiable temporaries as constants for this law.
    jit = false
    critical_energy_release_rate = 2
    softening_length = 0.5
    tangential_weight = 4
  []
[]
