[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'

    output_Vec_names = 'state/traction'
    output_Vec_values = 'T'

    # FD truncation error is large for the exponentially coupled traction; relax to match
    # the looser tolerance used by other exponential/coupled models (e.g. GeneralElasticity.i).
    derivative_abs_tol = 1e-6
    derivative_rel_tol = 1e-4
    parameter_derivative_abs_tol = 1e-6
    parameter_derivative_rel_tol = 1e-4

    check_values = true
    check_derivatives = true
    # Second derivatives are not declared by the model
    check_second_derivatives = false
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = 0.1
    tangential_gap_at_maximum_shear_traction = 0.2
    maximum_normal_traction = 50.0
    maximum_shear_traction = 30.0
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.05 0.03 0.04'
  []
  [T]
    type = Vec
    values = '39.94988624876583 4.361549555143547 5.815399406858063'
  []
[]
