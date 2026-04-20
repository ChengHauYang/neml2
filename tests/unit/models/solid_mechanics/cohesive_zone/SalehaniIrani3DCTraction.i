[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
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
  [traction]
    type = Vec
    # Computed: delta_n0=0.01, delta_t0=0.02 -> ds0=sqrt(2)*0.02
    # a_n=e*100, a_s=sqrt(2e)*80; g=(0.005,0.008,0.006)
    # x = 0.5 + 0.08^2/(0.02828^2) + 0.06^2/(0.02828^2) = 0.625
    values = '72.74957073 28.23990088 21.17992566'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = '0.01'
    tangential_gap_at_maximum_shear_traction = '0.02'
    maximum_normal_traction = '100.0'
    maximum_shear_traction = '80.0'
  []
[]
