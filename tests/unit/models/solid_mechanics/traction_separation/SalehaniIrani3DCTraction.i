[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    check_derivatives = true
    check_values = false
    derivative_rel_tol = 1e-2
    derivative_abs_tol = 1e-2
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.05 0.1 0.05'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = 0.1
    tangential_gap_at_maximum_shear_traction = 0.2
    maximum_normal_traction = 100
    maximum_shear_traction = 50
  []
[]
