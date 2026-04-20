[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump old_forces/interface_displacement_jump'
    input_Vec_values = 'jump jump_old'
    input_Scalar_names = 'old_state/internal/d forces/t old_forces/t'
    input_Scalar_values = 'd_old t t_old'
    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/d'
    output_Scalar_values = 'd_new'
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
  [jump_old]
    type = Vec
    values = '0.04 0.08 0.04'
  []
  [d_old]
    type = Scalar
    values = '0.1'
  []
  [t]
    type = Scalar
    values = '1.0'
  []
  [t_old]
    type = Scalar
    values = '0.0'
  []
  [traction]
    type = Vec
    values = '0 0 0'
  []
  [d_new]
    type = Scalar
    values = '0'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    normal_stiffness = 1000
    GI_c = 1.0
    GII_c = 2.0
    tensile_strength = 100
    shear_strength = 50
    exponent = 1.0
    viscosity = 0.0
    alpha = 1e-10
    jit = false
  []
[]
