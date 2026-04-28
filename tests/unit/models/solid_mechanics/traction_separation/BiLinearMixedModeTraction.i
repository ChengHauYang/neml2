[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'

    input_Vec_names = 'displacement_jump'
    input_Vec_values = 'jump'

    input_Scalar_names = 'damage~1'
    input_Scalar_values = '0.0'

    output_Vec_names = 'traction'
    output_Vec_values = 'T'

    output_Scalar_names = 'damage'
    output_Scalar_values = '0.5926889810380374'

    # FD against AD on a deeply branched forward needs a looser tolerance: the bilinear
    # damage law has a stiff slope (~140 per unit jump in some components) and forward FD
    # at h=1e-6 has visible truncation error.
    derivative_abs_tol = 1e-5
    derivative_rel_tol = 1e-3
    parameter_derivative_abs_tol = 1e-5
    parameter_derivative_rel_tol = 1e-3

    check_values = true
    check_derivatives = true
    # Second derivatives are not declared by the model
    check_second_derivatives = false
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    # AD-backed derivatives are incompatible with TorchScript tracing of the where()-heavy
    # forward graph (the tracer rejects grad-requiring tensors used as where-conditions).
    jit = false
    penalty_stiffness = 1.0e4
    mode_I_fracture_energy = 1.0
    mode_II_fracture_energy = 1.5
    normal_strength = 10.0
    shear_strength = 8.0
    power_law_exponent = 2.0
    criterion = 'BK'
  []
[]

[Tensors]
  # All components are well above the FD step (1e-6) so finite differencing of the
  # piecewise damage law is well-conditioned.
  [jump]
    type = Vec
    values = '2.0e-3 1.0e-3 5.0e-4'
  []
  [T]
    type = Vec
    values = '8.146220379239253 4.0731101896196265 2.0365550948098132'
  []
[]
