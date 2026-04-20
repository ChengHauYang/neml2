[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/dt old_state/damage'
    input_Scalar_values = 'dt damage_old'
    input_Vec_names = 'state/displacement_jump'
    input_Vec_values = 'jump'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = 'damage'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
    check_derivatives = true
    check_AD_parameter_derivatives = false
    derivative_rel_tol = 1e-4
    derivative_abs_tol = 1e-8
  []
[]

[Tensors]
  [jump]
    type = Vec
    batch_shape = '(1,4)'
    values = '0.02  0.01 0.0
              0.05  0.04 0.0
             -0.03  0.05 0.0
              0.08  0.04 0.0'
  []
  [dt]
    type = Scalar
    values = '1e-3'
  []
  [damage_old]
    type = Scalar
    values = '0.1 0.5 1'
    # broadcast: damage_old(i) = damage_old(i); all > d_trial≈0 so all frozen
    batch_shape = '(3,1)'
  []
  [K]
    type = Scalar
    values = '2e3'
  []
  [GIc]
    type = Scalar
    values = '1e3'
  []
  [GIIc]
    type = Scalar
    values = '1e2'
  []
  [N]
    type = Scalar
    values = '500'
  []
  [S]
    type = Scalar
    values = '300'
  []
  [eta]
    type = Scalar
    values = '1.2 2.2 4.0 6.0'
    batch_shape = '(1,4)'
  []
  [viscosity]
    type = Scalar
    values = '0 1e-4 1e-3 1e-2'
    batch_shape = '(1,4)'
  []
  [alpha]
    type = Scalar
    values = '1e-6'
  []
  [damage]
    type = Scalar
    # All test inputs are in elastic regime (delta_m << delta_init)
    # -> d_trial = 0, all frozen at d_old (0.1, 0.5, 1.0)
    batch_shape = '(3,4)'
    values = '0.1 0.1 0.1 0.1
              0.5 0.5 0.5 0.5
              1   1   1   1'
  []
  [traction]
    type = Vec
    # d_old=0.1: (1-0.1)*K=1800 on active, K=2000 on inactive
    # d_old=0.5: (1-0.5)*K=1000 on active, K=2000 on inactive
    # d_old=1.0: 0 on active, K=2000 on inactive (dn<0 only)
    batch_shape = '(3,4)'
    values = ' 36  18  0
               90  72  0
              -60  90  0
               144 72  0
               20  10  0
               50  40  0
              -60  50  0
               80  40  0
               0   0   0
               0   0   0
              -60  0   0
               0   0   0'
  []
[]

[Models]
  [model]
    type = BiLinearMixedModeTraction
    penalty_stiffness = 'K'
    GIc = 'GIc'
    GIIc = 'GIIc'
    normal_strength = 'N'
    shear_strength = 'S'
    eta = 'eta'
    viscosity = 'viscosity'
    alpha = 'alpha'
    lag_mode_mixity = false
    lag_displacement_jump = false
  []
[]
