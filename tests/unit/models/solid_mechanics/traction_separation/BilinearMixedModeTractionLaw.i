[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/dt old_state/damage'
    input_Scalar_values = 'dt damage_old'
    input_Vec_names = 'g'
    input_Vec_values = 'g'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = 'damage'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
    check_derivatives = false
    check_AD_parameter_derivatives = false
  []
[]

[Tensors]
  [g]
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
    values = '0 0.5 1'
    # broadcast: damage_old(i, j) = damage_old(i, 0)
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
    # flat = i * 4 + j
    # (i=0) j=0 1 2 3
    # (i=1) j=0 1 2 3
    # (i=2) j=0 1 2 3
    batch_shape = '(3,4)'
    values = '0   0   0   0
              0.5 0.5 0.5 0.5
              1   1   1   1'
  []
  [traction]
    type = Vec
    batch_shape = '(3,4)'
    values = '40  20  0
              100 80  0
              -60 100 0
              160 80  0
              20  10  0
              50  40  0
              -60 50  0
              80  40  0
              0   0   0
              0   0   0
              -60 0   0
              0   0   0'
  []
[]

[Models]
  [model]
    type = BilinearMixedModeTractionLaw
    jump = 'g'
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
