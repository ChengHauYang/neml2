[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
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
    batch_shape = '(3,1)'
    values = "1.2 0.0 0.0
              0.0 1.2 0.0
              0.0 0.0 1.2"
  []
  [gc]
    type = Scalar
    values = '0.5 0.4 1.8'
    batch_shape = '(3,1)'
  []
  [beta]
    type = Scalar
    values = '0.5 2.0 1.0'
    batch_shape = '(3,1)'
  []
  [delta0]
    type = Scalar
    values = '0.01 0.05 0.1 0.15 0.3'
    batch_shape = '(1,5)'
  []
  [traction]
    type = Vec
    # flat_index = i * 5 + j
    batch_shape = '(3,5)'
    values = "0.0 0.0 0.0
              -9.0603258229293698e-09 0.0 0.0
              -3.6865274120057018e-04 0.0 0.0
              -8.9456700774004361e-03 0.0 0.0
              -1.2210425925822813e-01 0.0 0.0
              0.0 0.0 0.0
              0.0 -3.4106051316484804e-13 0.0
              0.0 -2.0465364336530452e-06 0.0
              0.0 -2.6036196962309077e-04 0.0
              0.0 -1.8631942808779556e-02 0.0
              0.0 0.0 0.0
              0.0 0.0 -3.2617172962545731e-08
              0.0 0.0 -1.3271498683220526e-03
              0.0 0.0 -3.2204412278641570e-02
              0.0 0.0 -4.3957533332962129e-01"
  []
  [damage]
    type = Scalar
    # flat_index = i * 5 + j
    # for i in a: for j in b:
    batch_shape = '(3,5)'
    values = '1.0 0.99999999996224864 0.99999385578764666 0.99966453737209748 0.98168436111126578
              1.0 0.99999999999999822 0.9999999573638243 0.99998779553267392 0.99650651072335383
              1.0 0.99999999996224864 0.99999385578764666 0.99966453737209748 0.98168436111126578'
  []
[]

[Models]
  [model]
    type = ExponentialDamageTractionLaw
    jump = "g"
    fracture_energy = 'gc'
    tangential_opening_weight = 'beta'
    softening_length_scale = 'delta0'
  []
[]
