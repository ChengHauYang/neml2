[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/jump'
    input_Scalar_values = 'g'
    output_Scalar_names = 'state/traction'
    output_Scalar_values = 'traction'
    check_AD_parameter_derivatives = false
    derivative_abs_tol = 1e-7
  []
[]

[Tensors]
  [g]
    type = Vec
    batch_shape = '(3)'
    values = "1.2 0.0 0.0
              0.0 1.2 0.0
              0.0 0.0 1.2"
  []
  [n]
    type = Vec
    batch_shape = '(3)'
    values = "0.70710678 0.70710678 0.0
              0.0        0.70710678 0.70710678
              0.70710678 0.0        0.70710678"
  []
  [gc]
    type = Scalar
    values = '0.5 0.4 1.8'
    batch_shape = '(3)'
  []
  [beta]
    type = Scalar
    values = '0.5 2.0 1.0'
    batch_shape = '(3)'
  []
  [delta0]
    type = Scalar
    values = '0.01 0.1 0.15'
    batch_shape = '(3)'
  []
[]

[Models]
  [model]
    type = ExponentialDamageTractionLaw
    jump = "g"
    normal = 'n'
    fracture_energy = 'gc'
    tangential_opening_weight = 'beta'
    softening_length_scale = 'delta0'
  []
[]
