[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'g n'
    input_Vec_values = 'g n'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
    check_derivatives = false
    check_AD_parameter_derivatives = false
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
  [traction]
    type = Vec
    batch_shape = '(3)'
    values = "0.0 0.0 0.0
              0.0 -1.98811713e-05 0.0
              0.0 0.0 -3.22044125e-02"
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
