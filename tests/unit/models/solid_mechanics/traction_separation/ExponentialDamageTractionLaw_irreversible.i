[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'old_state/damage'
    input_Scalar_values = 'damage_old'
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
    batch_shape = '(3)'
    values = "1.2 0.0 0.0
              0.0 1.2 0.0
              0.0 0.0 1.2"
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
  [damage_old]
    type = Scalar
    batch_shape = '(3)'
    values = '0.0 9.9999999e-01 1.0'
  []
  [damage]
    type = Scalar
    batch_shape = '(3)'
    values = '1.0 9.9999999e-01 1.0'
  []
  [traction]
    type = Vec
    batch_shape = '(3)'
    values = "0.0 0.0 0.0
              0.0 -4.8e-07 0.0
              0.0 0.0 0.0"
  []
[]

[Models]
  [model]
    type = ExponentialDamageTractionLaw
    jump = "g"
    fracture_energy = 'gc'
    tangential_opening_weight = 'beta'
    softening_length_scale = 'delta0'
    irreversible = true
  []
[]
