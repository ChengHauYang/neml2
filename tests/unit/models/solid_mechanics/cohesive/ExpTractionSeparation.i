[Tensors]
  [delta]
    type = Vec
    values = '0.5 0.5 0.5'
  []
  [T_expected]
    type = Vec
    values = '0.2103100130270574 0.2103100130270574 0.2103100130270574'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    input_Scalar_names = 'old_state/internal/delta_eff_max'
    input_Scalar_values = '0.0'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
    output_Scalar_names = 'state/internal/delta_eff_max'
    output_Scalar_values = '0.8660254037844386'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 1.0
    softening_length = 1.0
    tangential_weight = 1.0
    irreversible_damage = true
  []
[]
