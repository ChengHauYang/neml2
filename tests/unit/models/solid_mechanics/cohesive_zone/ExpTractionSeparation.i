[Tensors]
  [delta]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [delta_max_old]
    type = Scalar
    values = '0.0'
  []
  [traction]
    type = Vec
    values = '0.1892455 0.3784910 0.5677365'
  []
  [delta_max]
    type = Scalar
    values = '0.3741657'
  []
  [damage]
    type = Scalar
    values = '0.5268991'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    input_Scalar_names = 'old_state/internal/delta_max'
    input_Scalar_values = 'delta_max_old'
    output_Vec_names = 'forces/traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/delta_max state/internal/d'
    output_Scalar_values = 'delta_max damage'
    value_rel_tol = 1e-3
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 1.0
    softening_length = 0.5
    weighting_factor = 1.0
  []
[]
