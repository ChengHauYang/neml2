[Tensors]
  [delta]
    type = Vec
    values = '0.5 0.5 0.5'
  []
  [T_expected]
    type = Vec
    values = '0.6420127083438707 0.38940039153570244 0.38940039153570244'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = 1.0
    tangential_gap_at_maximum_shear_traction = 1.0
    maximum_normal_traction = 1.0
    maximum_shear_traction = 1.0
  []
[]
