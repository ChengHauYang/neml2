[Tensors]
  [jump]
    type = Vec
    values = '0.2 -0.1 0.05'
  []
  [traction]
    type = Vec
    values = '8.947649017926986 -1.2210802788870787 0.6105401394435394'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/interface_traction'
    output_Vec_values = 'traction'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = 0.3
    tangential_gap_at_maximum_shear_traction = 0.4
    maximum_normal_traction = 10.0
    maximum_shear_traction = 6.0
  []
[]
