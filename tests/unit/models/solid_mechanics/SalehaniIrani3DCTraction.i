[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'traction'
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.1 0.12 -0.09'
  []
  [traction]
    type = Vec
    values = '7.2749570731 2.8239900883 -2.1179925662'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = 0.2
    tangential_gap_at_maximum_shear_traction = 0.3
    maximum_normal_traction = 10.0
    maximum_shear_traction = 8.0
  []
[]
