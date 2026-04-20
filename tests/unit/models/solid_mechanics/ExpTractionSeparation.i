[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    input_Scalar_names = 'old_state/internal/effective_jump_max'
    input_Scalar_values = '0.4'
    output_Vec_names = 'forces/interface_traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/internal/damage state/internal/effective_jump_max'
    output_Scalar_values = '0.5506710359 0.4'
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.2 0.1 -0.05'
  []
  [traction]
    type = Vec
    values = '0.7189263426 0.3594631713 -0.1797315856'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 2.0
    softening_length = 0.5
    beta = 4.0
  []
[]
