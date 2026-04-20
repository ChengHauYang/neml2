[Tensors]
  [jump]
    type = Vec
    values = '0.2 -0.1 0.3'
  []
  [effective_jump]
    type = Vec
    values = '0.2 -0.1224744871391589 0.36742346141747667'
  []
  [traction]
    type = Vec
    values = '2.1019315923815736 -1.0509657961907868 3.1528973885723603'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/interface_displacement_jump'
    input_Vec_values = 'jump'
    output_Scalar_names = 'state/internal/damage state/internal/effective_displacement_jump_scalar_max'
    output_Scalar_values = '0.6636909452189481 0.4358898943540675'
    output_Vec_names = 'state/internal/effective_jump state/interface_traction'
    output_Vec_values = 'effective_jump traction'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = 5.0
    softening_length = 0.4
    shear_weight = 1.5
    irreversible_damage = false
  []
[]
