[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'state/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'traction'
    output_Scalar_names = 'state/damage'
    output_Scalar_values = 'damage'
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-8
    parameter_derivative_rel_tol = 1e-3
    parameter_derivative_abs_tol = 1e-8
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.005 0.008 0.006'
  []
  [traction]
    type = Vec
    # Gc=0.5, delta0=0.01, beta=0.5; g=(0.005,0.008,0.006)
    # delta_eff = sqrt(0.005^2 + 0.5*(0.008^2+0.006^2)) = sqrt(7.5e-5) = 0.008660254
    # d = 0.57937997, c=5000; T = (1-d)*c*g
    values = '10.51550065 16.82480104 12.61860078'
  []
  [damage]
    type = Scalar
    values = '0.57937997'
  []
[]

[Models]
  [model]
    type = ExpTractionSeparation
    fracture_energy = '0.5'
    softening_length_scale = '0.01'
    tangential_opening_weight = '0.5'
  []
[]
