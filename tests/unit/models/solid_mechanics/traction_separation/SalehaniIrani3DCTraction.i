[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'jump'
    output_Vec_names = 'state/traction'
    output_Vec_values = 'T_expected'
    check_values = true
    check_derivatives = true
    # Large second derivative (O(1e6)) means FD error O(h*f'') ≈ 1e-6*1e6 ~ 1 for step=1e-6.
    # Loosen derivative tolerance to accommodate FD inaccuracy.
    derivative_rel_tol = 1e-3
    derivative_abs_tol = 1e-3
    parameter_derivative_rel_tol = 1e-3
  []
[]

[Tensors]
  [jump]
    type = Vec
    values = '0.005 0.003 0.002'
  []
  # du0_s = sqrt(2)*0.01 = 0.014142
  # b_n=0.25, b_s1=0.21213, b_s2=0.14142
  # x=0.25+0.04500+0.02000=0.31500, exp_x=0.72979
  # a_n=e*200=543.656, a_s=sqrt(2e)*100=233.164
  # T_n=543.656*0.25*0.72979=99.189
  # T_s1=233.164*0.21213*0.72979=36.097
  # T_s2=233.164*0.14142*0.72979=24.064
  [T_expected]
    type = Vec
    values = '99.188592 36.096553 24.064369'
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction = '0.02'
    tangential_gap_at_maximum_shear_traction = '0.01'
    maximum_normal_traction = '200.0'
    maximum_shear_traction = '100.0'
  []
[]
