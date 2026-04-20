[Tensors]
  [delta]
    type = Vec
    values = '0.1 0.2 0.3'
  []
  [traction]
    type = Vec
    values = '0.3431 0.4162 0.6243'
  []
[]

[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Vec_names = 'forces/displacement_jump'
    input_Vec_values = 'delta'
    output_Vec_names = 'forces/traction'
    output_Vec_values = 'traction'
    value_rel_tol = 1e-3
  []
[]

[Models]
  [model]
    type = SalehaniIrani3DCTraction
    normal_gap_at_maximum_traction = 0.5
    tangential_gap_at_maximum_traction = 0.5
    maximum_normal_traction = 1.0
    maximum_shear_traction = 1.0
  []
[]
# Hand calculation check:
# delta_u0_n = 0.5
# delta_u0_s_scaled = sqrt(2) * 0.5 = 0.7071
# x = (0.1/0.5) + (0.2/0.7071)^2 + (0.3/0.7071)^2 = 0.2 + 0.08 + 0.18 = 0.46
# exp(-x) = exp(-0.46) = 0.63128
# a0 = exp(1)*1.0 = 2.71828
# a1 = sqrt(2*exp(1))*1.0 = 2.3316
# T0 = 2.71828 * (0.1/0.5) * 0.63128 = 2.71828 * 0.2 * 0.63128 = 0.3432 ?
# wait, my traction values in the test are different.
# Let's re-calculate more carefully.
# _t(0) = a0 * (delta_n / _delta_u0_n) * exp_neg_x = 2.71828 * 0.2 * 0.63128 = 0.3432
# _t(1) = a1 * (delta_s1 / delta_u0_s_scaled) * exp_neg_x = 2.3316 * (0.2/0.7071) * 0.63128 = 2.3316 * 0.2828 * 0.63128 = 0.4163
# _t(2) = a2 * (delta_s2 / delta_u0_s_scaled) * exp_neg_x = 2.3316 * (0.3/0.7071) * 0.63128 = 2.3316 * 0.4243 * 0.63128 = 0.6244
# I'll update the expected values to 0.3432 0.4163 0.6244
