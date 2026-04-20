// Copyright 2024, UChicago Argonne, LLC
// All Rights Reserved
// Software Name: NEML2 -- the New Engineering material Model Library, version 2
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE (MIT)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "neml2/models/solid_mechanics/cohesive/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/stack.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() =
      "Salehani-Irani 3D coupled exponential traction-separation law. "
      "Normal component enters the exponent linearly (alpha=1); tangential components "
      "enter quadratically (alpha=2).";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Characteristic normal gap delta_n0 at which normal traction is maximised";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Characteristic tangential gap (raw input; effective gap is sqrt(2) times this value)";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Peak normal traction sigma_max";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Peak shear traction tau_max";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _u0n(declare_parameter<Scalar>("u0n", "normal_gap_at_maximum_normal_traction")),
    _u0s_raw(declare_parameter<Scalar>("u0s_raw", "tangential_gap_at_maximum_shear_traction")),
    _Tn_max(declare_parameter<Scalar>("Tn_max", "maximum_normal_traction")),
    _Ts_max(declare_parameter<Scalar>("Ts_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // Effective shear gap: delta_s0 = sqrt(2) * u0s_raw
  // Prefactors: a_n = e * Tn_max / u0n,  a_s = sqrt(e) * Ts_max / u0s_raw
  const auto scale_n = std::exp(1.0) * _Tn_max / _u0n;
  const auto scale_s = std::sqrt(std::exp(1.0)) * _Ts_max / _u0s_raw;

  const auto & delta = _jump();
  const auto delta_n = delta(0);
  const auto delta_s1 = delta(1);
  const auto delta_s2 = delta(2);

  // x = delta_n/u0n + (delta_s1^2 + delta_s2^2) / (2 * u0s_raw^2)
  const auto x = delta_n / _u0n +
                 (delta_s1 * delta_s1 + delta_s2 * delta_s2) / (2.0 * _u0s_raw * _u0s_raw);
  const auto exp_x = neml2::exp(-x);

  // Diagonal scaling matrix A = diag(scale_n, scale_s, scale_s)
  const auto A = R2::fill(scale_n, scale_s, scale_s);

  if (out)
    // T_i = a_i * delta_i * exp(-x)
    _traction = A * delta * exp_x;

  if (dout_din)
    if (_jump.is_dependent())
    {
      // J_ij = exp(-x) * (A_ij - (A*delta)_i * dx/ddelta_j)
      // dx/ddelta = [1/u0n, delta_s1/u0s_raw^2, delta_s2/u0s_raw^2]
      const auto dx_ddelta = Vec(base_stack(
          {Scalar::ones_like(delta_n) / _u0n, delta_s1 / (_u0s_raw * _u0s_raw),
           delta_s2 / (_u0s_raw * _u0s_raw)}));
      _traction.d(_jump) = exp_x * (A - outer(A * delta, dx_ddelta));
    }
}
} // namespace neml2
