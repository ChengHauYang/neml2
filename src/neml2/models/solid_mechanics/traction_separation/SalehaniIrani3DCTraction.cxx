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

#include "neml2/models/solid_mechanics/traction_separation/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " following the 3D exponential law of Salehani & Irani (2018). Normal and tangential "
      "directions use exponents alpha=1 and alpha=2, respectively.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Characteristic normal gap at maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Characteristic tangential gap at maximum shear traction (internally scaled by sqrt(2))";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Peak normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Peak shear traction (applied to both tangential directions)";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _du0_n(declare_parameter<Scalar>("du0_n", "normal_gap_at_maximum_normal_traction")),
    _du0_t(declare_parameter<Scalar>("du0_t", "tangential_gap_at_maximum_shear_traction")),
    _T_max_n(declare_parameter<Scalar>("T_max_n", "maximum_normal_traction")),
    _T_max_s(declare_parameter<Scalar>("T_max_s", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // Scaled characteristic tangential gap: delta_u0_s = sqrt(2) * tangential_gap
  const auto du0_s = std::sqrt(2.0) * _du0_t;

  // Traction prefactors: a_n = e * T_max_n, a_s = sqrt(2e) * T_max_s
  const auto a_n = std::exp(1.0) * _T_max_n;
  const auto a_s = std::sqrt(2.0 * std::exp(1.0)) * _T_max_s;

  const auto & jump = _displacement_jump();
  const auto delta_n = jump(0);
  const auto delta_s1 = jump(1);
  const auto delta_s2 = jump(2);

  // Normalised components
  const auto b_n = delta_n / _du0_n;
  const auto b_s1 = delta_s1 / du0_s;
  const auto b_s2 = delta_s2 / du0_s;

  // Exponent: x = b_n + b_s1^2 + b_s2^2
  const auto x = b_n + b_s1 * b_s1 + b_s2 * b_s2;
  const auto exp_x = exp(-x);

  if (out)
    _traction = Vec::fill(a_n * b_n * exp_x, a_s * b_s1 * exp_x, a_s * b_s2 * exp_x);

  if (dout_din)
    if (_displacement_jump.is_dependent())
    {
      // T = Vec::fill(a_n*b_n, a_s*b_s1, a_s*b_s2) * exp_x
      // dT/d(delta) = diag(a_n/du0_n, a_s/du0_s, a_s/du0_s)*exp_x - outer(T, dx/ddelta)
      // where dx/d(delta_n)  = 1/du0_n
      //       dx/d(delta_si) = 2*delta_si / du0_s^2
      const auto T_val = Vec::fill(a_n * b_n * exp_x, a_s * b_s1 * exp_x, a_s * b_s2 * exp_x);
      const auto diag_J = R2::fill(a_n * exp_x / _du0_n,
                                   a_s * exp_x / du0_s,
                                   a_s * exp_x / du0_s);
      const auto dx_ddelta = Vec::fill(
          Scalar::ones_like(delta_n) / _du0_n,
          2.0 * delta_s1 / (du0_s * du0_s),
          2.0 * delta_s2 / (du0_s * du0_s));
      _traction.d(_displacement_jump) = diag_J + (-1.0) * outer(T_val, dx_ddelta);
    }
}
} // namespace neml2
