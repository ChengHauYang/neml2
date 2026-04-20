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

#include "neml2/models/solid_mechanics/cohesive_zone/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/tensors.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/sqrt.h"

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential traction separation law from Salehani & Irani (2018).";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_traction");
  options.set("normal_gap_at_maximum_traction").doc() = "Normal gap at maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_traction");
  options.set("tangential_gap_at_maximum_traction").doc() = "Tangential gap at maximum shear traction";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Maximum shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _delta_u0_n(declare_parameter<Scalar>("delta_u0_n", "normal_gap_at_maximum_traction", false)),
    _delta_u0_s(declare_parameter<Scalar>("delta_u0_s", "tangential_gap_at_maximum_traction", false)),
    _max_tn(declare_parameter<Scalar>("max_tn", "maximum_normal_traction", false)),
    _max_ts(declare_parameter<Scalar>("max_ts", "maximum_shear_traction", false))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, [[maybe_unused]] bool d2out_din2)
{
  auto delta = _delta();
  auto delta_n = delta(0);
  auto delta_s1 = delta(1);
  auto delta_s2 = delta(2);

  // Scaling based on spec
  auto delta_u0_s_scaled = Scalar::full(std::sqrt(2.0)) * _delta_u0_s;

  auto x = (delta_n / _delta_u0_n) + (delta_s1 / delta_u0_s_scaled) * (delta_s1 / delta_u0_s_scaled) +
           (delta_s2 / delta_u0_s_scaled) * (delta_s2 / delta_u0_s_scaled);

  auto exp_neg_x = Scalar(exp(-x));

  auto a0 = Scalar::full(std::exp(1.0)) * _max_tn;
  auto a1 = Scalar::full(std::sqrt(2.0 * std::exp(1.0))) * _max_ts;
  auto a2 = a1;

  if (out)
  {
    _t = Tensor(Vec::fill(a0 * (delta_n / _delta_u0_n) * exp_neg_x,
                          a1 * (delta_s1 / delta_u0_s_scaled) * exp_neg_x,
                          a2 * (delta_s2 / delta_u0_s_scaled) * exp_neg_x));
  }

  if (dout_din)
  {
    if (_delta.is_dependent())
    {
      auto dx_ddelta = Vec::fill(1.0 / _delta_u0_n,
                                 2.0 * delta_s1 / (delta_u0_s_scaled * delta_u0_s_scaled),
                                 2.0 * delta_s2 / (delta_u0_s_scaled * delta_u0_s_scaled));

      auto b = Vec::fill(delta_n / _delta_u0_n,
                         delta_s1 / delta_u0_s_scaled,
                         delta_s2 / delta_u0_s_scaled);
      auto a = Vec::fill(a0, a1, a2);

      // dTi_duj = a_i * exp_x * (dbi_duj - b_i * dx_duj)
      // dbi_duj = diag(1/delta_u0_i)
      R2 db_du = R2::fill(1.0 / _delta_u0_n, 1.0 / delta_u0_s_scaled, 1.0 / delta_u0_s_scaled);

      _t.d(_delta) =
          Tensor(R2::fill(a(0) * exp_neg_x, a(1) * exp_neg_x, a(2) * exp_neg_x) *
                 (db_du - R2(outer(b, dx_ddelta))));
    }
  }
}
} // namespace neml2
