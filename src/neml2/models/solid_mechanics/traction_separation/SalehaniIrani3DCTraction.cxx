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
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/outer.h"
#include <cmath>

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential traction-separation law from Salehani & Irani (2018)";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction") = "1.0";
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal gap at maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction") = "1.0";
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Tangential gap at maximum shear traction";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction") = "1.0";
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction") = "1.0";
  options.set("maximum_shear_traction").doc() = "Maximum shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _delta_n0(declare_parameter<Scalar>("normal_gap_at_maximum_normal_traction", "normal_gap_at_maximum_normal_traction", true)),
    _delta_s0(declare_parameter<Scalar>("tangential_gap_at_maximum_shear_traction", "tangential_gap_at_maximum_shear_traction", true)),
    _T_n_max(declare_parameter<Scalar>("maximum_normal_traction", "maximum_normal_traction", true)),
    _T_s_max(declare_parameter<Scalar>("maximum_shear_traction", "maximum_shear_traction", true))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto dn = _jump()(0);
  const auto dt1 = _jump()(1);
  const auto dt2 = _jump()(2);

  const auto delta_u0_0 = _delta_n0;
  const auto delta_u0_1 = neml2::sqrt2 * _delta_s0;
  const auto delta_u0_2 = neml2::sqrt2 * _delta_s0;

  const auto x = (dn / delta_u0_0) + (dt1 / delta_u0_1) * (dt1 / delta_u0_1) +
                 (dt2 / delta_u0_2) * (dt2 / delta_u0_2);

  const auto exp_x = neml2::exp(-x);

  const auto exp_1 = std::exp(1.0);
  const auto a0 = exp_1 * _T_n_max;
  const auto a1 = std::sqrt(2.0 * exp_1) * _T_s_max;
  const auto a2 = std::sqrt(2.0 * exp_1) * _T_s_max;

  const auto b0 = dn / delta_u0_0;
  const auto b1 = dt1 / delta_u0_1;
  const auto b2 = dt2 / delta_u0_2;

  if (out)
  {
    _traction = Vec::fill(a0 * b0 * exp_x, a1 * b1 * exp_x, a2 * b2 * exp_x);
  }

  if (dout_din)
  {
    if (_jump.is_dependent())
    {
      const auto db0_du = Vec::fill(1.0 / delta_u0_0, Scalar::zeros(_jump.options()), Scalar::zeros(_jump.options()));
      const auto db1_du = Vec::fill(Scalar::zeros(_jump.options()), 1.0 / delta_u0_1, Scalar::zeros(_jump.options()));
      const auto db2_du = Vec::fill(Scalar::zeros(_jump.options()), Scalar::zeros(_jump.options()), 1.0 / delta_u0_2);

      const auto dx_du =
          Vec::fill(1.0 / delta_u0_0, 2.0 * dt1 / (delta_u0_1 * delta_u0_1), 2.0 * dt2 / (delta_u0_2 * delta_u0_2));

      _traction.d(_jump) =
          neml2::outer(Vec::fill(a0, Scalar::zeros(_jump.options()), Scalar::zeros(_jump.options())), exp_x * (db0_du - b0 * dx_du)) +
          neml2::outer(Vec::fill(Scalar::zeros(_jump.options()), a1, Scalar::zeros(_jump.options())), exp_x * (db1_du - b1 * dx_du)) +
          neml2::outer(Vec::fill(Scalar::zeros(_jump.options()), Scalar::zeros(_jump.options()), a2), exp_x * (db2_du - b2 * dx_du));
    }
  }
}
} // namespace neml2
