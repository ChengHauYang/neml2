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

#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/stack.h"

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Salehani-Irani 3DC exponential traction-separation law with linear normal and "
                  "quadratic tangential coupling in the exponent.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Characteristic normal gap at maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Characteristic tangential gap at maximum shear traction";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Maximum tangential traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _delta_n0(declare_parameter<Scalar>("delta_n0", "normal_gap_at_maximum_normal_traction")),
    _delta_t0(declare_parameter<Scalar>("delta_t0", "tangential_gap_at_maximum_shear_traction")),
    _max_n(declare_parameter<Scalar>("max_n", "maximum_normal_traction")),
    _max_t(declare_parameter<Scalar>("max_t", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto jump = _jump();
  const auto sqrt2 = std::sqrt(2.0);
  const auto exp1 = std::exp(1.0);
  const auto delta_t = sqrt2 * _delta_t0;
  const auto a0 = exp1 * _max_n;
  const auto a1 = std::sqrt(2.0 * exp1) * _max_t;

  const auto b0 = jump(0) / _delta_n0;
  const auto b1 = jump(1) / delta_t;
  const auto b2 = jump(2) / delta_t;
  const auto x = b0 + b1 * b1 + b2 * b2;
  const auto exp_x = exp(-x);

  if (out)
    _traction = Vec::fill(a0 * b0 * exp_x, a1 * b1 * exp_x, a1 * b2 * exp_x);

  if (dout_din && _jump.is_dependent())
  {
    const auto dx0 = 1.0 / _delta_n0;
    const auto dx1 = 2.0 * jump(1) / (delta_t * delta_t);
    const auto dx2 = 2.0 * jump(2) / (delta_t * delta_t);

    const auto row0 = Vec::fill(
        a0 * exp_x * (dx0 - b0 * dx0), a0 * exp_x * (-b0 * dx1), a0 * exp_x * (-b0 * dx2));
    const auto row1 = Vec::fill(a1 * exp_x * (-b1 * dx0),
                                a1 * exp_x * (1.0 / delta_t - b1 * dx1),
                                a1 * exp_x * (-b1 * dx2));
    const auto row2 = Vec::fill(a1 * exp_x * (-b2 * dx0),
                                a1 * exp_x * (-b2 * dx1),
                                a1 * exp_x * (1.0 / delta_t - b2 * dx2));
    _traction.d(_jump) = R2(base_stack({row0, row1, row2}, -2));
  }
}
} // namespace neml2
