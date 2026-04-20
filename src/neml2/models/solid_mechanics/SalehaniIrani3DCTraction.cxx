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

#include "neml2/models/solid_mechanics/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/exp.h"

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Three-dimensional exponential traction-separation law following Salehani and "
                  "Irani.";

  options.set<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal gap at the maximum normal traction";

  options.set<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Tangential gap at the maximum shear traction";

  options.set<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Maximum shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _delta_n0(declare_parameter<Scalar>("delta_n0", "normal_gap_at_maximum_normal_traction")),
    _delta_t0(declare_parameter<Scalar>("delta_t0", "tangential_gap_at_maximum_shear_traction")),
    _Tn_max(declare_parameter<Scalar>("Tn_max", "maximum_normal_traction")),
    _Ts_max(declare_parameter<Scalar>("Ts_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto jump = _jump();
  const auto sqrt2 = std::sqrt(2.0);
  const auto exp1 = std::exp(1.0);
  const Scalar delta_s0 = sqrt2 * _delta_t0;
  const Scalar b1 = jump(1) / delta_s0;
  const Scalar b2 = jump(2) / delta_s0;
  const Scalar x = jump(0) / _delta_n0 + b1 * b1 + b2 * b2;
  const Scalar expo = exp(-x);
  const Scalar Tn = exp1 * _Tn_max * jump(0) / _delta_n0 * expo;
  const Scalar Ts1 = std::sqrt(2.0 * exp1) * _Ts_max * jump(1) / delta_s0 * expo;
  const Scalar Ts2 = std::sqrt(2.0 * exp1) * _Ts_max * jump(2) / delta_s0 * expo;

  _traction = Vec(base_stack(std::vector<Tensor>{Tn, Ts1, Ts2}));
}
} // namespace neml2
