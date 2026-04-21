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

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Salehani-Irani three-dimensional exponential cohesive traction law.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal gap at maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Tangential gap at maximum shear traction";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Maximum shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _delta_n0(declare_parameter<Scalar>("delta_n0", "normal_gap_at_maximum_normal_traction")),
    _delta_t0(declare_parameter<Scalar>("delta_t0", "tangential_gap_at_maximum_shear_traction")),
    _T_n_max(declare_parameter<Scalar>("T_n_max", "maximum_normal_traction")),
    _T_s_max(declare_parameter<Scalar>("T_s_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::request_AD()
{
  _traction.request_AD(_jump);
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto delta_s0 = std::sqrt(2.0) * _delta_t0;
  const auto x = _jump()(0) / _delta_n0 + _jump()(1) * _jump()(1) / (delta_s0 * delta_s0) +
                 _jump()(2) * _jump()(2) / (delta_s0 * delta_s0);
  const auto decay = exp(-x);

  _traction = Vec::fill(std::exp(1.0) * _T_n_max * _jump()(0) / _delta_n0 * decay,
                        std::sqrt(2.0 * std::exp(1.0)) * _T_s_max * _jump()(1) / delta_s0 * decay,
                        std::sqrt(2.0 * std::exp(1.0)) * _T_s_max * _jump()(2) / delta_s0 * decay);
}
} // namespace neml2
