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

#include "neml2/models/solid_mechanics/traction_separation/PureElasticTractionSeparation.h"

namespace neml2
{
register_NEML2_object(PureElasticTractionSeparation);

OptionSet
PureElasticTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Purely elastic traction-separation law";

  options.set_parameter<TensorName<Scalar>>("normal_stiffness") = "1.0";
  options.set("normal_stiffness").doc() = "Normal stiffness";

  options.set_parameter<TensorName<Scalar>>("tangent_stiffness") = "1.0";
  options.set("tangent_stiffness").doc() = "Tangent stiffness";

  return options;
}

PureElasticTractionSeparation::PureElasticTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _kn(declare_parameter<Scalar>("normal_stiffness", "normal_stiffness", true)),
    _kt(declare_parameter<Scalar>("tangent_stiffness", "tangent_stiffness", true))
{
}

void
PureElasticTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  auto K = R2::fill(_kn, _kt, _kt);

  if (out)
    _traction = K * _jump();

  if (dout_din)
  {
    if (_jump.is_dependent())
      _traction.d(_jump) = K;

    if (auto p_kn = nl_param("normal_stiffness"))
      _traction.d(*p_kn) =
          R2::fill(Scalar::ones(_kn.options()), Scalar::zeros(_kn.options()), Scalar::zeros(_kn.options())) * _jump();

    if (auto p_kt = nl_param("tangent_stiffness"))
      _traction.d(*p_kt) =
          R2::fill(Scalar::zeros(_kt.options()), Scalar::ones(_kt.options()), Scalar::ones(_kt.options())) * _jump();
  }
}
} // namespace neml2
