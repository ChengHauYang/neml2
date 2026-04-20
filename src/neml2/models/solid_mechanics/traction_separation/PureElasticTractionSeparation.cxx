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
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/R2.h"

namespace neml2
{
register_NEML2_object(PureElasticTractionSeparation);

OptionSet
PureElasticTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Purely linear elastic traction-separation law.";

  options.set_parameter<TensorName<Scalar>>("normal_stiffness");
  options.set("normal_stiffness").doc() = "Elastic stiffness in the normal direction";

  options.set_parameter<TensorName<Scalar>>("tangent_stiffness");
  options.set("tangent_stiffness").doc() = "Elastic stiffness in the tangential direction";

  return options;
}

PureElasticTractionSeparation::PureElasticTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _kn(declare_parameter<Scalar>("Kn", "normal_stiffness")),
    _kt(declare_parameter<Scalar>("Kt", "tangent_stiffness"))
{
}

void
PureElasticTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  auto K = R2::fill(_kn, _kt, _kt);

  if (out)
    _traction = K * _displacement_jump;

  if (dout_din)
  {
    if (_displacement_jump.is_dependent())
      _traction.d(_displacement_jump) = K;

    if (const auto * const Kn = nl_param("Kn"))
      _traction.d(*Kn) =
          R2::fill(Scalar::fill(1.0, Kn->options()), Scalar::zeros_like(_kn), Scalar::zeros_like(_kn)) *
          _displacement_jump;

    if (const auto * const Kt = nl_param("Kt"))
      _traction.d(*Kt) =
          R2::fill(Scalar::zeros_like(_kt), Scalar::fill(1.0, Kt->options()), Scalar::fill(1.0, Kt->options())) *
          _displacement_jump;
  }
}
} // namespace neml2
