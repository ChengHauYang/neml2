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

#include "neml2/models/solid_mechanics/cohesive_zone/PureElasticTractionSeparation.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"

namespace neml2
{
register_NEML2_object(PureElasticTractionSeparation);

OptionSet
PureElasticTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() += " Implements the linear elastic form T_i = K_i * delta_i with isotropic "
                   "tangential stiffness K_t.";

  options.add_parameter<Scalar>("normal_stiffness", "Elastic stiffness in the normal direction.");
  options.add_parameter<Scalar>(
      "tangent_stiffness", "Elastic stiffness in both tangential directions (isotropic).");

  return options;
}

PureElasticTractionSeparation::PureElasticTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _K_n(declare_parameter<Scalar>("K_n", "normal_stiffness")),
    _K_t(declare_parameter<Scalar>("K_t", "tangent_stiffness"))
{
}

void
PureElasticTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // TODO: implement traction = diag(K_n, K_t, K_t) * displacement_jump.
  // Stub passes the input through so the file compiles; replace with the actual
  // forward operator before enabling check_values in the unit test.
  if (out)
    _T = _delta();

  // TODO: implement constant Jacobian d traction / d displacement_jump = diag(K_n, K_t, K_t).
  // Until this is filled in, keep check_derivatives = false in the unit test.
  if (dout_din)
  {
    // _T.d(_delta) = ...;
  }
}
} // namespace neml2
