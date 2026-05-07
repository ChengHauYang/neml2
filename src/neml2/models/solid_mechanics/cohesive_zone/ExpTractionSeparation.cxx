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

#include "neml2/models/solid_mechanics/cohesive_zone/ExpTractionSeparation.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() +=
      " Implements the reversible exponential softening form "
      "T = (1 - d) * (Gc / delta0^2) * delta with damage "
      "d = 1 - exp(-delta_eff / delta0). The MOOSE irreversibility branch "
      "(monotonically increasing delta_eff) is not yet wired here.";

  options.add_parameter<Scalar>("fracture_energy", "Fracture energy Gc.");
  options.add_parameter<Scalar>("characteristic_length", "Softening length scale delta0.");
  options.add_parameter<Scalar>("tangential_weight",
                                "Tangential weighting beta in the effective jump.");

  options.add_output("damage", "Scalar damage variable d in [0, 1).");

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "characteristic_length")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weight")),
    _d(declare_output_variable<Scalar>("damage"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // TODO: implement the forward operator from the spec, in this order:
  //   delta_eff = sqrt(delta_n^2 + beta * (delta_s1^2 + delta_s2^2) + eps)
  //   d         = 1 - exp(-delta_eff / delta0)
  //   T         = (1 - d) * (Gc / delta0^2) * delta
  // Use a small `eps` (e.g. 1e-16) to keep the sqrt differentiable at zero jump.
  if (out)
  {
    _T = _delta();
    // _d = ...;
  }

  // TODO: implement first derivatives of T (3x3) and d (3-vector) w.r.t. _delta.
  if (dout_din)
  {
    // _T.d(_delta) = ...;
    // _d.d(_delta) = ...;
  }
}
} // namespace neml2
