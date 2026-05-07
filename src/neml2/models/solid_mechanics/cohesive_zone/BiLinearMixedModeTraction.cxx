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

#include "neml2/models/solid_mechanics/cohesive_zone/BiLinearMixedModeTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() +=
      " Implements the Camanho-Davila bilinear mixed-mode law with a Benzeggagh-Kenane (BK) "
      "or power-law mixed-mode criterion. The current scaffolding handles only the reversible "
      "loading branch; irreversibility (history-based damage), mode-mixity / jump lagging, and "
      "viscous regularization are TODO.";

  options.add_parameter<Scalar>("penalty_stiffness", "Penalty elastic stiffness K.");
  options.add_parameter<Scalar>("GIc", "Mode I critical energy release rate.");
  options.add_parameter<Scalar>("GIIc", "Mode II critical energy release rate.");
  options.add_parameter<Scalar>("normal_strength", "Tensile (normal) strength N.");
  options.add_parameter<Scalar>("shear_strength", "Shear strength S.");
  options.add_parameter<Scalar>("criterion_exponent",
                                "Mixed-mode criterion exponent (BK or power law).");

  options.add_output("damage", "Scalar damage variable d in [0, 1].");

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "GIc")),
    _GIIc(declare_parameter<Scalar>("GIIc", "GIIc")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "criterion_exponent")),
    _d(declare_output_variable<Scalar>("damage"))
{
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // TODO: implement the bilinear law from the spec, in this order:
  //   1. mode mixity beta = sqrt(delta_s1^2 + delta_s2^2) / max(delta_n, eps_n)
  //   2. delta_init  from N, S, K, beta
  //   3. delta_final from GIc, GIIc, K, delta_init, beta (BK or POWER_LAW branch)
  //   4. effective mixed-mode jump delta_m using a regularized Macaulay split
  //   5. d_trial from the bilinear softening; clamp to [0, 1]
  //   6. T = (1 - d) * K * (delta_active) + K * (delta_inactive)
  // (The MOOSE source uses regularizedHeavyside(delta_n, alpha) for the
  // normal-component split; an `alpha` parameter and the irreversibility /
  // viscous-regularization branches will be added once they are designed.)
  if (out)
  {
    _T = _delta();
    // _d = ...;
  }

  // TODO: implement the 3x3 traction Jacobian and 3-vector damage Jacobian.
  if (dout_din)
  {
    // _T.d(_delta) = ...;
    // _d.d(_delta) = ...;
  }
}
} // namespace neml2
