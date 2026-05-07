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

#pragma once

#include "neml2/models/solid_mechanics/cohesive_zone/TractionSeparation.h"

namespace neml2
{
class Scalar;

/// Exponential softening cohesive law (reversible form):
///
///   delta_eff = sqrt(delta_n^2 + beta * (delta_s1^2 + delta_s2^2) + eps)
///   d         = 1 - exp(-delta_eff / delta0)
///   T         = (1 - d) * (Gc / delta0^2) * displacement_jump
///
/// The MOOSE source supports an `irreversible_damage` toggle that freezes
/// `delta_eff` at its historical maximum. That branch is intentionally NOT
/// scaffolded here: irreversibility in NEML2 wants either an additional
/// state variable for `max(delta_eff)` paired with a residual / time
/// integrator, or a hand-rolled stateful model. Add that wiring once the
/// design call is made.
class ExpTractionSeparation : public TractionSeparation
{
public:
  static OptionSet expected_options();

  ExpTractionSeparation(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Fracture energy.
  const Scalar & _Gc;

  /// Softening length scale.
  const Scalar & _delta0;

  /// Tangential weighting in the effective jump.
  const Scalar & _beta;

  /// Damage variable d in [0, 1).
  Variable<Scalar> & _d;
};
} // namespace neml2
