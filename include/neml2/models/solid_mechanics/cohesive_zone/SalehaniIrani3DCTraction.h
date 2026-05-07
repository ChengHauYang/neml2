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

/// Salehani-Irani 3D coupled exponential traction-separation law. The normal
/// component enters linearly (alpha=1) and the two tangential components
/// enter quadratically (alpha=2) in a single coupled exponential envelope:
///
///   x = (delta_n / d0_n) + (delta_s1 / d0_s)^2 + (delta_s2 / d0_s)^2
///   T_n  = e        * Tmax_n * (delta_n  / d0_n) * exp(-x)
///   T_si = sqrt(2e) * Tmax_s * (delta_si / d0_s) * exp(-x),  i = 1, 2
///
/// `d0_s` here is the *scaled* characteristic shear gap, i.e.
/// sqrt(2) * tangential_gap_at_maximum_shear_traction.
class SalehaniIrani3DCTraction : public TractionSeparation
{
public:
  static OptionSet expected_options();

  SalehaniIrani3DCTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Normal characteristic gap at peak normal traction.
  const Scalar & _delta0_n;

  /// Tangential characteristic gap at peak shear traction (raw, unscaled).
  const Scalar & _delta0_s;

  /// Maximum (peak) normal traction.
  const Scalar & _Tmax_n;

  /// Maximum (peak) shear traction.
  const Scalar & _Tmax_s;
};
} // namespace neml2
