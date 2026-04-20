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

#include "neml2/models/solid_mechanics/cohesive/TractionSeparationLaw.h"

namespace neml2
{
class Scalar;

/**
 * @brief Scalar-damage exponential traction-separation law.
 *
 * Damage: d = 1 - exp(-delta_eff / delta0)
 * Traction: T = (1 - d) * (Gc / delta0^2) * delta
 *
 * Supports optional irreversible damage via a max-effective-jump state variable.
 */
class ExpTractionSeparation : public TractionSeparationLaw
{
public:
  static OptionSet expected_options();

  ExpTractionSeparation(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Fracture energy G_c
  const Scalar & _Gc;

  /// Softening length scale delta_0
  const Scalar & _delta0;

  /// Tangential weighting factor beta
  const Scalar & _beta;

  /// Whether to enforce irreversible damage (no healing)
  const bool _irreversible;

  /// Historical maximum effective scalar jump (state output, updated when delta_eff advances)
  Variable<Scalar> & _delta_eff_max;

  /// Historical maximum effective scalar jump from previous step (state input)
  const Variable<Scalar> & _delta_eff_max_old;
};
} // namespace neml2
