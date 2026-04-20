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

#include "neml2/models/solid_mechanics/traction_separation/TractionSeparation.h"

namespace neml2
{
class Scalar;

/**
 * @brief Exponential traction-separation law.
 *
 * Implements an exponential softening law: T = (1 - d) * K * delta, where d = 1 - exp(-delta_eff / delta0).
 */
class ExpTractionSeparation : public TractionSeparation
{
public:
  static OptionSet expected_options();

  ExpTractionSeparation(const OptionSet & options);

  void request_AD() override;

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Fracture energy
  const Scalar & _Gc;

  /// Softening length scale
  const Scalar & _delta0;

  /// Tangential weighting factor
  const Scalar & _beta;

  /// Small regularizer for norm
  const Scalar & _eps;

  /// Irreversible damage
  bool _irreversible_damage;

  /// Maximum effective scalar jump (state)
  Variable<Scalar> & _delta_eff_max;

  /// Old maximum effective scalar jump (state input)
  const Variable<Scalar> & _delta_eff_max_old;

  /// Damage variable (derived output)
  Variable<Scalar> & _damage;
};
} // namespace neml2
