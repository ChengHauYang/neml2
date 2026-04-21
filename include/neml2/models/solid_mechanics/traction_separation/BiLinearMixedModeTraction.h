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

#include "neml2/models/solid_mechanics/traction_separation/TractionSeparationLaw.h"
#include "neml2/tensors/Scalar.h"

namespace neml2
{
/**
 * @brief Bilinear mixed-mode cohesive zone law with Benzeggagh-Kenane or power-law criterion.
 *
 * Implements mode-mixity-dependent damage initiation and propagation with optional
 * irreversibility, viscous regularization, and lagging of mode-mixity / displacement jump.
 */
class BiLinearMixedModeTraction : public TractionSeparationLaw
{
public:
  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Penalty stiffness K
  const Scalar & _K;

  /// Mode I critical energy release rate
  const Scalar & _GI_c;

  /// Mode II critical energy release rate
  const Scalar & _GII_c;

  /// Normal (tensile) strength
  const Scalar & _N;

  /// Shear strength
  const Scalar & _S;

  /// BK/power-law exponent
  const Scalar & _eta;

  /// Viscous regularization coefficient (0 = no regularization)
  const double _viscosity;

  /// Regularization parameter for the smooth Heaviside (tanh-based)
  const double _alpha;

  /// Whether to use the power-law criterion instead of BK
  const bool _use_power_law;

  /// Whether to lag the mode mixity computation using the old displacement jump
  const bool _lag_mode_mixity;

  /// Whether to lag the effective displacement jump using the old displacement jump
  const bool _lag_disp_jump;

  /// Current damage (state output)
  Variable<Scalar> & _damage;

  /// Old damage (old state input for irreversibility)
  const Variable<Scalar> & _damage_old;

  /// Old displacement jump (old forces, used when lagging is active)
  const Variable<Vec> & _displacement_jump_old;

  /// Current time (for viscous regularization)
  const Variable<Scalar> & _t;

  /// Old time (for viscous regularization)
  const Variable<Scalar> & _t_old;
};
} // namespace neml2
