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

/// Bilinear mixed-mode cohesive law (Camanho-Davila family).
///
/// Combines an elastic penalty branch with a damage-degraded softening
/// branch governed by a mixed-mode failure criterion (Benzeggagh-Kenane
/// or power law). The MOOSE implementation also exposes `lag_mode_mixity`
/// and `lag_displacement_jump` toggles; those, the damage state, and the
/// optional viscous regularization are NOT yet wired here. Once the
/// state-variable design is settled, declare the damage as a state input
/// (with `history_name` for `_d_old`) and produce a `residual_name(d)`
/// output for the integrator.
class BiLinearMixedModeTraction : public TractionSeparation
{
public:
  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Penalty stiffness K.
  const Scalar & _K;

  /// Mode I critical energy release rate.
  const Scalar & _GIc;

  /// Mode II critical energy release rate.
  const Scalar & _GIIc;

  /// Tensile (normal) strength N.
  const Scalar & _N;

  /// Shear strength S.
  const Scalar & _S;

  /// Mixed-mode criterion exponent (BK or power law).
  const Scalar & _eta;

  /// Damage variable d in [0, 1] (output, treated reversibly here for now).
  Variable<Scalar> & _d;
};
} // namespace neml2
