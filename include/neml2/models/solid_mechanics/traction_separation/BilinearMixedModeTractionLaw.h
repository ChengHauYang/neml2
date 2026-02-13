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

#include "CohesiveTractionLaw.h"

namespace neml2
{

/// Mixed-mode bilinear traction-separation law (damage-based)
class BilinearMixedModeTractionLaw : public CohesiveTractionLaw
{
public:
  static OptionSet expected_options();
  BilinearMixedModeTractionLaw(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// penalty stiffness
  const Scalar & _K;
  /// critical energy release rate in normal direction
  const Scalar & _GIc;
  /// critical energy release rate in shear direction
  const Scalar & _GIIc;
  /// normal strength
  const Scalar & _N;
  /// shear strength
  const Scalar & _S;
  /// power law parameter
  const Scalar & _eta;
  /// viscosity (0 -> off)
  const Scalar & _visc;
  /// smoothing_length
  const Scalar & _alpha;

  /// Criterion: "BK" or "POWER_LAW"
  const std::string _criterion;

  /// Whether to use old displacement jumps to compute the mode mixity
  const bool _lag_mode_mixity;
  /// Whether to use old displacement jumps to compute the effective displacement jump
  const bool _lag_disp_jump;

  /// time step size
  const Variable<Scalar> & _dt;
  /// old jump (for lagging)
  const Variable<Vec> & _jump_old;
  /// old damage
  const Variable<Scalar> & _damage_old;
  /// Output: damage variable (0: undamaged, 1: fully damaged)
  Variable<Scalar> & _damage;

  // Optional helpful outputs
  Variable<Scalar> & _beta_out;
  Variable<Scalar> & _delta_init_out;
  Variable<Scalar> & _delta_final_out;
  Variable<Scalar> & _delta_m_out;

  // -------------------------
  // Helpers (MOOSE-aligned)
  // -------------------------
  Scalar regularizedHeaviside(const Scalar & x, const Scalar & smoothing_length) const;
  Scalar regularizedHeavisideDerivative(const Scalar & x, const Scalar & smoothing_length) const;

  Scalar macauley_pos(const Scalar & x, const Scalar & smoothing_length) const;
};

} // namespace neml2
