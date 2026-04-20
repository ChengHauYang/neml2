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
#include <cstdint>

namespace neml2
{
class Scalar;
class Vec;

/**
 * @brief Bilinear mixed-mode cohesive traction-separation law.
 *
 * Implements the bilinear damage traction-separation model with:
 *   - Mode-mixity via the Benzeggagh-Kenane (BK) or POWER_LAW criterion
 *   - Optional irreversible damage (always enforced via max-damage state)
 *   - Optional viscous regularization: d = (d_trial + visc*d_old/dt) / (visc/dt + 1)
 *   - Optional lagging of mode-mixity and effective-jump quantities
 */
class BiLinearMixedModeTraction : public TractionSeparationLaw
{
public:
  enum class Criterion : std::uint8_t
  {
    BK,
    POWER_LAW
  };

  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Penalty elastic stiffness K
  const Scalar & _K;

  /// Mode I critical energy release rate G_Ic
  const Scalar & _GI_c;

  /// Mode II critical energy release rate G_IIc
  const Scalar & _GII_c;

  /// Tensile (normal) strength N
  const Scalar & _N;

  /// Shear strength S
  const Scalar & _S;

  /// Power law exponent eta
  const Scalar & _eta;

  /// Viscous regularization coefficient
  const Scalar & _viscosity;

  /// Mixed-mode propagation criterion
  const Criterion _criterion;

  /// Use previous-step jump when computing mode-mixity, delta_init, delta_final
  const bool _lag_mode_mixity;

  /// Use previous-step jump when computing effective mixed-mode jump delta_m
  const bool _lag_displacement_jump;

  /// Time step size (FORCES input)
  const Variable<Scalar> & _dt;

  /// Scalar damage variable d in [0, 1] (state output)
  Variable<Scalar> & _damage;

  /// Damage from previous step (old state input)
  const Variable<Scalar> & _damage_old;

  /// Stored displacement jump for lag computations (state output)
  Variable<Vec> & _jump_stored;

  /// Stored displacement jump from previous step (old state input)
  const Variable<Vec> & _jump_old;
};
} // namespace neml2
