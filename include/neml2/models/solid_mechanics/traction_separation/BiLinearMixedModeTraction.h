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

namespace neml2
{
/// Mixed-mode bilinear traction-separation law (Camanho-Davila) with BK or POWER_LAW criterion
class BiLinearMixedModeTraction : public TractionSeparationLaw
{
public:
  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

private:
  // ---- Parameters ----
  const Scalar & _K;      ///< Penalty stiffness
  const Scalar & _GIc;    ///< Mode-I critical energy release rate
  const Scalar & _GIIc;   ///< Mode-II critical energy release rate
  const Scalar & _N;      ///< Normal (tensile) strength
  const Scalar & _S;      ///< Shear strength
  const Scalar & _eta;    ///< Power-law exponent (BK or POWER_LAW)
  const Scalar & _visc;   ///< Viscous regularization coefficient (0 = off)
  const Scalar & _alpha;  ///< Smoothing length for regularized Heaviside

  // ---- Control ----
  const std::string _criterion;       ///< "BK" or "POWER_LAW"
  const bool _lag_mode_mixity;        ///< Use old jump for beta computation
  const bool _lag_disp_jump;          ///< Use old jump for effective jump computation

  // ---- State inputs ----
  const Variable<Scalar> & _dt;       ///< Time step size
  const Variable<Vec> * _jump_old;    ///< Old displacement jump (lagging)
  const Variable<Scalar> & _damage_old; ///< Old damage for irreversibility

  // ---- State output ----
  Variable<Scalar> & _damage;         ///< Damage variable

  // ---- Helpers ----
  Scalar reg_heaviside(const Scalar & x) const;
  Scalar reg_heaviside_deriv(const Scalar & x) const;
};
} // namespace neml2
