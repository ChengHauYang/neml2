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

#include "neml2/base/EnumSelection.h"
#include "neml2/models/solid_mechanics/traction_separation/TractionSeparationLaw.h"

namespace neml2
{
class Scalar;

class BiLinearMixedModeTraction : public TractionSeparationLaw
{
public:
  enum class Criterion : int
  {
    BK = 0,
    POWER_LAW = 1
  };

  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  static Scalar regularized_heaviside(const Scalar & x, const Scalar & alpha);
  static Scalar regularized_heaviside_derivative(const Scalar & x, const Scalar & alpha);

  const Scalar & _K;
  const Scalar & _GIc;
  const Scalar & _GIIc;
  const Scalar & _N;
  const Scalar & _S;
  const Scalar & _eta;
  const Scalar & _viscosity;
  const Variable<Scalar> & _time;
  const Variable<Scalar> & _time_old;

  const EnumSelection _criterion;
  const bool _lag_mode_mixity;
  const bool _lag_disp_jump;
  const double _alpha;

  Variable<Scalar> & _damage;
  const Variable<Scalar> & _damage_old;
  const Variable<Vec> & _jump_old;
  Variable<Scalar> * _mode_mixity;
  Variable<Scalar> * _delta_init;
  Variable<Scalar> * _delta_final;
  Variable<Scalar> * _delta_m;
};
} // namespace neml2
