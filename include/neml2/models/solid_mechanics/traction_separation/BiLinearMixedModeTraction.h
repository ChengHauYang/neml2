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
#include "neml2/base/EnumSelection.h"

namespace neml2
{
class Scalar;

/**
 * @brief Bilinear mixed-mode traction-separation law.
 *
 * Implements mixed-mode cohesive behavior with either BK or POWER_LAW criterion.
 */
class BiLinearMixedModeTraction : public TractionSeparation
{
public:
  enum class Criterion
  {
    BK,
    POWER_LAW
  };

  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

  void request_AD() override;

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Mode-mixity ratio criterion
  Criterion _criterion;

  /// Normal stiffness
  const Scalar & _K;

  /// Mode I critical energy release rate
  const Scalar & _GI_c;

  /// Mode II critical energy release rate
  const Scalar & _GII_c;

  /// Tensile (normal) strength
  const Scalar & _N;

  /// Shear strength
  const Scalar & _S;

  /// Power law exponent
  const Scalar & _eta;

  /// Viscous regularization coefficient
  const Scalar & _viscosity;

  /// Regularization parameter for Macaulay bracket
  const Scalar & _alpha;

  /// Damage variable (state)
  Variable<Scalar> & _d;

  /// Old damage variable (state input)
  const Variable<Scalar> & _d_old;

  /// Old displacement jump (state input)
  const Variable<Vec> & _displacement_jump_old;

  /// Time
  const Variable<Scalar> & _t;

  /// Old time
  const Variable<Scalar> & _t_old;

  /// Lag mode mixity
  bool _lag_mode_mixity;

  /// Lag displacement jump
  bool _lag_disp_jump;
};
} // namespace neml2
