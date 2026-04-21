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
#include "neml2/tensors/Scalar.h"

namespace neml2
{
/**
 * @brief Bilinear mixed-mode cohesive traction law with scalar damage.
 */
class BiLinearMixedModeTraction : public TractionSeparation
{
public:
  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void request_AD() override;
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  const Variable<Vec> & _jump_old;
  const Variable<Scalar> & _d_old;

  Variable<Scalar> & _d;
  Variable<Scalar> & _beta_out;
  Variable<Scalar> & _delta_init_out;
  Variable<Scalar> & _delta_final_out;
  Variable<Scalar> & _delta_m_out;

  const Scalar & _K;
  const Scalar & _G_Ic;
  const Scalar & _G_IIc;
  const Scalar & _N;
  const Scalar & _S;
  const Scalar & _eta;
  const Scalar & _viscosity;
  const Scalar & _dt;

  const std::string _criterion;
  const bool _lag_mode_mixity;
  const bool _lag_displacement_jump;
  const double _alpha;
};
} // namespace neml2
