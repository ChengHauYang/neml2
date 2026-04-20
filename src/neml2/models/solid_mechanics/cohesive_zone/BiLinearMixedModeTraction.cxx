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

#include "neml2/models/solid_mechanics/cohesive_zone/BiLinearMixedModeTraction.h"
#include "neml2/tensors/tensors.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/outer.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Bilinear mixed-mode traction separation law.";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty elastic stiffness K";

  options.set_parameter<TensorName<Scalar>>("mode_I_critical_energy");
  options.set("mode_I_critical_energy").doc() = "Mode I critical energy release rate GIc";

  options.set_parameter<TensorName<Scalar>>("mode_II_critical_energy");
  options.set("mode_II_critical_energy").doc() = "Mode II critical energy release rate GIIc";

  options.set_parameter<TensorName<Scalar>>("tensile_strength");
  options.set("tensile_strength").doc() = "Tensile strength N";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength S";

  options.set_parameter<TensorName<Scalar>>("power_law_exponent");
  options.set("power_law_exponent").doc() = "Power law exponent eta";
  options.set<TensorName<Scalar>>("power_law_exponent") = "1.0";

  options.set<std::string>("criterion") = "BK";
  options.set("criterion").doc() = "Mixed mode propagation criterion (BK or POWER_LAW)";

  options.set<double>("macaulay_regularization") = 1e-10;
  options.set("macaulay_regularization").doc() = "Regularization parameter for Macaulay bracket";

  options.set_output("damage") = VariableName(STATE, "internal", "d");
  options.set("damage").doc() = "Damage variable";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness", false)),
    _GIc(declare_parameter<Scalar>("GIc", "mode_I_critical_energy", false)),
    _GIIc(declare_parameter<Scalar>("GIIc", "mode_II_critical_energy", false)),
    _N(declare_parameter<Scalar>("N", "tensile_strength", false)),
    _S(declare_parameter<Scalar>("S", "shear_strength", false)),
    _eta(declare_parameter<Scalar>("eta", "power_law_exponent", false)),
    _criterion(options.get<std::string>("criterion")),
    _alpha(options.get<double>("macaulay_regularization")),
    _d(declare_output_variable<Scalar>("damage")),
    _d_old(declare_input_variable<Scalar>(_d.name().old()))
{
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, [[maybe_unused]] bool d2out_din2)
{
  auto delta = _delta();
  auto delta_n = delta(0);
  auto delta_s = sqrt(delta(1) * delta(1) + delta(2) * delta(2));

  auto zero = Scalar(0.0, delta.options());
  auto one = Scalar(1.0, delta.options());
  auto eps = Scalar(1e-16, delta.options());

  // 1. Mode mixity
  auto beta = neml2::where(delta_n > zero, delta_s / (delta_n + eps), zero);

  // 2. Critical displacement jump
  auto delta_n0 = _N / _K;
  auto delta_s0 = _S / _K;
  auto delta_mixed = sqrt(delta_s0 * delta_s0 + (beta * delta_n0) * (beta * delta_n0));
  auto delta_init_mixed = (delta_n0 * delta_s0 * sqrt(one + beta * beta)) / delta_mixed;
  auto delta_init = neml2::where(delta_n > zero, delta_init_mixed, delta_s0);

  // 3. Final displacement jump
  auto delta_final_pure_shear = 2.0 * _GIIc / _S;
  auto delta_final = delta_final_pure_shear;

  if (_criterion == "BK")
  {
    auto beta_sq_ratio = (beta * beta) / (one + beta * beta);
    auto delta_final_mixed =
        2.0 / (_K * delta_init) * (_GIc + (_GIIc - _GIc) * pow(beta_sq_ratio, _eta));
    delta_final = neml2::where(delta_n > zero, delta_final_mixed, delta_final_pure_shear);
  }
  else if (_criterion == "POWER_LAW")
  {
    auto Gc_mixed = pow(one / _GIc, _eta) + pow(beta * beta / _GIIc, _eta);
    auto delta_final_mixed = (2.0 + 2.0 * beta * beta) / (_K * delta_init) * pow(Gc_mixed, -one / _eta);
    delta_final = neml2::where(delta_n > zero, delta_final_mixed, delta_final_pure_shear);
  }

  // 4. Effective displacement jump
  // Regularized Macaulay bracket for normal component
  auto H = 0.5 * (one + delta_n / sqrt(delta_n * delta_n + Scalar(_alpha * _alpha, delta.options())));
  auto delta_n_pos = H * delta_n;
  auto delta_m = sqrt(delta_n_pos * delta_n_pos + delta_s * delta_s);

  // 5. Damage
  auto d_trial = delta_final * (delta_m - delta_init) / (delta_m * (delta_final - delta_init) + eps);
  d_trial = neml2::where(delta_m < delta_init, zero, d_trial);
  d_trial = neml2::where(delta_m > delta_final, one, d_trial);
  d_trial = neml2::where(d_trial < zero, zero, d_trial);

  auto d_old = _d_old();
  auto d = neml2::where(d_trial > d_old, d_trial, d_old);

  // 6. Traction
  auto delta_n_neg = delta_n - delta_n_pos;
  auto delta_active = Vec::fill(delta_n_pos, delta(1), delta(2));
  auto delta_inactive = Vec::fill(delta_n_neg, zero, zero);

  if (out)
  {
    _d = Tensor(d);
    _t = Tensor((one - d) * _K * delta_active + _K * delta_inactive);
  }

  if (dout_din)
  {
    if (_delta.is_dependent())
    {
      _t.d(_delta) = Tensor((one - d) * _K * R2::identity(_delta.options()));
    }
    if (_d_old.is_dependent())
    {
       _t.d(_d_old) = Tensor(-_K * delta_active * neml2::where(d_trial > d_old, zero, one));
       _d.d(_d_old) = Tensor(neml2::where(d_trial > d_old, zero, one));
    }
  }
}
} // namespace neml2
