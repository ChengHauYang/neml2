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

#include "neml2/models/solid_mechanics/traction_separation/BiLinearMixedModeTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/macaulay.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() +=
      " Camanho-Davila bilinear envelope with Benzeggagh-Kenane (BK) or power-law mixed-mode "
      "propagation criterion. Damage is irreversible; normal compression is recovered through a "
      "Macaulay split.";

  options.add_parameter<Scalar>("penalty_stiffness", "Penalty elastic stiffness K");
  options.add_parameter<Scalar>("mode_I_fracture_energy", "Mode I critical energy release rate");
  options.add_parameter<Scalar>("mode_II_fracture_energy", "Mode II critical energy release rate");
  options.add_parameter<Scalar>("normal_strength", "Tensile (normal) strength N");
  options.add_parameter<Scalar>("shear_strength", "Shear strength S");
  options.add_parameter<Scalar>("power_law_exponent",
                                "Mixed-mode criterion exponent (BK or power-law)");

  EnumSelection criterion({"BK", "POWER_LAW"}, "BK");
  options.add<EnumSelection>(
      "criterion",
      criterion,
      "Mixed-mode propagation criterion. Options are: " + criterion.join());
  options.add<double>("epsilon",
                      1e-16,
                      "Small regularizer added inside sqrt() to keep its derivative bounded at "
                      "zero displacement jump");

  options.add_output("damage", "Scalar damage variable (irreversible)");

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _d(declare_output_variable<Scalar>("damage")),
    _d_old(declare_input_variable<Scalar>(history_name(_d.name(), /*nstep=*/1))),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "mode_I_fracture_energy")),
    _GIIc(declare_parameter<Scalar>("GIIc", "mode_II_fracture_energy")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "power_law_exponent")),
    _criterion(options.get<EnumSelection>("criterion")),
    _eps(options.get<double>("epsilon"))
{
}

void
BiLinearMixedModeTraction::request_AD()
{
  // Hand-rolling the analytical Jacobian through the mode-mixity branch, the
  // criterion branch, the damage-regime branch, the Macaulay split, and the
  // irreversibility max-clamp is too brittle to be worth the maintenance cost.
  Model::request_AD(_T, _delta);
  Model::request_AD(_T, _d_old);
  Model::request_AD(_d, _delta);
  Model::request_AD(_d, _d_old);
}

void
BiLinearMixedModeTraction::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto dn = _delta()(0);
  const auto ds1 = _delta()(1);
  const auto ds2 = _delta()(2);

  const auto zero_s = Scalar::zeros_like(dn);
  const auto one_s = Scalar::ones_like(dn);

  // Macaulay split on the normal jump
  const auto delta_n_pos = neml2::macaulay(dn);
  const auto delta_n_neg = dn - delta_n_pos;

  // Tangential magnitude (regularized so its derivative is finite at zero)
  const auto delta_s_sq = ds1 * ds1 + ds2 * ds2;
  const auto delta_s = neml2::sqrt(delta_s_sq + _eps);

  // Mode mixity ratio beta = delta_s / delta_n is only defined for opening (delta_n > 0).
  // Use a safe denominator in the false branch so torch's autograd does not see NaN
  // gradients in the unused path.
  const auto pos_mask = dn > 0.0;
  const auto safe_delta_n_pos = neml2::where(pos_mask, delta_n_pos, one_s);
  const auto beta = neml2::where(pos_mask, delta_s / safe_delta_n_pos, zero_s);
  const auto beta_sq = beta * beta;

  // Initiation displacement jump (mixed-mode, falls back to pure-shear when no opening)
  const auto delta_normal0 = _N / _K;
  const auto delta_shear0 = _S / _K;
  const auto delta_mixed_init =
      neml2::sqrt(delta_shear0 * delta_shear0 + beta_sq * delta_normal0 * delta_normal0 + _eps);
  const auto delta_init_mixed =
      delta_normal0 * delta_shear0 * neml2::sqrt(1.0 + beta_sq) / delta_mixed_init;
  const auto delta_init = neml2::where(pos_mask, delta_init_mixed, delta_shear0);

  // Final (full-degradation) displacement jump
  Scalar delta_final_mixed = zero_s;
  if (_criterion == "BK")
  {
    const auto beta_sq_ratio = beta_sq / (1.0 + beta_sq);
    // Adding a tiny constant to the base avoids 0^eta NaN gradients in the unused branch
    const auto term = _GIc + (_GIIc - _GIc) * neml2::pow(beta_sq_ratio + _eps, _eta);
    delta_final_mixed = 2.0 / (_K * delta_init) * term;
  }
  else // POWER_LAW
  {
    const auto Gc_mixed =
        neml2::pow(1.0 / _GIc, _eta) + neml2::pow(beta_sq / _GIIc + _eps, _eta);
    delta_final_mixed = (2.0 + 2.0 * beta_sq) / (_K * delta_init) * neml2::pow(Gc_mixed, -1.0 / _eta);
  }
  const auto delta_final_default =
      Scalar::full_like(dn, std::sqrt(2.0)) * 2.0 * _GIIc / _S;
  const auto delta_final = neml2::where(pos_mask, delta_final_mixed, delta_final_default);

  // Effective mixed-mode displacement jump
  const auto delta_m =
      neml2::sqrt(delta_n_pos * delta_n_pos + delta_s_sq + _eps);

  // Bilinear damage (with safe denominator so the unused branch gradient is finite)
  const auto safe_denom = neml2::where(delta_final > delta_init, delta_final - delta_init, one_s);
  const auto bilinear_d = delta_final * (delta_m - delta_init) / (delta_m * safe_denom);
  const auto d_trial =
      neml2::where(delta_m < delta_init,
                   zero_s,
                   neml2::where(delta_m < delta_final, bilinear_d, one_s));

  // Irreversibility
  const auto d = neml2::where(d_trial > _d_old(), d_trial, _d_old());

  // Traction
  const auto T_active = _K * (1.0 - d) * Vec::fill(delta_n_pos, ds1, ds2);
  const auto T_inactive = _K * Vec::fill(delta_n_neg, zero_s, zero_s);

  _T = T_active + T_inactive;
  _d = d;
}
} // namespace neml2
