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
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/clamp.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/imap.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Bilinear mixed-mode cohesive traction law with irreversible damage, optional "
                  "lagged mode mixity, and viscous regularization.";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty stiffness \\f$ K \\f$";

  options.set_parameter<TensorName<Scalar>>("critical_energy_release_rate_I");
  options.set("critical_energy_release_rate_I").doc() = "Mode-I critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("critical_energy_release_rate_II");
  options.set("critical_energy_release_rate_II").doc() = "Mode-II critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Maximum tensile traction";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Maximum tangential traction";

  options.set_parameter<TensorName<Scalar>>("eta") = "1.0";
  options.set("eta").doc() = "Mixed-mode fracture criterion exponent";

  options.set_parameter<TensorName<Scalar>>("viscosity") = "0.0";
  options.set("viscosity").doc() = "Viscous regularization coefficient";

  options.set_input("time") = VariableName(FORCES, "t");
  options.set("time").doc() = "Current time";

  EnumSelection criterion({"BK", "POWER_LAW"}, "BK");
  options.set<EnumSelection>("criterion") = criterion;
  options.set("criterion").doc() = "Mixed-mode failure criterion";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() =
      "Use the previous-step jump to compute the mode mixity and critical jumps";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use the previous-step jump to compute the effective displacement jump";

  options.set<double>("alpha") = 1e-10;
  options.set("alpha").doc() = "Half-width of the smooth Macaulay regularization interval";

  options.set_output("damage") = VariableName(STATE, "internal", "damage");
  options.set("damage").doc() = "Damage variable";

  options.set<VariableName>("mode_mixity") = VariableName(STATE, "internal", "mode_mixity");
  options.set("mode_mixity").doc() = "Optional mode mixity ratio";

  options.set<VariableName>("effective_jump_init") =
      VariableName(STATE, "internal", "effective_jump_init");
  options.set("effective_jump_init").doc() = "Optional jump at damage initiation";

  options.set<VariableName>("effective_jump_final") =
      VariableName(STATE, "internal", "effective_jump_final");
  options.set("effective_jump_final").doc() = "Optional jump at full degradation";

  options.set<VariableName>("effective_jump") = VariableName(STATE, "internal", "effective_jump");
  options.set("effective_jump").doc() = "Optional effective mixed-mode displacement jump";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "critical_energy_release_rate_I")),
    _GIIc(declare_parameter<Scalar>("GIIc", "critical_energy_release_rate_II")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "eta")),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity")),
    _time(declare_input_variable<Scalar>("time")),
    _time_old(declare_input_variable<Scalar>(_time.name().old())),
    _criterion(options.get<EnumSelection>("criterion")),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump")),
    _alpha(options.get<double>("alpha")),
    _damage(declare_output_variable<Scalar>("damage")),
    _damage_old(declare_input_variable<Scalar>(_damage.name().old())),
    _jump_old(declare_input_variable<Vec>(_jump.name().old())),
    _mode_mixity(options.get<VariableName>("mode_mixity").empty()
                     ? nullptr
                     : &declare_output_variable<Scalar>("mode_mixity")),
    _delta_init(options.get<VariableName>("effective_jump_init").empty()
                    ? nullptr
                    : &declare_output_variable<Scalar>("effective_jump_init")),
    _delta_final(options.get<VariableName>("effective_jump_final").empty()
                     ? nullptr
                     : &declare_output_variable<Scalar>("effective_jump_final")),
    _delta_m(options.get<VariableName>("effective_jump").empty()
                 ? nullptr
                 : &declare_output_variable<Scalar>("effective_jump"))
{
}

Scalar
BiLinearMixedModeTraction::regularized_heaviside(const Scalar & x, const Scalar & alpha)
{
  const auto t = clamp((x + alpha) / (2 * alpha), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

Scalar
BiLinearMixedModeTraction::regularized_heaviside_derivative(const Scalar & x, const Scalar & alpha)
{
  const auto t = clamp((x + alpha) / (2 * alpha), 0.0, 1.0);
  return 3.0 * t * (1.0 - t) / alpha;
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto jump = _jump();
  const auto zero = Scalar::zeros_like(jump(0));
  const auto one = Scalar::ones_like(jump(0));
  const auto zvec = Vec::zeros_like(jump);

  const auto mix_jump = _lag_mode_mixity ? _jump_old() : jump;
  const auto eff_jump = _lag_disp_jump ? _jump_old() : jump;

  const auto dn_mix = mix_jump(0);
  const auto ds1_mix = mix_jump(1);
  const auto ds2_mix = mix_jump(2);
  const auto ds_mix = sqrt(ds1_mix * ds1_mix + ds2_mix * ds2_mix);
  const auto opening_mix = dn_mix > zero;
  const auto tiny_shear = ds_mix <= Scalar::full_like(ds_mix, 1e-16);

  const auto beta_open = ds_mix / dn_mix;
  const auto dbeta_full = Vec::fill(
      -ds_mix / (dn_mix * dn_mix), ds1_mix / (ds_mix * dn_mix), ds2_mix / (ds_mix * dn_mix));
  const auto dbeta_axis = Vec::fill(-ds_mix / (dn_mix * dn_mix), zero, zero);
  const auto dbeta_open = where(tiny_shear, dbeta_axis, dbeta_full);
  const auto beta = where(opening_mix, beta_open, zero);
  const auto dbeta = (_lag_mode_mixity ? zvec : where(opening_mix, dbeta_open, zvec));

  const auto delta_n0 = _N / _K;
  const auto delta_s0 = _S / _K;
  const auto delta_mixed = sqrt(delta_s0 * delta_s0 + beta * beta * delta_n0 * delta_n0);
  const auto delta_init_open = delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed;
  const auto delta_init = where(opening_mix, delta_init_open, delta_s0);

  auto ddelta_init = Vec::zeros_like(dbeta);
  if (!_lag_mode_mixity)
  {
    const auto ddelta_init_dbeta =
        delta_init_open * beta *
        (1.0 / (1.0 + beta * beta) - delta_n0 * delta_n0 / (delta_mixed * delta_mixed));
    ddelta_init = where(opening_mix, ddelta_init_dbeta * dbeta, zvec);
  }

  const auto pure_shear_final = std::sqrt(2.0) * 2.0 * _GIIc / _S;
  Scalar delta_final;
  Vec ddelta_final = Vec::zeros_like(dbeta);
  if (_criterion == "BK")
  {
    const auto beta_sq_ratio = beta * beta / (1.0 + beta * beta);
    const auto delta_final_open =
        2.0 / (_K * delta_init) * (_GIc + (_GIIc - _GIc) * pow(beta_sq_ratio, _eta));
    delta_final = where(opening_mix, delta_final_open, pure_shear_final);
    if (!_lag_mode_mixity)
    {
      const auto ddelta_final_ddelta_init = -delta_final_open / delta_init;
      const auto dbeta_sq_ratio_dbeta = 2.0 * beta / pow(1.0 + beta * beta, 2.0);
      const auto ddelta_final_dbeta = 2.0 / (_K * delta_init) * (_GIIc - _GIc) * _eta *
                                      pow(beta_sq_ratio, _eta - 1.0) * dbeta_sq_ratio_dbeta;
      ddelta_final = where(
          opening_mix, ddelta_final_ddelta_init * ddelta_init + ddelta_final_dbeta * dbeta, zvec);
    }
  }
  else
  {
    const auto Gc_mixed = pow(1.0 / _GIc, _eta) + pow(beta * beta / _GIIc, _eta);
    const auto prefactor = (2.0 + 2.0 * beta * beta) / (_K * delta_init);
    const auto delta_final_open = prefactor * pow(Gc_mixed, -1.0 / _eta);
    delta_final = where(opening_mix, delta_final_open, pure_shear_final);
    if (!_lag_mode_mixity)
    {
      const auto ddelta_final_ddelta_init = -delta_final_open / delta_init;
      const auto dGc_mixed_dbeta =
          _eta * pow(beta * beta / _GIIc, _eta - 1.0) * (2.0 * beta / _GIIc);
      const auto dprefactor_dbeta = 4.0 * beta / (_K * delta_init);
      const auto dGc_term_dbeta =
          (-1.0 / _eta) * pow(Gc_mixed, -1.0 / _eta - 1.0) * dGc_mixed_dbeta;
      const auto ddelta_final_dbeta =
          dprefactor_dbeta * pow(Gc_mixed, -1.0 / _eta) + prefactor * dGc_term_dbeta;
      ddelta_final = where(
          opening_mix, ddelta_final_ddelta_init * ddelta_init + ddelta_final_dbeta * dbeta, zvec);
    }
  }

  const auto Hm = regularized_heaviside(eff_jump(0), Scalar(_alpha, eff_jump(0).options()));
  const auto dHm =
      regularized_heaviside_derivative(eff_jump(0), Scalar(_alpha, eff_jump(0).options()));
  const auto delta_n_pos_eff = Hm * eff_jump(0);
  const auto delta_m = sqrt(eff_jump(1) * eff_jump(1) + eff_jump(2) * eff_jump(2) +
                            delta_n_pos_eff * delta_n_pos_eff);
  auto ddelta_m = Vec::zeros_like(jump);
  if (!_lag_disp_jump)
  {
    const auto ddelta_n_pos = Hm + eff_jump(0) * dHm;
    ddelta_m = Vec::fill(delta_n_pos_eff * ddelta_n_pos, eff_jump(1), eff_jump(2)) / delta_m;
  }

  const auto denom = delta_m * (delta_final - delta_init);
  const auto d_trial_mid = delta_final * (delta_m - delta_init) / denom;
  const auto below = delta_m < delta_init;
  const auto above = delta_m > delta_final;
  const auto d_trial = where(below, zero, where(above, one, d_trial_mid));

  auto dd_trial = Vec::zeros_like(jump);
  const auto num1 = ddelta_final * (delta_m - delta_init) + delta_final * (ddelta_m - ddelta_init);
  const auto num2 = ddelta_m * (delta_final - delta_init) + delta_m * (ddelta_final - ddelta_init);
  const auto dd_trial_mid =
      num1 / denom - delta_final * (delta_m - delta_init) * num2 / (denom * denom);
  dd_trial = where(below, zvec, where(above, zvec, dd_trial_mid));

  auto d_no_visc = where(d_trial < _damage_old(), _damage_old(), d_trial);
  auto dd_no_visc = where(d_trial < _damage_old(), zvec, dd_trial);

  const auto dt = _time - _time_old;
  const auto visc_denom = _viscosity / dt + 1.0;
  const auto damage = (d_no_visc + _viscosity * _damage_old / dt) / visc_denom;
  const auto ddamage = dd_no_visc / visc_denom;

  const auto H = regularized_heaviside(jump(0), Scalar(_alpha, jump(0).options()));
  const auto dH = regularized_heaviside_derivative(jump(0), Scalar(_alpha, jump(0).options()));
  const auto delta_n_pos = H * jump(0);
  const auto delta_n_neg = jump(0) - delta_n_pos;
  const auto ddelta_n_pos = H + jump(0) * dH;
  const auto ddelta_n_neg = 1.0 - ddelta_n_pos;
  const auto active = Vec::fill(delta_n_pos, jump(1), jump(2));
  const auto inactive = Vec::fill(delta_n_neg, zero, zero);
  const auto Kdiag = R2::fill(_K, _K, _K);

  if (out)
  {
    _damage = damage;
    if (_mode_mixity)
      *_mode_mixity = beta;
    if (_delta_init)
      *_delta_init = delta_init;
    if (_delta_final)
      *_delta_final = delta_final;
    if (_delta_m)
      *_delta_m = delta_m;
    _traction = (1.0 - damage) * (Kdiag * active) + Kdiag * inactive;
  }

  if (dout_din && _jump.is_dependent())
  {
    _damage.d(_jump) = ddamage;
    if (_mode_mixity)
      _mode_mixity->d(_jump) = dbeta;
    if (_delta_init)
      _delta_init->d(_jump) = ddelta_init;
    if (_delta_final)
      _delta_final->d(_jump) = ddelta_final;
    if (_delta_m)
      _delta_m->d(_jump) = ddelta_m;

    const auto dactive = R2::fill(ddelta_n_pos, one, one);
    const auto dinactive = R2::fill(ddelta_n_neg, zero, zero);
    _traction.d(_jump) = (1.0 - damage) * Kdiag * dactive + Kdiag * dinactive -
                         Kdiag * neml2::outer(active, ddamage);
  }
}
} // namespace neml2
