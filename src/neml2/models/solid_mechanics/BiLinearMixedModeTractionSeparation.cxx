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

#include "neml2/models/solid_mechanics/BiLinearMixedModeTractionSeparation.h"
#include "neml2/base/EnumSelection.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/abs.h"
#include "neml2/tensors/functions/cosh.h"
#include "neml2/tensors/functions/norm.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/tanh.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTractionSeparation);

OptionSet
BiLinearMixedModeTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Bilinear mixed-mode traction-separation law with optional lagged mixity and "
                  "viscous damage regularization.";

  EnumSelection criterion_selection({"BK", "POWER_LAW"}, "BK");
  options.set<EnumSelection>("criterion") = criterion_selection;
  options.set("criterion").doc() =
      "Mixed-mode propagation criterion: " + criterion_selection.join();

  options.set<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty elastic stiffness";

  options.set<TensorName<Scalar>>("mode_I_fracture_energy");
  options.set("mode_I_fracture_energy").doc() = "Mode I critical energy release rate";

  options.set<TensorName<Scalar>>("mode_II_fracture_energy");
  options.set("mode_II_fracture_energy").doc() = "Mode II critical energy release rate";

  options.set<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Maximum normal traction";

  options.set<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Maximum shear traction";

  options.set<TensorName<Scalar>>("power_law_exponent");
  options.set("power_law_exponent").doc() = "Exponent used by the mixed-mode criterion";

  options.set<TensorName<Scalar>>("viscosity") = TensorName<Scalar>("0");
  options.set("viscosity").doc() = "Damage viscosity coefficient";

  options.set_input("old_displacement_jump") =
      VariableName(OLD_FORCES, "interface_displacement_jump");
  options.set("old_displacement_jump").doc() = "Previous displacement jump";

  options.set_input("old_damage") = VariableName(OLD_STATE, "internal", "damage");
  options.set("old_damage").doc() = "Previous cohesive damage";

  options.set_input("time") = VariableName(FORCES, "t");
  options.set("time").doc() = "Current time";

  options.set_input("old_time") = VariableName(OLD_FORCES, "t");
  options.set("old_time").doc() = "Previous time";

  options.set_output("damage") = VariableName(STATE, "internal", "damage");
  options.set("damage").doc() = "Updated cohesive damage";

  options.set_output("mode_mixity");
  options.set("mode_mixity").doc() = "Optional mixed-mode ratio";

  options.set_output("critical_displacement_jump");
  options.set("critical_displacement_jump").doc() =
      "Optional effective displacement jump at damage initiation";

  options.set_output("final_displacement_jump");
  options.set("final_displacement_jump").doc() =
      "Optional effective displacement jump at full degradation";

  options.set_output("effective_displacement_jump");
  options.set("effective_displacement_jump").doc() =
      "Optional current effective displacement jump magnitude";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() =
      "Use the previous displacement jump when computing the mixed-mode ratio";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use the previous displacement jump when computing the effective displacement jump";

  options.set<double>("alpha") = 1e-10;
  options.set("alpha").doc() = "Regularization scale for the normal opening split";

  options.set<double>("eps") = 1e-16;
  options.set("eps").doc() = "Small regularizer for divisions by near-zero quantities";

  return options;
}

BiLinearMixedModeTractionSeparation::BiLinearMixedModeTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "mode_I_fracture_energy")),
    _GIIc(declare_parameter<Scalar>("GIIc", "mode_II_fracture_energy")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "power_law_exponent")),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity")),
    _jump_old(declare_input_variable<Vec>("old_displacement_jump")),
    _damage_old(declare_input_variable<Scalar>("old_damage")),
    _time(declare_input_variable<Scalar>("time")),
    _time_old(declare_input_variable<Scalar>("old_time")),
    _damage(declare_output_variable<Scalar>("damage")),
    _mode_mixity(options.get("mode_mixity").user_specified()
                     ? &declare_output_variable<Scalar>("mode_mixity")
                     : nullptr),
    _delta_init(options.get("critical_displacement_jump").user_specified()
                    ? &declare_output_variable<Scalar>("critical_displacement_jump")
                    : nullptr),
    _delta_final(options.get("final_displacement_jump").user_specified()
                     ? &declare_output_variable<Scalar>("final_displacement_jump")
                     : nullptr),
    _delta_m(options.get("effective_displacement_jump").user_specified()
                 ? &declare_output_variable<Scalar>("effective_displacement_jump")
                 : nullptr),
    _criterion(options.get<EnumSelection>("criterion").as<MixedModeCriterion>()),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_jump(options.get<bool>("lag_displacement_jump")),
    _alpha(options.get<double>("alpha")),
    _eps(options.get<double>("eps"))
{
}

void
BiLinearMixedModeTractionSeparation::request_AD()
{
  TractionSeparationLaw::request_AD();
  _damage.request_AD(_jump);
  if (_mode_mixity)
    _mode_mixity->request_AD(_jump);
  if (_delta_init)
    _delta_init->request_AD(_jump);
  if (_delta_final)
    _delta_final->request_AD(_jump);
  if (_delta_m)
    _delta_m->request_AD(_jump);
}

Scalar
BiLinearMixedModeTractionSeparation::regularized_heaviside(const Scalar & x, double alpha)
{
  return 0.5 * (tanh(x / alpha) + 1.0);
}

void
BiLinearMixedModeTractionSeparation::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto jump = _jump();
  const auto jump_mix = _lag_mode_mixity ? _jump_old() : jump;
  const auto jump_eff = _lag_jump ? _jump_old() : jump;

  const Scalar zero = Scalar::zeros_like(jump(0));
  const Scalar one = Scalar::ones_like(jump(0));

  const Scalar delta_n_mix = jump_mix(0);
  const Scalar delta_s_mix = sqrt(jump_mix(1) * jump_mix(1) + jump_mix(2) * jump_mix(2));
  const Scalar opening = delta_n_mix > 0.0;
  const Scalar safe_delta_n_mix = where(opening, delta_n_mix, one);
  const Scalar beta = where(opening, delta_s_mix / safe_delta_n_mix, zero);

  const Scalar delta_n0 = _N / _K;
  const Scalar delta_s0 = _S / _K;
  const Scalar delta_mixed0 = sqrt(delta_s0 * delta_s0 + beta * beta * delta_n0 * delta_n0 + _eps);
  const Scalar delta_init =
      where(opening, delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed0, delta_s0);

  const Scalar pure_shear_delta_final = std::sqrt(2.0) * 2.0 * _GIIc / _S;
  Scalar delta_final(pure_shear_delta_final);
  if (_criterion == MixedModeCriterion::BK)
  {
    const Scalar beta_sq_ratio = beta * beta / (1.0 + beta * beta);
    delta_final =
        where(opening,
              2.0 / (_K * delta_init) * (_GIc + (_GIIc - _GIc) * pow(beta_sq_ratio, _eta)),
              pure_shear_delta_final);
  }
  else
  {
    const Scalar Gc_mixed = pow(1.0 / _GIc, _eta) + pow(beta * beta / _GIIc, _eta);
    delta_final = where(opening,
                        (2.0 + 2.0 * beta * beta) / (_K * delta_init) * pow(Gc_mixed, -1.0 / _eta),
                        delta_final);
  }

  const Scalar H_mix = regularized_heaviside(jump_eff(0), _alpha);
  const Scalar delta_n_pos = H_mix * jump_eff(0);
  const Scalar delta_m = sqrt(delta_n_pos * delta_n_pos + jump_eff(1) * jump_eff(1) +
                              jump_eff(2) * jump_eff(2) + _eps);

  const Scalar below = delta_m < delta_init;
  const Scalar above = delta_m > delta_final;
  const Scalar denom = delta_m * (delta_final - delta_init);
  const Scalar safe_denom = where(abs(denom) > _eps, denom, one);
  const Scalar d_soft = delta_final * (delta_m - delta_init) / safe_denom;
  const Scalar d_trial = where(below, zero, where(above, one, d_soft));
  const Scalar d_hist = where(d_trial < _damage_old(), _damage_old(), d_trial);

  const Scalar dt = _time() - _time_old();
  const Scalar safe_dt = where(abs(dt) > _eps, dt, one);
  const Scalar visc_factor = one + _viscosity / safe_dt;
  const Scalar damage = (d_hist + _viscosity * _damage_old() / safe_dt) / visc_factor;

  const Scalar H = regularized_heaviside(jump(0), _alpha);
  const Scalar delta_n_active = H * jump(0);
  const Scalar delta_n_neg = jump(0) - delta_n_active;
  const Vec delta_active = Vec::fill(delta_n_active, jump(1), jump(2));
  const Vec delta_inactive = Vec::fill(delta_n_neg, zero, zero);

  _damage = damage;
  _traction = (1.0 - damage) * _K * delta_active + _K * delta_inactive;

  if (_mode_mixity)
    *_mode_mixity = beta;
  if (_delta_init)
    *_delta_init = delta_init;
  if (_delta_final)
    *_delta_final = delta_final;
  if (_delta_m)
    *_delta_m = delta_m;
}
} // namespace neml2
