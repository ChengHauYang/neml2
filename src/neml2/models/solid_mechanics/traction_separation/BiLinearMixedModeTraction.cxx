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
#include "neml2/misc/assertions.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/tanh.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

namespace
{
Scalar
regularized_heaviside(const Scalar & x, double alpha)
{
  return 0.5 * (1.0 + tanh(x / alpha));
}

Scalar
positive_part(const Scalar & x, double alpha)
{
  return regularized_heaviside(x, alpha) * x;
}
} // namespace

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Bilinear mixed-mode cohesive traction law with irreversible damage.";

  options.set_output("damage") = VariableName(STATE, "internal", "interface_damage");
  options.set("damage").doc() = "Scalar interface damage";

  options.set_output("mode_mixity") = VariableName(STATE, "internal", "mode_mixity");
  options.set("mode_mixity").doc() = "Mode mixity ratio";

  options.set_output("damage_initiation_jump") =
      VariableName(STATE, "internal", "damage_initiation_jump");
  options.set("damage_initiation_jump").doc() = "Effective displacement jump at damage initiation";

  options.set_output("full_degradation_jump") =
      VariableName(STATE, "internal", "full_degradation_jump");
  options.set("full_degradation_jump").doc() = "Effective displacement jump at full degradation";

  options.set_output("effective_displacement_jump") =
      VariableName(STATE, "internal", "effective_displacement_jump");
  options.set("effective_displacement_jump").doc() = "Current effective mixed-mode jump";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty elastic stiffness";

  options.set_parameter<TensorName<Scalar>>("mode_I_critical_energy_release_rate");
  options.set("mode_I_critical_energy_release_rate").doc() = "Mode I critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("mode_II_critical_energy_release_rate");
  options.set("mode_II_critical_energy_release_rate").doc() =
      "Mode II critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Tensile normal strength";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength";

  options.set_parameter<TensorName<Scalar>>("power_law_exponent") = "1.0";
  options.set("power_law_exponent").doc() = "Mixed-mode criterion exponent";

  options.set_parameter<TensorName<Scalar>>("viscosity") = "0.0";
  options.set("viscosity").doc() = "Viscous damage regularization coefficient";

  options.set_parameter<TensorName<Scalar>>("time_step") = "1.0";
  options.set("time_step").doc() = "Time step used by viscous damage regularization";

  options.set<std::string>("criterion") = "BK";
  options.set("criterion").doc() = "Mixed-mode propagation criterion, either BK or POWER_LAW";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() = "Use the old jump to compute mode mixity";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use the old jump to compute the effective mixed-mode jump";

  options.set<double>("regularization") = 1e-10;
  options.set("regularization").doc() = "Regularization parameter for positive normal opening";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _jump_old(declare_input_variable<Vec>(_jump.name().old())),
    _d_old(declare_input_variable<Scalar>(options.get<VariableName>("damage").old())),
    _d(declare_output_variable<Scalar>("damage")),
    _beta_out(declare_output_variable<Scalar>("mode_mixity")),
    _delta_init_out(declare_output_variable<Scalar>("damage_initiation_jump")),
    _delta_final_out(declare_output_variable<Scalar>("full_degradation_jump")),
    _delta_m_out(declare_output_variable<Scalar>("effective_displacement_jump")),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _G_Ic(declare_parameter<Scalar>("G_Ic", "mode_I_critical_energy_release_rate")),
    _G_IIc(declare_parameter<Scalar>("G_IIc", "mode_II_critical_energy_release_rate")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "power_law_exponent")),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity")),
    _dt(declare_parameter<Scalar>("dt", "time_step")),
    _criterion(options.get<std::string>("criterion")),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_displacement_jump(options.get<bool>("lag_displacement_jump")),
    _alpha(options.get<double>("regularization"))
{
  neml_assert(_criterion == "BK" || _criterion == "POWER_LAW",
              "criterion must be either BK or POWER_LAW");
}

void
BiLinearMixedModeTraction::request_AD()
{
  const std::vector<const VariableBase *> inputs = {&_jump, &_jump_old, &_d_old};
  _traction.request_AD(inputs);
  _d.request_AD(inputs);
  _beta_out.request_AD(inputs);
  _delta_init_out.request_AD(inputs);
  _delta_final_out.request_AD(inputs);
  _delta_m_out.request_AD(inputs);
}

void
BiLinearMixedModeTraction::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto eps = machine_precision(_jump.scalar_type());
  const auto zero = Scalar::zeros_like(_jump()(0));
  const auto one = Scalar::ones_like(_jump()(0));

  const auto mix_jump = _lag_mode_mixity ? _jump_old() : _jump();
  const auto mix_dn = mix_jump(0);
  const auto mix_ds = sqrt(mix_jump(1) * mix_jump(1) + mix_jump(2) * mix_jump(2) + eps);
  const auto safe_mix_dn = where(mix_dn > eps, mix_dn, one);
  const auto beta = where(mix_dn > 0.0, mix_ds / safe_mix_dn, zero);

  const auto delta_normal0 = _N / _K;
  const auto delta_shear0 = _S / _K;
  const auto delta_mixed =
      sqrt(delta_shear0 * delta_shear0 + beta * beta * delta_normal0 * delta_normal0 + eps);
  const auto delta_init_open = delta_normal0 * delta_shear0 * sqrt(1.0 + beta * beta) / delta_mixed;
  const auto delta_init = where(mix_dn > 0.0, delta_init_open, delta_shear0);

  const auto beta_sq_ratio = beta * beta / (1.0 + beta * beta);
  const auto delta_final_bk =
      2.0 / (_K * delta_init) * (_G_Ic + (_G_IIc - _G_Ic) * pow(beta_sq_ratio + eps, _eta));
  const auto Gc_mixed = pow(1.0 / _G_Ic, _eta) + pow(beta * beta / _G_IIc + eps, _eta);
  const auto delta_final_power =
      (2.0 + 2.0 * beta * beta) / (_K * delta_init) * pow(Gc_mixed, -1.0 / _eta);
  const auto delta_final_open = _criterion == "BK" ? delta_final_bk : delta_final_power;
  const auto delta_final_shear = std::sqrt(2.0) * 2.0 * _G_IIc / _S;
  const auto delta_final = where(mix_dn > 0.0, delta_final_open, delta_final_shear);

  const auto eff_jump = _lag_displacement_jump ? _jump_old() : _jump();
  const auto delta_n_pos = positive_part(eff_jump(0), _alpha);
  const auto delta_m =
      sqrt(delta_n_pos * delta_n_pos + eff_jump(1) * eff_jump(1) + eff_jump(2) * eff_jump(2) + eps);

  const auto delta_range = delta_final - delta_init;
  const auto d_mid =
      delta_final * (delta_m - delta_init) / (delta_m * where(delta_range > eps, delta_range, one));
  auto d_trial = where(delta_m < delta_init, zero, where(delta_m > delta_final, one, d_mid));
  d_trial = where(d_trial < _d_old, _d_old(), d_trial);
  const auto visc_factor = 1.0 + _viscosity / _dt;
  const auto d = (d_trial + _viscosity * _d_old / _dt) / visc_factor;

  const auto H = regularized_heaviside(_jump()(0), _alpha);
  const auto delta_active = Vec::fill(H * _jump()(0), _jump()(1), _jump()(2));
  const auto delta_inactive = Vec::fill(_jump()(0) - H * _jump()(0), zero, zero);

  _beta_out = beta;
  _delta_init_out = delta_init;
  _delta_final_out = delta_final;
  _delta_m_out = delta_m;
  _d = d;
  _traction = (1.0 - d) * _K * delta_active + _K * delta_inactive;
}
} // namespace neml2
