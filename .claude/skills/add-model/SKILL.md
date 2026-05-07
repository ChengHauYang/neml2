---
name: add-model
description: Scaffold a new NEML2 `: public Model` subclass — header, registered source, and `.i` ModelUnitTest input. Use this skill whenever the user wants to add, create, or scaffold any Model-family class: constitutive models (yield functions, hardening laws, flow rules, elasticity, eigenstrains), or non-constitutive Model subclasses (converters, parameter maps, smoothing functions, geometric helpers). Trigger on phrases like "add a new isotropic hardening model", "create a model that maps X to Y", "scaffold a YieldFunction subclass", "add a Foo model under solid_mechanics", even when the user doesn't say the word "Model" — if it derives from `Model`, this skill applies. Spec-driven: begins with a Step 0 that consumes design specs from `design/` (per-developer, gitignored) when present, drafts a spec from textbook references for simple cases, and stops to request a reference for complex / research-level requests.
---

# add-model

Scaffold a new `: public Model` subclass in NEML2 — header, source with the registration macro, and a `ModelUnitTest` input file. CMake auto-discovers all three via `file(GLOB_RECURSE)`, so no CMake edits are needed when adding under an existing submodule.

The workflow has four steps: **Step 0** — settle the spec; **Step 1** — confirm scaffolding parameters; **Step 2** — generate the files; **Step 3** — hand off to the user. Do not skip Step 0: scaffolding code without a settled spec produces typed-out guesses that compile but encode the wrong math.

## When to use

Use this skill for any class that derives from `Model`: constitutive models, but also converters, parameter maps, smoothing functions, geometric helpers — anything in `include/neml2/models/`. Typical phrasing: "add a new model that …", "scaffold a Foo model under solid_mechanics", "create a YieldFunction subclass".

Do **not** use this skill for:
- Adding a new tensor class (under `include/neml2/tensors/`) — different patterns and different test binary (`tensor_tests`).
- Adding a `Driver`, `Solver`, or `WorkScheduler` — different base class, different test setup.
- Creating a brand-new top-level submodule library (a sibling of `models/`, `solvers/`, `drivers/`, etc.). That's a rare event that needs a `neml2_add_submodule(...)` call, a `link_anchor.cxx`, and a force-link symbol declared in `include/neml2/neml2.h`. Surface that to the user instead of attempting it.

Adding a brand-new *domain directory* under an existing submodule (e.g. `models/electrochemistry/`) is fine — `models/` is already a submodule and its CMake glob is recursive, so any new subdirectory is picked up automatically.

## Step 0: Design Spec (do this first)

Before scaffolding, decide what spec drives the implementation. Specs are the source of truth for math, sign conventions, units, and variable interfaces — without one, scaffolded code is a typed-out guess. Supported spec formats: `.md` (preferred), `.pdf`, `.txt`, `.html`, `.htm`. Convert non-Markdown specs to Markdown before writing any code.

This skill **consumes** specs; it does not own the spec template. The bundled sibling skill `spec-generate` (in `.claude/skills/spec-generate/`) owns the spec template, output path convention (`design/<domain>/<Name>_spec.md`), and writing protocol. Whenever this skill needs to produce or refresh a spec — cases A, D, E below — invoke `spec-generate` via the Skill tool rather than improvising. The Spec format summary further down is for orientation only; `spec-generate` is the authoritative source.

### Decision tree

**A. User attached spec file(s) in this turn** — read all of them before writing any code. Unsupported extension → stop and ask for `.md` / `.pdf` / `.txt` / `.html` / `.htm`. Markdown is source of truth; for non-Markdown that isn't trivially mappable, invoke `spec-generate` to normalize it to `design/<domain>/<Name>_spec.md` first, then implement from that.

**B. No attached spec, but `design/` exists** — search recursively for `.md` files (`spec-generate` always writes Markdown, so other extensions only appear if a contributor put them there manually; treat them as case A inputs needing normalization).
- Single relevant match → use it as source of truth. If the spec gives analytical Jacobians, implement them exactly as written.
- Multiple specs in the same directory → read **all** of them before writing code; treat the directory as a candidate model family. Compare their Input/Output sections:
  - Same primary input type **and** same primary output type → that shared interface justifies an abstract base class. Introduce the base first, then implement each variant in the same scope by default.
  - Divergent inputs/outputs → implement separate independent classes.
- Specs in different directories → implement only the closest match; don't sprawl.

**C. No spec found, request is straightforward** — single model, clear math, unambiguous interface (e.g., user inlined the formula, named the input/output types, gave a class name). Proceed from the prompt + the nearest existing similar model. No spec file needs to be created.

**D. No spec found, request is simple but math/convention isn't inlined** — invoke `spec-generate` to draft a Markdown spec from your own knowledge **only if the formula is textbook / industry consensus** (linear elasticity, J2, Voce hardening, linear-softening CZM, etc.). `spec-generate` writes to `design/<domain>/<Name>_spec.md` so case B picks it up on a future run. Cite the reference (textbook section, well-known paper, Wikipedia entry) inline in the spec's `## Reference` section. **Stop and show the draft path + contents to the user — wait for explicit confirmation before writing any code.** If you can't cite a definite source for the form you're using, escalate to case E.

**E. No spec found, request is complex or research-flavored** — multiple plausible variants (e.g., "traction separation law", "viscoelastic kernel", "yield surface family"), or a model family, or non-trivial math with subtle sign / mode-mixity / reference-configuration conventions. **Stop. Do not guess.** Ask the user to either:
- place a supported spec file under `design/<domain>/`, or
- provide a reference (paper + equation number, textbook section) so `spec-generate` can produce a spec against a fixed source.

LLM knowledge is unreliable on these — sign conventions, normal-vs-tangential coupling, and reference-configuration choices vary between papers and silently encode the wrong physics. For a model family in case E, do not pick one variant unilaterally; surface the choices and let the user pick (or commission specs for the full family via `spec-generate`).

### `design/` directory convention

`design/` lives at the **repo root** and is **per-developer, gitignored** — each contributor keeps their own working specs locally; specs are not committed. Only mention this to the user when `design/` is **not** in `.gitignore` (so they can add it before any spec lands in a commit). Otherwise stay silent on the gitignore status.

### Spec format

See `spec-generate` for the authoritative template. Do not improvise from memory. The one piece this skill relies on is the mandatory `## NEML2 wiring` block (class name, domain, base class) — Step 1 reads it directly to skip re-prompting.

### Scope rule

Scope narrowing is **opt-in, never implicit**. If multiple same-directory specs or a complex request implies a model family, the default scope is the **full family**. Do not silently stop after the easiest variant. To defer a variant, mark it `DEFERRED` explicitly and get user approval before writing code.

### Before moving to Step 1, state to the user:

- which spec file(s) were used (including any draft generated in this turn), or that none was needed because the request is straightforward (case C);
- the planned class hierarchy — base + variants, or independent classes;
- the implementation scope for this turn — full family or an explicitly-named subset/MVP;
- the exact list of classes that will be scaffolded this turn, named explicitly;
- any `DEFERRED` variant with a one-line reason, or state that there are none.

## Step 1: Confirm scaffolding parameters

Once Step 0 has produced or pointed at a settled spec (or established that the request is straightforward enough not to need one), gather the parameters below.

**If a spec exists**, the four parameters below are already in the spec's `## NEML2 wiring` block (class name, domain, base class) and `# Input` / `# Output` sections (variables and types). Read them from there — do not re-prompt the user. Only ask if the spec's `## NEML2 wiring` block is missing or ambiguous.

**If no spec was needed (case C)**, confirm:

1. **Class name** in PascalCase (e.g., `MyVoceHardening`).
2. **Domain** — must be one of the existing directories under `include/neml2/models/`. Verify with `ls include/neml2/models/` if unsure. Common choices: `solid_mechanics`, `solid_mechanics/elasticity`, `solid_mechanics/crystal_plasticity`, `crystallography`, `chemical_reactions`, `phase_field_fracture`, `porous_flow`, `finite_volume`, `common`. The directory may be one or two levels deep.
3. **Base class** — default `Model`. Many solid-mechanics classes subclass a domain interface instead so that consumers can swap implementations: yield functions subclass `YieldFunction`, isotropic-hardening laws subclass `IsotropicHardening`, flow rules subclass `FlowRule`, etc. If the user describes the model in those terms, prefer the domain base. When unsure, grep `class .* : public ` in the target header directory and follow the dominant pattern.
4. **Input/output variables** — what tensor types come in and go out. Common types: `Scalar`, `Vec`, `Rot`, `R2`, `SR2`, `WR2`, `R3`, `R4`, `SSR4`, `WSR4`. Users often specify these implicitly ("maps strain to stress" → input `SR2`, output `SR2`).

Use the variable names the user gives, or the names that match the domain convention (e.g., `stress`/`strain` in solid mechanics) — input and output names appear in user-facing input files, so guessing leads to renames later. Ask if you can't tell what the convention should be.

## Step 2: Generate the files

For a class `<Name>` in domain `<Domain>` (which may contain `/`, e.g. `solid_mechanics/elasticity`), generate three files using the templates below. The templates cover the common shape: `: public Model`, one input, one output, no parameters.

For anything more specialized — multiple inputs/outputs of different tensor kinds, parameterized models, models that subclass a domain interface like `YieldFunction` or `Elasticity`, or domains with their own conventions like `crystal_plasticity` — read one nearby existing model in the same directory first and adapt its pattern. The references at the bottom of this skill list good starting points.

### License header (both `.h` and `.cxx`)

Every C++ file must begin with the project's standard MIT copyright header — `scripts/check_copyright.py` does an exact textual check, so paraphrasing fails CI. The header is the contents of the repo-root `LICENSE` file with each line prefixed by `// ` (and `//` for the blank lines between paragraphs). Read `LICENSE` and reproduce it that way. Match the year in `LICENSE` exactly; don't bump it for new files.

If unsure of the format, mirror the top of `src/neml2/models/common/RotationMatrix.cxx` (or any other `.cxx`/`.h` in the tree).

### 1. Header — `include/neml2/models/<Domain>/<Name>.h`

```cpp
// <project MIT header — see "License header" section above>

#pragma once

#include "neml2/models/Model.h"   // or e.g. ".../solid_mechanics/IsotropicHardening.h" for a non-Model base

namespace neml2
{
class <InputType>;   // forward declarations of tensor types used only by reference/pointer
class <OutputType>;

/**
 * @brief <One-line summary the user gives. The first sentence becomes the syntax-doc summary.>
 */
class <Name> : public <Base>
{
public:
  static OptionSet expected_options();

  <Name>(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  const Variable<<InputType>> & _<input_name>;

  Variable<<OutputType>> & _<output_name>;
};
} // namespace neml2
```

### 2. Source — `src/neml2/models/<Domain>/<Name>.cxx`

```cpp
// <project MIT header — see "License header" section above>

#include "neml2/models/<Domain>/<Name>.h"
#include "neml2/tensors/<InputType>.h"
#include "neml2/tensors/<OutputType>.h"

namespace neml2
{
register_NEML2_object(<Name>);

OptionSet
<Name>::expected_options()
{
  OptionSet options = <Base>::expected_options();
  options.doc() = "<Same one-line summary used in the header docstring.>";

  options.add_input("<input_name>", "<Short description of the input.>");
  options.add_output("<output_name>", "<Short description of the output.>");

  return options;
}

<Name>::<Name>(const OptionSet & options)
  : <Base>(options),
    _<input_name>(declare_input_variable<<InputType>>("<input_name>")),
    _<output_name>(declare_output_variable<<OutputType>>("<output_name>"))
{
}

void
<Name>::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // TODO: implement the forward operator.
  if (out)
    _<output_name> = /* function of _<input_name>() */;

  // TODO: implement the first derivative. Until this is filled in,
  // the unit test must set check_first_derivatives = false.
  if (dout_din)
    _<output_name>.d(_<input_name>) = /* d output / d input */;
}
} // namespace neml2
```

Notes on the source template:
- `register_NEML2_object(<Name>);` is mandatory. Without it the type is silently unknown to the factory at runtime.
- `_<input_name>()` (with parentheses) returns the underlying tensor of the input variable; `_<output_name> = …` and `_<output_name>.d(_<input_name>) = …` write the forward value and the derivative respectively.
- If the model has parameters (calibratable scalars/tensors), declare them with `declare_parameter<...>("name")` and add a matching `options.add<...>("name")` in `expected_options()`. Look at `LinearIsotropicHardening` or `FredrickArmstrongPlasticHardening` for parameterized examples.
- If `set_value` provides analytic first derivatives, NEML2 will cross-check them against finite differencing in the unit test (see step 3). If you can't yet provide them, set `check_first_derivatives = false` in the `.i` file.

### Naming convention for rates, histories, and residuals — DERIVE, don't request

When a variable's name is the *rate*, *history*, or *residual* form of another variable, **derive the name from the base variable using the helper functions** (`rate_name`, `history_name`, `residual_name` from `NEML2Object`) — do **not** add a separate `add_input`/`add_output` option that lets the user name it directly.

These helpers respect the global `[Settings]` block, where users can override `rate_prefix`, `rate_suffix` (default `_rate`), `history_separator` (default `~`), `residual_prefix`, and `residual_suffix` (default `_residual`). A model that hardcodes the literal name `viscous_strain_rate` works correctly only by accident with default settings — change `rate_suffix = '.dot'` and the model and the rest of the framework no longer agree. Worse, a model that requests the rate name as a separate input option lets the user *type a different name from the base variable* — the model and `BackwardEulerTimeIntegration` would silently disagree on which buffer is the rate, and the integration would degenerate (the rate input falls back to zero) without any error.

Pattern from `BackwardEulerTimeIntegration`:

```cpp
// In expected_options(): only add_input the BASE variable.
options.add_input("variable", "Variable being integrated");
options.add_input("time", "t", "Time");

// In the constructor: derive the rate, history, and residual names from the base.
_s(declare_input_variable<T>("variable")),
_sn(declare_input_variable<T>(history_name(_s.name(), /*nstep=*/1))),     // x~1
_rate(declare_input_variable<T>(rate_name(_s.name()))),                   // x_rate
_r(declare_output_variable<T>(residual_name(_s.name())))                  // x_residual
```

Same pattern for outputs. From `IsotropicHardeningStaticRecovery`:

```cpp
// expected_options(): only add the base variable.
options.add_input("isotropic_hardening", "Isotropic hardening variable");
// no add_output for the rate

// constructor: rate output name is derived.
_h(declare_input_variable<Scalar>("isotropic_hardening")),
_h_dot(declare_output_variable<Scalar>(rate_name(_h.name())))
```

**Don't declare a base variable as input just to give yourself a name to derive from.** If your model's equation does not actually use the value of the underlying state variable, declaring it as an input lies about the model's dependency graph — the dependency resolver will think this model needs the variable when it doesn't, which can force a spurious order on `ComposedModel` or block evaluation when the variable isn't yet defined. In that case, expose the rate as a freestanding `add_output` option and let the user name it explicitly. They are responsible for choosing a name that lines up with `rate_name` of whatever state variable the time integrator owns.

`MaxwellViscoelasticity` is the canonical example: the dashpot equation `dε_v/dt = σ/η` depends only on stress. Declaring `viscous_strain` as an unused input would falsely couple Maxwell to that variable. So Maxwell takes only `stress`, and `viscous_strain_rate` is a freestanding `add_output` option with a default name that happens to match `rate_name(viscous_strain)` at default settings — but the option's docstring should make explicit that the user is on the hook for keeping it consistent with their integrator. `PlasticFlowRate` is similar: `flow_rate` is the rate of a "consistency parameter" that doesn't exist as a state variable anywhere, so it's a freestanding output option.

Quick checklist when declaring an input or output:

- Does this variable's *value* enter the equation, AND is it the rate/history/residual of another input variable that's also used? → derive the name with `rate_name`/`history_name`/`residual_name`, no `add_<input/output>` for it.
- Is this freestanding — either no underlying base variable exists (e.g. `flow_rate`, `stress`, `strain`), or the model genuinely doesn't use the base variable's value (e.g. Maxwell's `viscous_strain_rate`)? → `add_<input/output>` with the literal name, document the user's responsibility in the option docstring if it's a rate-like quantity that needs to align with an integrator.

Don't add a literal-name `add_<input/output>` just because it's easier — that's the failure mode this convention exists to prevent. Reserve it for genuine standalone variables.

### 3. Unit-test input — `tests/unit/models/<Domain>/<Name>.i`

The test binary `unit_tests` auto-discovers every `.i` file under `tests/unit/models/`. There is no need to add the file anywhere else.

```hit
[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_<TENSOR_KIND>_names = '<input_name>'
    input_<TENSOR_KIND>_values = '<input_value_tensor>'
    output_<TENSOR_KIND>_names = '<output_name>'
    output_<TENSOR_KIND>_values = '<expected_output_tensor>'
    # Set to false if analytic derivatives are not yet implemented:
    check_first_derivatives = true
    derivative_abs_tol = 1e-6
    derivative_rel_tol = 1e-4
  []
[]

[Tensors]
  [<input_value_tensor>]
    type = Fill<TensorType>
    values = '<comma-or-space-separated values>'
  []
  [<expected_output_tensor>]
    type = Fill<TensorType>
    values = '<comma-or-space-separated values>'
  []
[]

[Models]
  [model]
    type = <Name>
    <input_name> = '<input_name>'
    <output_name> = '<output_name>'
    # Add any parameters the model declares here.
  []
[]
```

`<TENSOR_KIND>` is the suffix used by `ModelUnitTest`'s option groups. The common ones are `Scalar`, `SR2`, `R2`, `Rot`, `Vec`, `SSR4`. If the model has multiple inputs/outputs of different kinds, repeat the `input_<KIND>_names` / `input_<KIND>_values` pair (and likewise for outputs) for each kind.

If you don't have known reference values yet, leave `<expected_output_tensor>` filled with placeholder zeros and tell the user to substitute real numbers before committing.

## Step 3: Hand off to the user

Print a short next-steps message to the user:

```
Scaffolded:
  include/neml2/models/<Domain>/<Name>.h
  src/neml2/models/<Domain>/<Name>.cxx
  tests/unit/models/<Domain>/<Name>.i

Build:
  cmake --build --preset dev

Run only this test:
  ./build/dev/tests/unit/unit_tests models -c <Domain>/<Name>.i

CMake auto-discovers the new files via file(GLOB_RECURSE) — no CMakeLists.txt edits are needed.

Before committing:
  1. Fill in the set_value forward op and derivatives in <Name>.cxx.
  2. Replace the placeholder tensor values in <Name>.i with actual reference data.
  3. Keep check_first_derivatives = true once derivatives are correct (NEML2 cross-checks them against finite differencing).
  4. clang-format will rewrite the files on pre-commit; that's expected.
```

## Easy mistakes to verify against

These are the failure modes worth a second look — most of them compile fine and surface only at runtime or test time:

- **Missing `register_NEML2_object(<Name>);`** in the `.cxx`. Without it the type loads but the factory doesn't know about it; the `.i` test fails with "unknown type". The macro must be at namespace scope, outside any function.
- **Tensor types used by value but only forward-declared.** Anything used as a value (constructed, returned, arithmetic on) needs a full `#include "neml2/tensors/<T>.h"` in the `.cxx`. Forward declaration in the `.h` is fine for member types held by reference.
- **`.i` options that don't match `expected_options()`** — every input/output declared in `expected_options()` must be set in the `[Models]` block, and vice versa. The HIT parser fails fast and the message is clear, but it's easy to miss when copying templates.

## Reference: canonical example to follow

A minimal-but-complete `: public Model` implementation in the codebase to mirror is `RotationMatrix`:
- `include/neml2/models/common/RotationMatrix.h`
- `src/neml2/models/common/RotationMatrix.cxx`

For a parameterized Model with multiple inputs and a derivative w.r.t. parameters, see:
- `include/neml2/models/solid_mechanics/LinearIsotropicHardening.h`

For a Model that subclasses a domain base (not `Model` directly), see:
- `include/neml2/models/solid_mechanics/elasticity/LinearIsotropicElasticity.h`

The official documentation walkthrough is `doc/content/tutorials/contributing.md` and the [extension](doc/content/tutorials/extension/) tutorial set.
