---
name: spec-generate
description: Generate a compact, implementation-oriented Markdown spec for a NEML2 `Model` subclass (or a family of related Models). The spec follows a fixed three-section template — `# Input`, `# Algorithm`, `# Output` — plus a small `## NEML2 wiring` block that downstream scaffolders consume. Use this skill (a) when the user explicitly asks for a spec ("write a spec for X", "/spec-generate <Model>", "document the math of Y"), (b) when called from another skill's Step 0 (`add-model`, `add-domain`) to produce or refresh the spec that drives scaffolding, or (c) when reverse-engineering an existing Model's math into a spec for review or porting. Output path is always `design/<domain>/<Name>_spec.md` under the repo-root `design/` directory (gitignored, per-developer).
---

# spec-generate

Produce a compact spec that is the source of truth for a NEML2 `Model` subclass. The spec is implementation-oriented, not theoretical — its job is to let `add-model` (or a human) write `set_value` without going back to papers or guessing sign conventions.

The output has three core sections — `# Input`, `# Algorithm`, `# Output` — plus a `## NEML2 wiring` block and a few tail sections. Do not write theory notes, motivation, or literature reviews. Preserve the actual computation flow.

## When to use

- User explicitly asks for a spec (`/spec-generate`, "write a spec for `LinearIsotropicHardening`", "generate an implementation spec for the Voce law").
- Called from `add-model` Step 0 case A (user-attached non-Markdown spec needs to be normalized), case D (textbook formula needs a draft), or case E (user supplied a paper reference).
- Called from `add-domain` Step 0 to draft specs for each variant in a family.
- Reverse-engineering an existing Model: read its header + `.cxx` + base classes and write the spec retrospectively.

## When NOT to use

- The user inlined the formula and interface in the prompt and wants `add-model` to scaffold immediately — that's case C in `add-model`, no spec file needed.
- The request is for narrative / user-facing documentation — use `add-doc` instead.

## Step 1 — Read source

Locate and read the relevant code. For NEML2:

- `include/neml2/models/<domain>/` — header(s); look at `expected_options()` declaration and member layout
- `src/neml2/models/<domain>/` — implementation(s); read `expected_options()` body and the full `set_value(out, dout_din, d2out_din2)` including helper free functions it calls
- `tests/unit/models/<domain>/` — existing `.i` ModelUnitTest(s); they document the observable input/output contract
- `include/neml2/models/<domain>/<Base>.h` — base-class header if the model subclasses a domain interface (`YieldFunction`, `IsotropicHardening`, `FlowRule`, `Elasticity`, `Eigenstrain`, …); the base often declares variables the subclass inherits
- `include/neml2/tensors/<T>.h` — when an unfamiliar tensor type is involved, confirm batched semantics and component layout
- `doc/content/` — only if narrative docs already exist for the family

When the spec is for a *new* model with no implementation yet, read the **nearest existing model in the same domain** as a template, plus the user-provided reference (paper section, textbook chapter, equation number). Quote the reference inline in the `## Reference` tail section.

Always read the base-class chain. Do not guess `expected_options()`, member layout, or which derivatives `set_value` is expected to write.

## Step 2 — Output path

Write the spec to:

```
design/<domain>/<Name>_spec.md
```

- `<domain>` mirrors the relative path under `include/neml2/models/`. It may be nested, e.g. `solid_mechanics/elasticity/LinearIsotropicElasticity_spec.md`.
- `design/` lives at the **repo root** and is **gitignored** (per-developer, not committed). If `design/` is missing, create it. If the gitignore entry is missing, surface that to the user before writing — they may want to add it before the file lands.
- Create intermediate directories as needed.
- If a spec already exists at the target path, read it first and present a diff / proposed update rather than blindly overwriting.

## Step 3 — Classification rules

### `# Input`

Everything required before `set_value` runs. Split into four subsections:

#### Current inputs
- variables declared via `declare_input_variable<T>(...)`
- include the tensor type (`Scalar`, `SR2`, `Vec`, `Rot`, `R2`, `WR2`, `R3`, `R4`, `SSR4`, `WSR4`)
- include component order, units, and forced/state/old classification

#### State inputs
- previous-step / history quantities (anything declared on the `OLD` sub-axis or retrieved via the `_old` accessor)
- mark explicitly as `state input`

#### Parameters
- declared via `declare_parameter<T>(...)` and the matching `options.add<T>(...)`
- include type, default (if any), units, calibratable y/n

#### Control inputs
- booleans / enums / mode selectors set in `expected_options()` via plain `options.add<bool>` / `options.add<MultiEnumSelection>` / etc.
- e.g. lag flags, criterion selection, regularization switches, smoothing on/off

Rules:
- If it is read before compute → it is input
- If it affects branching → it is input
- Do NOT create a separate "Options" section — fold options into Parameters or Control inputs

### `# Algorithm`

Write the full computation flow of `set_value` in execution order. For each step:

1. what is computed (name the intermediate)
2. formula or pseudocode (TeX-style fine — `$...$`)
3. branch conditions
4. state update (if any)
5. source function (the primary is always `set_value`; helper free functions in `tensors/functions/` if used)

Must include:
- all intermediate tensors required for implementation
- all branch logic
- all guards (division by zero, sign split, eigenvalue clamping, …)
- all regularization (smoothing, log/exp clamping)
- all irreversibility / history logic
- first derivatives (`_out.d(_in)`) — spec out what `dout_din` should compute
- second derivatives (`d2out_din2`) only if implemented
- parameter derivatives if the model is calibratable

Do NOT collapse multi-step logic into a single equation. Do NOT omit guards.

### `# Output`

Anything produced by `set_value`:

#### Primary outputs
- variables declared via `declare_output_variable<T>(...)`
- include tensor type and units

#### Updated state outputs
- state variables this model writes (mark explicitly as `updated state output`)

#### Optional derived outputs
- only if exposed as separate output variables for debugging / downstream consumption

Do NOT include internal temporaries (those belong in `# Algorithm`, not `# Output`).

## Required structure

```markdown
# Input

## Current inputs
...

## State inputs
...

## Parameters
...

## Control inputs
...

---

# Algorithm

1. Compute ...
   - formula:
   - branch:
   - state update:
   - source:

2. Compute ...
...

---

# Output

## Primary outputs
...

## Updated state outputs
...

## Optional derived outputs
...

---

## NEML2 wiring
- Class name: <PascalCase>
- Domain (path under `include/neml2/models/`): <domain>
- Base class: <Model | YieldFunction | IsotropicHardening | FlowRule | Elasticity | Eigenstrain | ...>
- Test file: tests/unit/models/<domain>/<Name>.i
- Force-link: not needed (umbrella `models/` already linked)

## Source files read
- ...

## Reference
- paper / textbook / equation number / URL — required for new models

## Open issues
- ...
```

The `## NEML2 wiring` block is mandatory — `add-model` and `add-domain` consume it to scaffold without re-prompting the user for class name, domain, and base class. Without it, downstream skills will stop and ask, defeating the purpose of having a spec.

## Hard rules

- Prefer existing NEML2 variable names (`equivalent_plastic_strain`, `mandel_stress`, `cauchy_stress`, …) over renamings; only rename when the user explicitly does
- Preserve execution order of `set_value`
- Do not hide branch logic
- Do not omit intermediate variables
- Treat anything on the `OLD` axis or accessed via `_old` as a state input
- Treat updated history as state output
- Treat ALL flags / enums / modes as input (under Control inputs)
- Always include the `## NEML2 wiring` block
- Always include `## Reference` for new models — `set_value` derivations must be traceable to a citation. For reverse-engineered specs (existing implementation), `## Reference` may say "reverse-engineered from src/neml2/models/<domain>/<Name>.cxx" plus any in-code comments / paper citations found there.

## After writing

Print a short message to the user:

```
Wrote spec to: design/<domain>/<Name>_spec.md

Next:
- Review the spec — especially the `# Algorithm` section and the `## NEML2 wiring` block.
- If correct, run `/add-model <Name>` (or `/add-domain <domain>`); Step 0 case B will pick up this spec automatically.
- The `design/` directory is gitignored and per-developer — this spec is not committed.
```

If `design/` was not in `.gitignore` when the spec was written, add a line to the message asking the user to add `design/` to `.gitignore` before any commit.
