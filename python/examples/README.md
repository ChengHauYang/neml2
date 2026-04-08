# Python examples and reader tests

Use one Python environment for both the notebooks in `python/examples` and the tests in `python/tests`.

## Recommended setup

Run these commands from the repository root:

### `venv`

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install 'torch==2.8.0'
python -m pip install -e '.[dev]' --no-deps
python -m pip install pyzag==1.1.1 pytest pytest-xdist nbmake pybind11-stubgen PyYAML matplotlib numpy scipy pandas loguru pyro-ppl tqdm nbconvert nbformat livereload jupytext pre-commit
python -m pip install jupyterlab ipykernel requests
python -m ipykernel install --user --name neml2 --display-name neml2
```

### `conda`

```bash
conda create -n neml2 python=3.12 -y
conda activate neml2
python -m pip install --upgrade pip
python -m pip install 'torch==2.8.0'
python -m pip install -e '.[dev]' --no-deps
python -m pip install pyzag==1.1.1 pytest pytest-xdist nbmake pybind11-stubgen PyYAML matplotlib numpy scipy pandas loguru pyro-ppl tqdm nbconvert nbformat livereload jupytext pre-commit
python -m pip install jupyterlab ipykernel requests
python -m ipykernel install --user --name neml2 --display-name neml2
```

Notes:

- `torch==2.8.0` matches the version used by this repository's CI and avoids binary compatibility problems with the compiled extension.
- `-e '.[dev]' --no-deps` installs `neml2` in editable mode without upgrading `torch` behind your back.
- `jupyterlab`, `ipykernel`, and `requests` are added explicitly for the notebook workflow.
- `openai` is needed only if you want to use the OpenAI API example below.
- The editable install builds the Python extension and the bundled tools such as `neml2-syntax`.

## Activate the environment

### `venv`

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
source .venv/bin/activate
```

### `conda`

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
conda activate neml2
```

## Run the reader test

Run tests from the repository root so imports and relative paths behave consistently:

```bash
pytest -q python/tests/reader/test_explain_workflow.py
```

Do not run the test as `python test_explain_workflow.py`; it is a `pytest` test module, not a standalone script.

## Run the notebook

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
jupyter lab
```

Then open [explain_custom_input.md](/Users/chenghau.yang/Documents/Package/neml2/python/examples/explain_custom_input.md), set `INPUT_FILE` to your `.i` file, and run the cells.

## User guide for `explain_custom_input.md`

1. Activate the environment.

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
source .venv/bin/activate
```

2. Start Jupyter.

```bash
jupyter lab
```

3. Open [explain_custom_input.md](/Users/chenghau.yang/Documents/Package/neml2/python/examples/explain_custom_input.md).

4. In the `Choose the input file` cell, replace `INPUT_FILE` with your actual `.i` file path. Example:

```python
INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)
```

5. Select the `neml2` kernel if Jupyter asks for a kernel.

6. Run the cells from top to bottom:

- import cell
- `Choose the input file`
- `Build the syntax database`
- `Inspect the generated prompt`

7. If you only want the structured prompt, stop there.

8. If you want the LLM explanation too, configure the Argo cell and then uncomment the final cell.

Expected result:

- the prompt cell prints the system prompt and the first part of the user prompt
- the optional Argo cell returns a natural-language explanation of the input file

## Runnable Python example for LLM explanation

Use [run_explain.py](/Users/chenghau.yang/Documents/Package/neml2/python/examples/run_explain.py).

Run this from the repository root after activating the environment:

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
source .venv/bin/activate
export ARGO_USER='your-anl-username'
python python/examples/run_explain.py
```

What to edit in [run_explain.py](/Users/chenghau.yang/Documents/Package/neml2/python/examples/run_explain.py):

```python
INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)
```

Notes:

- Replace `"your-anl-username"` with the value required by your Argo deployment.
- Set `ARGO_USER` before running the script.
- If you only want to inspect the prompt, modify [run_explain.py](/Users/chenghau.yang/Documents/Package/neml2/python/examples/run_explain.py) to print `system_prompt` and `user_prompt` instead of calling `client.complete(...)`.
- Run the script from the repository root so the local editable install and tools are picked up consistently.

## Runnable Python example with OpenAI

OpenAI recommends the Responses API for new integrations. This repo now includes [run_explain_openai.py](/Users/chenghau.yang/Documents/Package/neml2/python/examples/run_explain_openai.py), which uses that API.

Install the client once in your active environment:

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
source .venv/bin/activate
# python -m pip install openai # do it once
```

Run it like this:

```bash
cd /Users/chenghau.yang/Documents/Package/neml2
source .venv/bin/activate
# export OPENAI_API_KEY='your-openai-api-key' # we put this in .zshrc
python python/examples/run_explain_openai.py
```

What to edit in [run_explain_openai.py](/Users/chenghau.yang/Documents/Package/neml2/python/examples/run_explain_openai.py):

```python
INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)

MODEL = "gpt-5"
```

Model choice:

- Use `"gpt-5"` as the default general-purpose model.
- Use `"gpt-5-codex"` if you specifically want the Codex-optimized model in the Responses API.

Notes:

- `OPENAI_API_KEY` is required.
- If you only want to inspect the prompt, print `system_prompt` and `user_prompt` instead of calling `client.responses.create(...)`.
- Run the script from the repository root so the local editable install and tools are picked up consistently.

## Version control for Jupyter notebook examples

See [these instructions](../../doc/content/tutorials/contributing.md#jupyter-notebooks) for details.
