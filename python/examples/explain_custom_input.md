---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: neml2
  language: python
  name: python3
---

# Explain a custom NEML2 input file

Use this notebook when you want to explain a `.i` file that is not stored next to the example notebook.

## How to use this notebook

1. Activate your `neml2` environment.
2. Start Jupyter from the repository root with `jupyter lab`.
3. Open this notebook and select the `neml2` kernel.
4. Change `INPUT_FILE` below to the `.i` file you want to explain.
5. Run the cells from top to bottom.
6. Stop after `Inspect the generated prompt` if you only want the prompt.
7. Use the Argo section only if you want a live LLM response.

```{code-cell} ipython3
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import requests

import neml2
from neml2.reader import describe
from neml2.reader._syntax import SyntaxDB
from neml2.reader._llm import LLMClient
```

## Choose the input file

Set `INPUT_FILE` to the file you want to explain. An absolute path is the safest choice.

```{code-cell} ipython3
INPUT_FILE = Path("/absolute/path/to/your/model.i")
INCLUDE_PARAMS = True

if not INPUT_FILE.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

INPUT_FILE
```

Example:

```python
INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)
```

## Build the syntax database

This uses the bundled `neml2-syntax` tool. The notebook prefers the binary inside the imported `neml2` package and only falls back to `PATH` if needed.

```{code-cell} ipython3
syntax_exe = find_syntax_exe()
print(f"Using syntax tool: {syntax_exe}")

result = subprocess.run([str(syntax_exe)], check=True, capture_output=True, text=True)
with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as handle:
    handle.write(result.stdout)
    syntax_path = Path(handle.name)

syntax_db = SyntaxDB(syntax_path)
syntax_db.available
```

## Inspect the generated prompt

Run this first to verify parsing works before calling any LLM endpoint.

```{code-cell} ipython3
system_prompt, user_prompt = describe(INPUT_FILE, syntax_db, include_params=INCLUDE_PARAMS)

print("SYSTEM PROMPT")
print("=" * 80)
print(system_prompt)
print()
print("USER PROMPT")
print("=" * 80)
print(user_prompt[:4000])
```

## Optional: call Argo

This cell is only for use on Argonne's internal network. Leave it alone if you only want the prompt.

> The Argo request body requires the `"user"` field to contain your ANL username or API identity expected by your deployment.

```{code-cell} ipython3
ARGO_URL = "https://apps-dev.inside.anl.gov/argoapi/api/v1/resource/chat/"


def find_syntax_exe() -> Path:
    candidates = [
        Path(neml2.__file__).resolve().parent / "bin" / "neml2-syntax",
        Path.cwd() / "python" / "neml2" / "bin" / "neml2-syntax",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    found = shutil.which("neml2-syntax")
    if found:
        return Path(found)

    raise FileNotFoundError("Could not find neml2-syntax.")


class ArgoClient(LLMClient):
    def __init__(self, model: str, user: str | None = None):
        self._model = model
        self._user = user or os.environ["ARGO_API_KEY"]

    def complete(self, system: str, user: str) -> str:
        payload = {
            "user": self._user,
            "model": self._model,
            "messages": [
                {"role": "system", "content": system},
                {"role": "user", "content": user},
            ],
            "stop": [],
            "temperature": 0.1,
            "top_p": 0.9,
            "max_completion_tokens": 64000,
        }
        response = requests.post(ARGO_URL, data=json.dumps(payload), headers={"Content-Type": "application/json"})
        if response.status_code != 200:
            raise RuntimeError(
                f"Argo API request failed with status {response.status_code}: {response.text}"
            )
        return response.json()["response"]
```

```{code-cell} ipython3
# Uncomment to use the live Argo endpoint.
# client = ArgoClient(model="gpt54", user="your-anl-username")
# response = client.complete(system_prompt, user_prompt)
#
# from IPython.display import Markdown, display
# display(Markdown(response))
```
