import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import requests

import neml2
from neml2.reader import describe
from neml2.reader._llm import LLMClient
from neml2.reader._syntax import SyntaxDB

INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)

INCLUDE_PARAMS = True
ARGO_URL = "https://apps-dev.inside.anl.gov/argoapi/api/v1/resource/chat/"
ARGO_TIMEOUT_SECONDS = float(os.environ.get("ARGO_TIMEOUT_SECONDS", "120"))


def find_syntax_exe() -> Path:
    candidates = [
        Path(neml2.__file__).resolve().parent / "bin" / "neml2-syntax",
        Path(__file__).resolve().parents[2] / "python" / "neml2" / "bin" / "neml2-syntax",
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
        try:
            response = requests.post(
                ARGO_URL,
                data=json.dumps(payload),
                headers={"Content-Type": "application/json"},
                timeout=ARGO_TIMEOUT_SECONDS,
            )
        except requests.Timeout as exc:
            raise RuntimeError(
                "Argo API request timed out after "
                f"{ARGO_TIMEOUT_SECONDS} seconds. "
                "The server may be queued or taking too long to generate a response. "
                "Try increasing ARGO_TIMEOUT_SECONDS or reducing the prompt size."
            ) from exc
        except requests.RequestException as exc:
            raise RuntimeError(f"Argo API request failed: {exc}") from exc

        if response.status_code != 200:
            raise RuntimeError(
                f"Argo API request failed with status {response.status_code}: {response.text}"
            )

        body = response.json()
        if "response" not in body:
            raise RuntimeError(f"Argo API response did not contain 'response': {body}")

        return body["response"]


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    syntax_exe = find_syntax_exe()
    print(f"Using syntax tool: {syntax_exe}")

    result = subprocess.run([str(syntax_exe)], check=True, capture_output=True, text=True)
    with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as handle:
        handle.write(result.stdout)
        syntax_path = Path(handle.name)

    try:
        syntax_db = SyntaxDB(syntax_path)

        system_prompt, user_prompt = describe(INPUT_FILE, syntax_db, include_params=INCLUDE_PARAMS)
        print(
            "Prepared prompt "
            f"(system={len(system_prompt)} chars, user={len(user_prompt)} chars)."
        )

        argo_user = os.environ.get("ARGO_USER", "").strip()
        if not argo_user or argo_user.lower() in {
            "your-anl-username",
            "your-argo-username",
            "set-me",
        }:
            raise RuntimeError(
                "Set ARGO_USER to your actual Argo username before running the LLM call."
            )

        client = ArgoClient(model="gpt54", user=argo_user)
        print(
            f"Sending request to Argo at {ARGO_URL} "
            f"with timeout={ARGO_TIMEOUT_SECONDS}s..."
        )
        response = client.complete(system_prompt, user_prompt)
        print(response)
    finally:
        syntax_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
