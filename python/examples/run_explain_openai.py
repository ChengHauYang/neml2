import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import neml2
from openai import OpenAI

from neml2.reader import describe
from neml2.reader._syntax import SyntaxDB

INPUT_FILE = Path(
    "/Users/chenghau.yang/Documents/Package/app/SBM-Poisson/test/tests/"
    "sbm_solid_mechanics_interface/CPFE-easy-NEML2/approx_kinematics_neml2.i"
)

INCLUDE_PARAMS = True
MODEL = "gpt-5"


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


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    if "OPENAI_API_KEY" not in os.environ:
        raise RuntimeError("Set OPENAI_API_KEY before running this script.")

    syntax_exe = find_syntax_exe()
    print(f"Using syntax tool: {syntax_exe}")
    print(f"Using OpenAI model: {MODEL}")

    result = subprocess.run([str(syntax_exe)], check=True, capture_output=True, text=True)
    with tempfile.NamedTemporaryFile("w", suffix=".yml", delete=False) as handle:
        handle.write(result.stdout)
        syntax_path = Path(handle.name)

    try:
        syntax_db = SyntaxDB(syntax_path)
        system_prompt, user_prompt = describe(INPUT_FILE, syntax_db, include_params=INCLUDE_PARAMS)

        client = OpenAI()
        response = client.responses.create(
            model=MODEL,
            instructions=system_prompt,
            input=user_prompt,
        )
        print(response.output_text)
    finally:
        syntax_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
