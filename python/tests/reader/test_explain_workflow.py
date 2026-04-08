# Copyright 2024, UChicago Argonne, LLC
# All Rights Reserved
# Software Name: NEML2 -- the New Engineering material Model Library, version 2
# By: Argonne National Laboratory
# OPEN SOURCE LICENSE (MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import textwrap
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

from neml2.reader import describe, explain
from neml2.reader._syntax import SyntaxDB

SYNTAX_CONTENT = textwrap.dedent(
    """\
    neml2::LinearIsotropicElasticity:
      section: Models
      doc: |-
        Relate elastic strain to stress for linear isotropic material.
      youngs_modulus:
        type: double
        ftype: PARAMETER
        doc: Young's modulus.
        suppressed: 0
        value:
      poisson_ratio:
        type: double
        ftype: PARAMETER
        doc: Poisson's ratio.
        suppressed: 0
        value:
""",
)

INPUT_CONTENT = textwrap.dedent(
    """\
    [Models]
      [elastic]
        type = LinearIsotropicElasticity
        youngs_modulus = 210000
        poisson_ratio = 0.3
      []
    []
""",
)


class MockClient:
    def __init__(self):
        self.calls = []

    def complete(self, system: str, user: str) -> str:
        self.calls.append((system, user))
        return "Mock explanation"


@pytest.fixture
def syntax_db(tmp_path):
    p = tmp_path / "syntax.yml"
    p.write_text(SYNTAX_CONTENT)
    return SyntaxDB(p)


@pytest.fixture
def input_file(tmp_path):
    p = tmp_path / "custom_model.i"
    p.write_text(INPUT_CONTENT)
    return p


def test_describe_accepts_path_objects(input_file: Path, syntax_db: SyntaxDB):
    system, user = describe(input_file, syntax_db, include_params=True)
    assert "NEML2 computational materials science library" in system
    assert "custom_model.i" not in user
    assert 'Model "elastic"' in user
    assert "youngs_modulus = 210000" in user


def test_explain_uses_mock_client(input_file: Path, syntax_db: SyntaxDB):
    client = MockClient()
    response = explain(input_file, syntax_db, client, include_params=False)

    assert response == "Mock explanation"
    assert len(client.calls) == 1

    system, user = client.calls[0]
    assert "NEML2 computational materials science library" in system
    assert "youngs_modulus = 210000" not in user
    assert "LinearIsotropicElasticity" in user
