"""Every ``run:`` block in a workflow has to be valid shell.

Nothing checks a workflow's shell until the job runs. For the workflows that run
on a schedule, that means a quoting mistake sits there silently and costs a whole
night -- and PR CI never exercises them, so review will not catch it either.

Not hypothetical: a stray backslash-quote at the end of one line in
``parity-nightly.yml`` killed the job at its licence step.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from ruamel.yaml import YAML

WORKFLOWS = sorted((Path(__file__).resolve().parents[1] / ".github" / "workflows").glob("*.yml"))


def read(path: Path) -> dict:
    """Parse a workflow, with CRLF normalised away.

    Git stores these files with LF, which is what a runner gets. A Windows
    working copy has CRLF, and bash reports a stray carriage return as a syntax
    error -- so without this the check would fail on every Windows machine and
    pass in CI, which is the wrong way round.
    """
    return YAML(typ="safe").load(path.read_text(encoding="utf-8").replace("\r\n", "\n"))


def run_blocks(path: Path):
    """Yield (step label, script) for every ``run:`` in the workflow."""
    document = read(path)
    for job_name, job in (document.get("jobs") or {}).items():
        for index, step in enumerate(job.get("steps") or []):
            script = step.get("run")
            if script:
                yield f"{path.name}:{job_name}:{step.get('name') or index}", script


def strip_expressions(script: str) -> str:
    """Replace ``${{ ... }}`` with a literal, so this checks shell not templating.

    Actions substitutes these before bash sees them. Leaving them in would make
    bash parse ``{{`` as shell; dropping them entirely would let an unquoted
    empty expansion look like valid syntax when it is not.
    """
    lines = []
    for line in script.split("\n"):
        while "${{" in line and "}}" in line:
            head, rest = line.split("${{", 1)
            line = head + "PLACEHOLDER" + rest.split("}}", 1)[1]
        lines.append(line)
    return "\n".join(lines)


def test_the_workflow_directory_was_actually_found():
    """Otherwise every test below passes by having nothing to check."""
    assert WORKFLOWS, "no workflow files found"


@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_workflow_is_valid_yaml(path: Path):
    assert read(path)


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")
@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_every_run_block_parses_as_shell(path: Path):
    for label, script in run_blocks(path):
        # Bytes, not text: in text mode Python translates "\n" to "\r\n" on the
        # way into stdin under Windows, and bash then rejects the carriage
        # returns it just received.
        result = subprocess.run(
            ["bash", "-n"],
            input=strip_expressions(script).encode("utf-8"),
            capture_output=True,
        )
        assert result.returncode == 0, f"{label}\n{result.stderr.decode('utf-8', 'replace')}"
