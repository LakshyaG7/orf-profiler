"""End-to-end test of the Snakemake workflow on a tiny dataset."""

import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).parent.parent
TINY = Path(__file__).parent / "data" / "tiny.fa"

pytestmark = pytest.mark.skipif(
    shutil.which("snakemake") is None, reason="snakemake is not installed"
)


def test_workflow_runs_end_to_end(tmp_path):
    # Split the test FASTA into a protein-coding and a lncRNA file named like GENCODE's.
    records = TINY.read_text().split(">")[1:]
    data = tmp_path / "data"
    data.mkdir()
    groups = {
        "pc": [r for r in records if "CDS:" in r],
        "lncRNA": [r for r in records if "CDS:" not in r],
    }
    for dataset, group in groups.items():
        with gzip.open(data / f"gencode.vtest.{dataset}_transcripts.fa.gz", "wt") as handle:
            handle.write("".join(">" + record for record in group))

    results = tmp_path / "results"
    command = [
        "snakemake",
        "--snakefile", str(REPO / "workflow" / "Snakefile"),
        "--directory", str(REPO),
        "--cores", "1",
        "--config", "gencode_release=test", f"data_dir={data}", f"results_dir={results}",
    ]
    # Make sure the workflow calls the orf-profiler installed in this environment.
    env_path = str(Path(sys.executable).parent)
    completed = subprocess.run(
        command, capture_output=True, text=True,
        env={**os.environ, "PATH": env_path + os.pathsep + os.environ["PATH"]},
    )
    assert completed.returncode == 0, completed.stderr

    assert (results / "pc" / "orfs.tsv.gz").exists()
    assert (results / "lncRNA" / "orfs.tsv.gz").exists()
    assert (results / "figures" / "orf_length_by_class.png").stat().st_size > 0
    assert (results / "figures" / "kozak_by_class.png").stat().st_size > 0
