# Maintaining the KEGG data artefacts

This guide is for the **package maintainer** who rebuilds raven-toolbox's KEGG
artefacts once per KEGG release. End users never do this — they download the
published, version-pinned artefacts. The build has two implemented steps so far:
**3b.1 download** (`reconstruction/kegg/download.py`) and **3b.2 parse**
(`reconstruction/kegg/parse.py`); see PLAN.md §2.3b for the full pipeline.

## Prerequisites

### A paid KEGG FTP subscription
The bulk KEGG dump is licensed. You need an active subscription to
`ftp.kegg.net`, which gives you a **username and password**.

### Credentials in `~/.netrc`
The download reads your KEGG username and password from a `~/.netrc` file — it
never takes them on the command line, so they stay out of your shell history and
out of `ps` output. Create the file (readable only by you) and add a `machine`
line for the KEGG host:

```bash
touch ~/.netrc && chmod 600 ~/.netrc
```

Then add this single line to `~/.netrc`, substituting your subscription
credentials:

```
machine ftp.kegg.net login YOUR_KEGG_USER password YOUR_KEGG_PASSWORD
```

Notes:
- The host **must be `ftp.kegg.net`** — that is the machine name the downloader
  looks up. A `machine` line for any other host is ignored.
- The file **must be mode `600`** (owner read/write only). Python's `netrc`
  parser refuses a `.netrc` that other users can read.
- `~/.netrc` is the same convention `curl`, `wget` and `git` use, so if you
  already have one, just add the `ftp.kegg.net` line to it.

If you keep secrets somewhere other than `$HOME`, point the downloader at a
different file with `netrc_path=...` (see below); the format is identical.

## Step 3b.1 — download and arrange the dump

With `~/.netrc` in place, no credentials need to be passed in code:

```python
from raven_toolbox.reconstruction.kegg import download_kegg_dump

# Reads ~/.netrc, fetches the KEGG archives, extracts and arranges them.
download_kegg_dump("keggdb")
```

This fetches the reaction / compound / glycan / ko archives, the eukaryote and
prokaryote proteomes, and the taxonomy file; extracts them; and arranges the
flat layout the parser expects (`reaction`, `reaction.lst`,
`reaction_mapformula.lst`, `compound` = compound + glycan, `compound.inchi`,
`ko`, `genes.pep` = both proteomes, `taxonomy`).

Credential alternatives:

```python
# A .netrc in a non-default location:
download_kegg_dump("keggdb", netrc_path="/run/secrets/kegg_netrc")

# Pass credentials explicitly (only when they come from a secret manager at
# runtime — never hardcode literals in committed code):
download_kegg_dump("keggdb", auth=("YOUR_KEGG_USER", "YOUR_KEGG_PASSWORD"))
```

Already-downloaded files are skipped; pass `force=True` to re-fetch (for a new
KEGG release).

Per-file download and extraction progress bars are shown by default; pass
`progress=False` for non-interactive runs (e.g. logging to a file).

## Step 3b.2 — parse into the published artefacts

```python
from raven_toolbox.reconstruction.kegg import parse_kegg_dump

parse_kegg_dump("keggdb", "artefacts", version="kegg116")
```

This writes the gene-free reference model (`reference_model.yml.gz`, gzipped
RAVEN/cobra YAML) and the relational tables as gzipped TSV. With `version=` set,
every output filename is version-prefixed (e.g. `kegg116_organism_gene_ko.tsv.gz`),
matching the published release assets. See
[kegg_data_format.md](kegg_data_format.md) for what those tables contain and the
format rationale.

Pass `progress=True` to report each parse stage and show a progress bar over the
large `organism_gene_ko` (ko) streaming pass.

## Step 3b.3 — build the HMM libraries

Build the per-domain profile-HMM libraries that the de-novo query path (3b.5)
searches. This needs **HMMER** (`hmmbuild`, `hmmpress`), **MAFFT**, and
**CD-HIT** on `PATH` (or set `RAVEN_PYTHON_HMMBUILD` / `RAVEN_PYTHON_MAFFT` /
`RAVEN_PYTHON_CDHIT`, etc.); install e.g. `conda install -c bioconda hmmer mafft cd-hit`.

> **OS note:** these three tools run on Linux and macOS but **not native
> Windows** — on Windows, run this step inside WSL2. See the native-OS-support
> matrix in [maintaining_binaries.md](maintaining_binaries.md#native-os-support-per-tool).

```python
from raven_toolbox.reconstruction.kegg import build_hmm_library, read_kegg_table

organism_gene_ko = read_kegg_table("artefacts/kegg116_organism_gene_ko.tsv.gz")
for domain in ("prokaryotes", "eukaryotes"):
    build_hmm_library(
        organism_gene_ko,
        "keggdb/genes.pep",      # proteomes from 3b.1
        "keggdb/taxonomy",       # domain split, from 3b.1
        f"hmms/{domain}",
        domain=domain,
        progress=True,           # show an "N of M KOs" bar for the long build
    )
```

For each KO in the domain it gathers the member sequences, dereplicates with
CD-HIT (~90 % identity), aligns with MAFFT, trains a profile with `hmmbuild`, and
finally concatenates them into a single `library.hmm`. This is the slowest step
(hours, once per KEGG release); it skips KOs whose `.hmm` already exists, so it is
resumable. The concatenated library is published **gzipped** as
`<version>_<domain>.hmm.gz` (e.g. `kegg116_prokaryotes.hmm.gz`); end users
decompress it and search it directly with `hmmsearch` (see `ensure_kegg_hmm_library`),
which keeps the download ~10× smaller than a binary index, stays portable across
HMMER versions, and lets the same artefact serve MATLAB RAVEN.

## Building and publishing in one go

[`scripts/build_kegg_artefacts.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/README.md) runs 3b.2 (+ 3b.3 with
`--hmms`) and lays the output out as publishable, version-prefixed assets: the core
model files (reference model + KO/reaction tables) bundled into `<version>_core.tar.gz`
(which `ensure_kegg_data` extracts), and `<version>_<domain>.hmm.gz` per domain (named
for `ensure_kegg_hmm_library`). It also publishes
`<version>_taxonomy.gz` — the domain split plus the source for
[`phyl_dist`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/taxonomy.py),
which regenerates RAVEN's `keggPhylDist` (used by GECKO) with no `.mat` file:

```bash
python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts \
    --version kegg116 --hmms --threads 8
```

The build is **idempotent**: if it fails partway (the HMM step can run for hours),
re-run the *same command* — each stage is skipped when its output already exists
(parsed tables, taxonomy, per-domain HMM library, core bundle) and the per-KO HMM
build resumes where it left off, so finished work is not repeated. Pass `--force`
to rebuild everything from scratch.

Upload the contents of `artefacts/` to the release, then record the artefacts in
both the shared `data/manifest.json` and `raven_toolbox.data._DATA_REGISTRY` with
[`scripts/make_registry_snippet.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/README.md)
(it computes each file's SHA256 + size):

```bash
# shared source of truth (read by raven-toolbox and MATLAB RAVEN):
# --tag is the release tag; the script builds the …/releases/download/<tag> URL itself.
python scripts/make_registry_snippet.py manifest --manifest data/manifest.json \
    --target data --dataset kegg --version kegg116 --dir artefacts --tag v0.1.0
# in-code registry, so end users auto-fetch with no env var (paste into _DATA_REGISTRY):
python scripts/make_registry_snippet.py data --dataset kegg --version kegg116 \
    --dir artefacts --tag v0.1.0
```

From then on `ensure_data` fetches and verifies the artefacts for end users
automatically.

## End-user paths (3b.4 / 3b.5)

End users do **not** run the steps above — the published artefacts are fetched and
cached automatically by `ensure_data` (`raven_toolbox.data`) under
`~/.cache/raven-toolbox/data/kegg-<version>/` on first use, so the entry points below
can be called with no local paths at all (pass an explicit `artefact_dir=`/
`library=` to use your own build instead). Two runtime entry points build a draft
model from the artefacts:

- **3b.4 — species in KEGG** (`get_kegg_model_for_organism_from_artefacts`): no
  binaries needed; uses the organism's KEGG gene↔KO annotations. Fully
  cross-platform. `organism_id="prokaryotes"`/`"eukaryotes"` builds a whole-domain
  model (pass `taxonomy=`).
- **3b.5 — organism not in KEGG** (`get_kegg_model_from_sequences`): `hmmscan`-es a
  proteome FASTA against the domain library (decompressed and `hmmpress`-ed from
  the published `.hmm.gz` on first use), so it needs **HMMER** (`hmmpress`,
  `hmmscan`) — Linux/macOS or WSL2 (see the OS matrix). Tune assignment with
  `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`.
