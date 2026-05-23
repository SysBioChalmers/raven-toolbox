# Maintaining the KEGG data artefacts

This guide is for the **package maintainer** who rebuilds ravengem's KEGG
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
from ravengem.reconstruction.kegg import download_kegg_dump

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

## Step 3b.2 — parse into the published artefacts

```python
from ravengem.reconstruction.kegg import parse_kegg_dump

parse_kegg_dump("keggdb", "artefacts")
```

This writes the gene-free reference model (`reference_model.yml.gz`, gzipped
RAVEN/cobra YAML) and the relational tables as gzipped TSV. See
[kegg_data_format.md](kegg_data_format.md) for what those tables contain and the
format rationale.

## Step 3b.3 — build the HMM libraries

Build the per-domain profile-HMM libraries that the de-novo query path (3b.5)
searches. This needs **HMMER** (`hmmbuild`, `hmmpress`), **MAFFT**, and
**CD-HIT** on `PATH` (or set `RAVENGEM_HMMBUILD` / `RAVENGEM_MAFFT` /
`RAVENGEM_CDHIT`, etc.); install e.g. `conda install -c bioconda hmmer mafft cd-hit`.

> **OS note:** these three tools run on Linux and macOS but **not native
> Windows** — on Windows, run this step inside WSL2. See the native-OS-support
> matrix in [maintaining_binaries.md](maintaining_binaries.md#native-os-support-per-tool).

```python
from ravengem.reconstruction.kegg import build_hmm_library, read_kegg_table

organism_gene_ko = read_kegg_table("artefacts/organism_gene_ko.tsv.gz")
for domain in ("prokaryotes", "eukaryotes"):
    build_hmm_library(
        organism_gene_ko,
        "keggdb/genes.pep",      # proteomes from 3b.1
        "keggdb/taxonomy",       # domain split, from 3b.1
        f"hmms/{domain}",
        domain=domain,
    )
```

For each KO in the domain it gathers the member sequences, dereplicates with
CD-HIT (~90 % identity), aligns with MAFFT, trains a profile with `hmmbuild`, and
finally concatenates and `hmmpress`-es them into a single `library.hmm` for fast
`hmmscan` querying. This is the slowest step (hours, once per KEGG release); it
skips KOs whose `.hmm` already exists, so it is resumable. The resulting
libraries are published as version-pinned artefacts alongside the reference model
and tables.

## Building and publishing in one go

[`scripts/build_kegg_artefacts.py`](../scripts/README.md) runs 3b.2 (+ 3b.3 with
`--hmms`) and lays the output out as publishable assets (`<domain>.hmm` named for
`ensure_kegg_hmm_library`):

```bash
python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts --hmms --threads 8
```

Upload the contents of `artefacts/` to a release, then emit the registry entry for
`ravengem.data._DATA_REGISTRY` with [`scripts/make_registry_snippet.py`](../scripts/README.md):

```bash
python scripts/make_registry_snippet.py data --dataset kegg --version kegg116 \
    --dir artefacts --base-url https://github.com/ORG/ravengem/releases/download/kegg-data-kegg116
```

Paste the printed block into `_DATA_REGISTRY`; from then on `ensure_data` fetches
and verifies the artefacts for end users automatically.

## End-user paths (3b.4 / 3b.5)

End users do **not** run the steps above — the published artefacts are fetched and
cached automatically by `ensure_data` (`ravengem.data`) under
`~/.cache/ravengem/data/kegg-<version>/` on first use, so the entry points below
can be called with no local paths at all (pass an explicit `artefact_dir=`/
`library=` to use your own build instead). Two runtime entry points build a draft
model from the artefacts:

- **3b.4 — species in KEGG** (`get_kegg_model_for_organism_from_artefacts`): no
  binaries needed; uses the organism's KEGG gene↔KO annotations. Fully
  cross-platform. `organism_id="prokaryotes"`/`"eukaryotes"` builds a whole-domain
  model (pass `taxonomy=`).
- **3b.5 — organism not in KEGG** (`get_kegg_model_from_sequences`): `hmmscan`-es a
  proteome FASTA against the pressed `library.hmm`, so it needs **HMMER**
  (`hmmscan`) — Linux/macOS or WSL2 (see the OS matrix). Tune assignment with
  `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`.
