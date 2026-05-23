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

This writes the gene-free reference model (`reference_model.xml`) and the
relational tables as gzipped TSV. See [kegg_data_format.md](kegg_data_format.md)
for what those tables contain and why they use gzipped TSV.
