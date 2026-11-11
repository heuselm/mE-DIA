# mE-DIA from source and on the command line

This page is for running mE-DIA from a source checkout: on Linux, in WSL, or with
your own Python on Windows. If you only want the GUI on Windows, the
[README](../README.md#installation) route (`mEDIA.exe` + DIA-NN) needs none of this.

- [Setup](#setup)
- [Launching the GUI from source](#launching-the-gui-from-source)
- [Running from WSL](#running-from-wsl)
- [Running the pipeline from the command line](#running-the-pipeline-from-the-command-line)
- [Modifications](#modifications)
- [Enzymes and missed cleavages](#enzymes-and-missed-cleavages)
- [All options](#all-options)
- [Output](#output)
- [Tests](#tests)

---

## Setup

Requirements: **Python 3.10+**, **DIA-NN 1.8.1**, and a FASTA protein database. Node.js is
needed only if you rebuild the GUI (`bin/ui/`), and the built GUI ships with the repository.

```bash
sudo apt install python3-venv          # Debian/Ubuntu/WSL, once; without it venv fails ("ensurepip is not available")
git clone https://github.com/heuselm/mE-DIA.git
cd mE-DIA
python3 -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
python -m pip install -r bin/requirements.txt
```

`bin/requirements.txt` installs `fastapi`, `uvicorn[standard]`, `websockets`, `pandas` and
`matplotlib`.

Use a virtual environment as above. On Debian, Ubuntu and WSL the system Python is
"externally managed" (PEP 668) and refuses `pip install` with
`error: externally-managed-environment`. A conda/mamba environment works too. Always install
with `python -m pip` from the environment you will run mE-DIA with. A bare `pip` on `PATH`
may belong to another interpreter.

DIA-NN is found, in order, from `--diann`, the `MEDIA_DIANN_EXE` environment variable,
and `diann` on `PATH`. In the GUI, set it under *DIA-NN Settings*.

## Launching the GUI from source

```bash
python bin/start.py                    # installs missing deps, starts the server, opens the browser
```

The GUI is served at **http://localhost:8765/**.

| `start.py` flag | Effect |
|---|---|
| `--port 8765` | Listening port |
| `--no-browser` | Do not open a browser |
| `--build-ui` | Rebuild the React UI from `bin/ui/` (needs Node.js 18+) |

`python bin/server.py --port 8765` starts the server directly, without the dependency check.
Settings saved in the GUI go to `bin/config.json`. To start with your DIA-NN path filled in,
`cp bin/config.example.json bin/config.json` and edit it.

## Running from WSL

The server and the pipeline can run in WSL2 (Ubuntu) while DIA-NN stays the Windows
`diann.exe`.

- **Start** with `python3 bin/start.py` from a WSL shell. The GUI opens in the Windows
  browser.
- **Keep data on a Windows drive.** Raw files, FASTA, libraries and the output folder must be
  under `/mnt/<drive>/`, because a Windows program cannot read the Linux filesystem. The run
  stops before calling DIA-NN and names the offending path.
- **DIA-NN path:** give `/mnt/c/DIA-NN/1.8.1/diann.exe` or `C:\DIA-NN\1.8.1\diann.exe`.
- **Paths** may be pasted in Windows form (`E:\data\a.raw`). They are converted to POSIX on
  input, and `pipeline.settings` records POSIX paths. Paths are converted back to Windows form
  only for `diann.exe` and `dia-nn-plotter.exe`.
- **File dialogs** open on Windows and return `/mnt/...` paths. Without Windows interop
  (`powershell.exe` unreachable), type or paste paths instead.
- **Rebuilding the UI** (`--build-ui`): if `bin/ui/node_modules` was installed by Windows npm,
  it lacks the Linux build tools. The launcher then builds with Windows npm, which needs the
  repository on a Windows drive.

## Running the pipeline from the command line

`bin/run_pipeline.py` is the pipeline. The GUI is a front end for exactly this script, and
`python bin/run_pipeline.py --help` prints the full reference.

Minimal run (four enzymes, fixed carbamidomethyl-Cys by default):

```bash
python bin/run_pipeline.py \
  --raw sample_TCGE_DIA.raw \
  --fasta proteome.fasta \
  --out-dir results/ \
  --enzymes TCGE \
  --threads 14 \
  --diann /path/to/diann
```

Single-enzyme data work the same way, e.g. `--enzymes T`.

With the semi-specific extension (passes 3 and 4) and initiator-methionine excision:

```bash
python bin/run_pipeline.py \
  --raw run1.raw run2.raw \
  --fasta proteome.fasta \
  --out-dir results/ \
  --enzymes TCGE --semi-ext --met-excision \
  --diann "C:/DIA-NN/1.8.1/diann.exe"
```

Reusing previously predicted libraries skips the slowest step:

```bash
python bin/run_pipeline.py \
  --raw sample.raw --fasta proteome.fasta --out-dir results/ --enzymes TCGE \
  --speclib T=libs/T.predicted.speclib --speclib C=libs/C.predicted.speclib \
  --speclib G=libs/G.predicted.speclib --speclib E=libs/E.predicted.speclib
# or let the pipeline find *.predicted.speclib files in one folder:
#   --speclib-dir libs/
```

A library is reused only if its enzyme definition (cleavage rule, missed cleavages) matches.
Otherwise the run refuses rather than mixing definitions.

Skipping pass 1 with a ready-made combined library:

```bash
python bin/run_pipeline.py \
  --raw sample.raw --fasta proteome.fasta --out-dir results/ --enzymes TCGE \
  --skip-pass1 --lib tcge_combined.tsv
```

**Resuming.** Rerunning into an existing output folder skips completed steps. A step that was
interrupted (cancel, crash) is detected via `pipeline.inflight` and repeated. A resume with
different enzyme definitions is refused, so start a new output folder instead.

## Modifications

Modifications apply identically to every DIA-NN call and are recorded in `pipeline.settings`.

| Option | Effect |
|---|---|
| `--mods LIST` | Comma-separated presets, or `none`. Presets: `carbamidomethyl` (fixed C, **default**), `oxidation` (M), `acetyl-nterm` (protein N-term), `phospho` (STY), `gg` (K). |
| `--fixed-mod NAME,MASS,SITES` | Additional fixed modification, DIA-NN syntax, e.g. `UniMod:4,57.021464,C` |
| `--var-mod NAME,MASS,SITES` | Additional variable modification. Sites may include `n` (peptide N-term) and `*n` (protein N-term). |
| `--var-mods N` | Maximum variable modifications per peptide (default 1) |
| `--monitor-mod NAME` | Score this variable modification as a PTM of interest (DIA-NN `--monitor-mod`; adds `PTM.Q.Value` / `PTM.Site.Confidence`). The `acetyl-nterm`, `phospho` and `gg` presets set it. |
| `--no-cut-after-mod NAME` | The modified residue blocks cleavage (e.g. GlyGly-K against trypsin) |
| `--met-excision` | Variable protein N-terminal methionine excision |

Each variable modification multiplies the library size. Two missed cleavages plus two variable
modifications on a human proteome can exceed the memory of a 32 GB workstation during library
prediction.

## Enzymes and missed cleavages

| Code | Enzyme | Cleaves at | Blocked by | Sense | Missed cleavages |
|------|--------|-----------|------------|-------|------------------|
| `T`  | Trypsin | K, R | P | C-term | global |
| `T2` | Trypsin/P | K, R | – | C-term | global |
| `C`  | Chymotrypsin | F, W, Y, L | P | C-term | global |
| `G`  | GluC | D, E | P | C-term | global |
| `E`  | Elastase | A, V | – | C-term | global |
| `EX` | Elastase, extended (Promega V1891) | A, V, S, G, L, I | – | C-term | **3** |
| `L`  | LysC | K | P | C-term | global |
| `LN` | LysN | K | – | N-term | global |
| `R`  | ArgC | R | P | C-term | global |
| `A`  | AspN | D | – | N-term | global |
| `M`  | CNBr | M | – | C-term | global |
| `F`  | PepsinA | F, L | P | C-term | global |

Each enzyme runs with, in order of precedence:

1. `--enzyme-mc CODE=N` (repeatable; in the GUI, the field next to the enzyme);
2. the enzyme's own default (`EX`: 3);
3. the global `--missed-cleavages` (default 1).

The values apply to pass-1 prediction and to the pass-3 digest. They are recorded per enzyme in
`pipeline.settings` (`enzymeDefinitions`).

More enzymes can be added at runtime with `--enzyme-config myenzymes.json`: a JSON map of
`code → {code, name, cut, sense, cut_res, noclv, supports_semi, missed_cleavages}`, where the
last two are optional. See [`cleavage_rules.md`](cleavage_rules.md).

## All options

| Option | Default | Description |
|--------|---------|-------------|
| `--raw FILE [FILE ...]` | *required* | DIA input files (`.raw`, `.mzML`). Space-separated or repeated. |
| `--out-dir DIR` | *required* | Output directory |
| `--fasta FILE` | – | Protein database; required for pass-1 prediction and annotation |
| `--enzymes CODES` | `TCGE` | Enzyme codes, e.g. `TCGE` or `T,C,G,E`. One code for single-enzyme data. |
| `--threads N` | `14` | CPU threads passed to DIA-NN |
| `--semi-ext` | off | Semi-specific extension (pass 3 + pass 4) |
| `--speclib CODE=PATH` | – | Per-enzyme predicted speclib; repeat per enzyme |
| `--speclib-dir DIR` | – | Fallback: scan a directory for `*.predicted.speclib` |
| `--lib FILE` | – | Combined library TSV for pass 2 (skips pass 1) |
| `--skip-pass1` | off | Skip pass 1 (use with `--lib`) |
| `--diann PATH` | `$MEDIA_DIANN_EXE`, else `diann` | DIA-NN executable |
| `--min-pep-len N` | `7` | Minimum peptide length |
| `--max-pep-len N` | `30` | Maximum peptide length |
| `--missed-cleavages N` | `1` | Missed cleavages for enzymes without their own value |
| `--enzyme-mc CODE=N` | – | Missed cleavages for one enzyme; repeatable |
| `--enzyme-config FILE` | – | JSON file extending the enzyme registry |
| `--mods LIST` | `carbamidomethyl` | Modification presets, or `none` (see [Modifications](#modifications)) |
| `--fixed-mod` / `--var-mod` | – | Additional modifications, DIA-NN syntax; repeatable |
| `--var-mods N` | `1` | Maximum variable modifications per peptide |
| `--monitor-mod NAME` | – | Score a variable modification as PTM of interest; repeatable |
| `--no-cut-after-mod NAME` | – | Modified residue blocks cleavage; repeatable |
| `--met-excision` | off | Protein N-terminal Met excision |
| `--diann-arg=ARG` | – | Raw argument appended to every DIA-NN call; repeatable |
| `--no-isd-filter` | filter on | Proteoform step: keep in-source-decay products |
| `--isd-rt-window MIN` | `0.05` | Proteoform step: max apex RT difference (min) to a longer, more intense parent for an in-source-decay flag |

## Output

| Folder | Contents |
|--------|----------|
| `pass1_<E>/` | Per-enzyme search: `report.tsv`, predicted and empirical libraries |
| `pass2_<ENZYMES>/` | Combined fully specific search and quantification matrices |
| `pass3_semiExt_<E>/` | Per-enzyme semi-specific search (only with `--semi-ext`) |
| `pass4_semiExt_<ENZYMES>/` | Final semi-extended search (only with `--semi-ext`) |

The **final pass folder** (`pass4_semiExt_<ENZYMES>/` with `--semi-ext`, otherwise
`pass2_<ENZYMES>/`) contains:

| File | Contents |
|---|---|
| `report_annotated.tsv` | DIA-NN report with `Enzyme` and `Specificity` columns per precursor |
| `report.pr_matrix.tsv`, `report.pg_matrix.tsv`, `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv` | DIA-NN quantification matrices |
| `report.used.fasta` | The protein sequences the search used, for the Peptide Mapper |
| `PeptideMapper.html` | Interactive Peptide Mapper; open it in a browser |
| `report.pf_matrix.tsv` | Candidate proteoform events, with `isd_flag` / `isd_parent` for in-source-decay products |
| `report.pf_isd.tsv` | All in-source-decay flags |
| `report.pf_log.txt`, `report.pf_plot.pdf` | Proteoform step log and summary plot |

At the top level of the output directory:

- `pipeline.log`: DIA-NN output for every pass.
- `pipeline.settings`: the exact settings, reloadable in the GUI. It also records `buildInfo`
  (mE-DIA commit) and `runtime` (Python and library versions, DIA-NN version banner).
- `pipeline.inflight`: present only while a DIA-NN call runs or after one was interrupted.

DIA-NN also writes `.quant` files next to each raw file and overwrites them on every search of
that file.

## Tests

Tests are plain scripts. Each `bin/test_*.py` prints PASS/FAIL per check and exits non-zero on
failure:

```bash
for t in bin/test_*.py; do python "$t" > /dev/null && echo "ok   $t" || echo "FAIL $t"; done
```
