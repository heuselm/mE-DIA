#!/usr/bin/env python3
"""
Multi-enzyme DIA Proteomics Pipeline
=====================================
Unified pipeline for DIA-NN-based multi-enzyme proteomics analysis.

Enzyme registry matches cleavage_rules.md (sense 0 = N-terminal, 1 = C-terminal).

Usage:
  python bin/run_pipeline.py \\
    --raw sample1.raw [sample2.raw ...] \\
    --fasta database.fasta \\
    --out-dir results/ \\
    --enzymes TCGE \\
    --threads 14 \\
    [--semi-ext] \\
    [--speclib T=/path/T.predicted.speclib] \\
    [--speclib C=/path/C.predicted.speclib] \\
    [--speclib-dir fallback/scan/dir] \\
    [--lib combined.tsv]

Output structure:
  results/
    pass1_T/               per-enzyme predict + search
    pass1_C/
    pass2_TCGE/            combined-library search + matrices
    pass3_semiExt_T/       (if --semi-ext) inline-reannotated semi search
    pass4_semiExt_TCGE/    (if --semi-ext) combined semi-ext search + matrices
"""

import argparse
import dataclasses
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

# Local module: build provenance, so every result folder records the code that
# produced it (bin/build_info.py)
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import build_info as _build_info  # noqa: E402
    BUILD_INFO = _build_info.get()
    RUNTIME_INFO = _build_info.runtime
except Exception:
    _build_info = None
    BUILD_INFO = {"commit": "unknown", "branch": "unknown", "source": "unavailable"}

    def RUNTIME_INFO(diann_exe=None):
        return {"python": "unknown", "packages": {}, "diann": {"exe": str(diann_exe or "unknown"),
                                                              "banner": "unknown"}}

# Local module: WSL detection and the path boundary to Windows executables
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wslpaths  # noqa: E402

# Local module: end-of-pipeline proteoform analysis (sibling of this file)
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import proteoform_analysis as _pf  # noqa: E402
except Exception as _pf_err:
    _pf = None
    print(f"  NOTE: proteoform_analysis unavailable ({_pf_err}); will skip.")

# ─── Enzyme Registry ─────────────────────────────────────────────────────────
# Mirrors cleavage_rules.md:
#   sense=1  → C-terminal cutter (cleave AFTER cut_res, blocked by noclv as NEXT residue)
#   sense=0  → N-terminal cutter (cleave BEFORE cut_res)
# supports_semi=False means skip virtual-FASTA semi-extension for this enzyme
# (all C-terminal cutters support it; N-terminal cutters also supported now)

@dataclasses.dataclass
class EnzymeDef:
    code:          str          # short identifier used in folder names
    name:          str          # display name
    cut:           str          # DIA-NN --cut pattern
    sense:         int          # 1=C-terminal, 0=N-terminal
    cut_res:       str          # string of cleavage residues
    noclv:         str          # residues that block cleavage (next AA for sense=1)
    supports_semi: bool = True  # False for no-enzyme / special cases
    # The enzyme's own missed-cleavage allowance; None = the run's global
    # --missed-cleavages. Last, with a default, so positional construction and
    # --enzyme-config JSON written before it keep working.
    missed_cleavages: Optional[int] = None

    @property
    def stem(self) -> str:
        """Filename-safe form of the name, for the libraries named after it."""
        return file_stem(self.name)

    @property
    def cterm(self) -> str:
        """Residues that appear at the C-terminus of fully-specific peptides."""
        return self.cut_res if self.sense == 1 else ""

    @property
    def nterm(self) -> str:
        """Residues that appear at the N-terminus of fully-specific peptides."""
        return self.cut_res if self.sense == 0 else ""


def file_stem(name: str) -> str:
    """'Elastase (extended)' -> 'elastase_extended'; existing names map to their
    plain lower case ('Elastase' -> 'elastase'), except 'Trypsin/P' -> 'trypsin_p'."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


ENZYME_REGISTRY: Dict[str, EnzymeDef] = {
    # code  name            DIA-NN --cut                   s  cut_res  noclv
    "T":  EnzymeDef("T",  "Trypsin",       "K*,R*,!*P",       1, "KR",   "P"),
    "T2": EnzymeDef("T2", "Trypsin/P",     "K*,R*",            1, "KR",   ""),  # no pro rule
    "C":  EnzymeDef("C",  "Chymotrypsin",  "F*,W*,Y*,L*,!*P", 1, "FWYL", "P"),
    "G":  EnzymeDef("G",  "GluC",          "D*,E*,!*P",       1, "DE",   "P"),
    "E":  EnzymeDef("E",  "Elastase",      "A*,V*",            1, "AV",   ""),
    # Promega Elastase V1891: C-terminal to A, V, S, G, L, I (docs/cleavage_rules.md,
    # rule 11b). Six residues make far more internal sites than A/V, so it carries
    # its own missed-cleavage allowance (#35).
    "EX": EnzymeDef("EX", "Elastase (extended)", "A*,V*,S*,G*,L*,I*", 1, "AVSGLI", "",
                    missed_cleavages=3),
    "L":  EnzymeDef("L",  "LysC",          "K*,!*P",           1, "K",    "P"),
    "LN": EnzymeDef("LN", "LysN",          "*K",               0, "K",    ""),
    "R":  EnzymeDef("R",  "ArgC",          "R*,!*P",           1, "R",    "P"),
    "A":  EnzymeDef("A",  "AspN",          "*D",               0, "D",    ""),
    "M":  EnzymeDef("M",  "CNBr",          "M*",               1, "M",    ""),
    "F":  EnzymeDef("F",  "PepsinA",       "F*,L*,!*P",       1, "FL",   "P"),
}


# ─── Modifications ───────────────────────────────────────────────────────────
# DIA-NN 1.8.1 flags (README): --fixed-mod / --var-mod NAME,MASS,SITES[,label],
# --var-mods N (max variable mods per peptide), --met-excision. Sites are amino
# acid letters, 'n' (peptide N-term) or '*n' (protein N-term). --unimod4 is not
# a documented 1.8.1 option and is deliberately not used.
#
# The same configuration must reach every DIA-NN call: a library predicted under
# one configuration and searched under another is a silent mass mismatch.
# Masses: monoisotopic, from the UniMod elemental compositions.

MOD_PRESETS: Dict[str, dict] = {
    "carbamidomethyl": {"name": "UniMod:4",  "mass": 57.021464, "sites": "C",  "type": "fixed"},
    "oxidation":       {"name": "UniMod:35", "mass": 15.994915, "sites": "M",  "type": "variable"},
    # The DIA-NN 1.8.1 GUI adds --monitor-mod for a PTM of interest, which is what produces
    # the PTM.Site.Confidence / PTM.Q.Value columns; it does not do so for oxidation. The
    # presets below reproduce that, verified against GUI reference runs (see docs/plan.md).
    "acetyl-nterm":    {"name": "UniMod:1",  "mass": 42.010565, "sites": "*n", "type": "variable",
                        "monitor": True},
    "phospho":         {"name": "UniMod:21", "mass": 79.966331, "sites": "STY", "type": "variable",
                        "monitor": True},
    "gg":              {"name": "UniMod:121", "mass": 114.042927, "sites": "K", "type": "variable",
                        "monitor": True, "noCutAfter": True},
}
DEFAULT_MOD_PRESETS = ["carbamidomethyl"]
_MOD_RESIDUES = set("ACDEFGHIKLMNPQRSTUVWY")


class ModConfigError(ValueError):
    """Malformed modification, or a library built under another configuration."""


def _mod_sites(sites: str) -> str:
    """Validate a DIA-NN sites string: residues plus optional 'n' / '*n'."""
    rest = sites.strip()
    if not rest:
        raise ModConfigError("modification sites are empty")
    while rest:
        if rest.startswith("*n"):
            rest = rest[2:]
        elif rest[0] == "n" or rest[0] in _MOD_RESIDUES:
            rest = rest[1:]
        else:
            raise ModConfigError(
                f"unknown modification site '{rest[0]}' in '{sites}' "
                f"(allowed: amino-acid letters, 'n' = peptide N-term, '*n' = protein N-term)")
    return sites.strip()


def make_mod(name, mass, sites, mod_type, monitor=False, no_cut_after=False) -> dict:
    """Validated, normalised {name, mass, sites, type}; raises ModConfigError.

    monitor      -> DIA-NN --monitor-mod: score this PTM's site localisation.
    no_cut_after -> DIA-NN --no-cut-after-mod: the modified residue blocks cleavage.
    Both keys are omitted when false, so a record stays comparable to one written
    before they existed.
    """
    name = str(name).strip()
    if not name or any(c in name for c in ", \t"):
        raise ModConfigError(f"bad modification name '{name}' (no commas or spaces)")
    try:
        mass_f = float(mass)
    except (TypeError, ValueError):
        raise ModConfigError(f"bad mass '{mass}' for modification {name}: not a number")
    if not math.isfinite(mass_f) or mass_f == 0 or abs(mass_f) > 2000:
        raise ModConfigError(f"bad mass '{mass}' for modification {name}: "
                             f"must be a finite, non-zero delta mass below 2000 Da")
    if mod_type not in ("fixed", "variable"):
        raise ModConfigError(f"modification type must be 'fixed' or 'variable', not '{mod_type}'")
    mod = {"name": name, "mass": round(mass_f, 6), "sites": _mod_sites(str(sites)),
           "type": mod_type}
    if monitor:
        if mod_type != "variable":
            raise ModConfigError(
                f"only a variable modification can be monitored, and {name} is {mod_type}: "
                f"DIA-NN's --monitor-mod scores the localisation of a PTM it searches")
        mod["monitor"] = True
    if no_cut_after:
        mod["noCutAfter"] = True
    return mod


def parse_mod_spec(spec: str, mod_type: str) -> dict:
    """'NAME,MASS,SITES' (the DIA-NN --fixed-mod / --var-mod syntax) → mod dict."""
    parts = [p.strip() for p in str(spec).split(",")]
    if len(parts) != 3:
        raise ModConfigError(f"modification '{spec}' must be NAME,MASS,SITES "
                             f"(e.g. UniMod:4,57.021464,C)")
    return make_mod(parts[0], parts[1], parts[2], mod_type)


def make_mod_config(presets=(), fixed=(), variable=(), max_var_mods=1,
                    met_excision=False, extra=(), monitor=(), no_cut_after=()) -> dict:
    """The modification settings recorded in pipeline.settings and applied everywhere."""
    mods: List[dict] = []
    for p in presets:
        if p not in MOD_PRESETS:
            raise ModConfigError(f"unknown modification preset '{p}' "
                                 f"(known: {', '.join(MOD_PRESETS)})")
        m = MOD_PRESETS[p]
        mods.append(make_mod(m["name"], m["mass"], m["sites"], m["type"],
                             m.get("monitor", False), m.get("noCutAfter", False)))
    mods += [parse_mod_spec(s, "fixed") for s in fixed]
    mods += [parse_mod_spec(s, "variable") for s in variable]
    seen = set()
    for m in mods:
        key = (m["name"], m["sites"], m["type"])
        if key in seen:
            raise ModConfigError(f"modification {m['name']} on {m['sites']} given twice")
        seen.add(key)
    # --monitor-mod / --no-cut-after-mod name a modification that must already be declared;
    # a typo would otherwise silently drop PTM scoring, which is the failure being fixed.
    by_name = {m["name"]: m for m in mods}
    for names, flag in ((monitor, "--monitor-mod"), (no_cut_after, "--no-cut-after-mod")):
        for raw in names:
            nm = str(raw).strip()
            if not nm:
                continue
            if nm not in by_name:
                raise ModConfigError(
                    f"{flag} {nm}: no modification named {nm} is configured "
                    f"(configured: {', '.join(by_name) or 'none'})")
            if flag == "--monitor-mod":
                if by_name[nm]["type"] != "variable":
                    raise ModConfigError(
                        f"--monitor-mod {nm}: only a variable modification can be monitored")
                by_name[nm]["monitor"] = True
            else:
                by_name[nm]["noCutAfter"] = True
    if int(max_var_mods) < 0:
        raise ModConfigError("--var-mods must be >= 0")
    return {"modifications": mods, "maxVarMods": int(max_var_mods),
            "metExcision": bool(met_excision), "extraDiannArgs": list(extra)}


def _fmt_mass(m: float) -> str:
    return f"{m:.6f}".rstrip("0").rstrip(".")


def mod_flags(cfg: Optional[dict], protein_level: bool = True,
              include_extra: bool = True) -> List[str]:
    """DIA-NN flags for a modification config.

    protein_level=False is for the pass-3 virtual FASTA, where every entry is a
    single candidate peptide: protein-N-term ('*n') modifications and
    --met-excision would apply to every peptide there, so they are left out.
    """
    if not cfg:
        return []
    out: List[str] = []
    n_var = 0
    for m in cfg.get("modifications", []):
        sites = m["sites"]
        if not protein_level and "*n" in sites:
            sites = sites.replace("*n", "")
            if not sites:
                continue
        flag = "--fixed-mod" if m["type"] == "fixed" else "--var-mod"
        out += [flag, f"{m['name']},{_fmt_mass(m['mass'])},{sites}"]
        # emitted per modification, immediately after it, as the DIA-NN GUI does; skipped
        # together with a modification that was dropped above
        if m.get("monitor"):
            out += ["--monitor-mod", m["name"]]
        if m.get("noCutAfter"):
            out += ["--no-cut-after-mod", m["name"]]
        n_var += m["type"] == "variable"
    if n_var:
        out += ["--var-mods", str(cfg.get("maxVarMods", 1))]
    if protein_level and cfg.get("metExcision"):
        out.append("--met-excision")
    if include_extra:
        out += list(cfg.get("extraDiannArgs", []))
    return out


def mod_config_to_cli(cfg: dict) -> List[str]:
    """run_pipeline.py arguments that reproduce cfg exactly (used by the GUI server)."""
    args = ["--mods", "none"]
    for m in cfg.get("modifications", []):
        flag = "--fixed-mod" if m["type"] == "fixed" else "--var-mod"
        args += [flag, f"{m['name']},{_fmt_mass(m['mass'])},{m['sites']}"]
        if m.get("monitor"):
            args += ["--monitor-mod", m["name"]]
        if m.get("noCutAfter"):
            args += ["--no-cut-after-mod", m["name"]]
    args += ["--var-mods", str(cfg.get("maxVarMods", 1))]
    if cfg.get("metExcision"):
        args.append("--met-excision")
    for a in cfg.get("extraDiannArgs", []):
        args.append(f"--diann-arg={a}")
    return args


def describe_mods(flags: Optional[List[str]]) -> str:
    return " ".join(flags) if flags else "(no modifications)"


def settings_for_lib(lib_path: str) -> Optional[str]:
    """pipeline.settings of the mE-DIA run that built lib_path, if there is one.

    The pipeline writes every library one folder below its output folder
    (pass1_T/, pass2_TCGE/, pass3_semiExt_T/ ...), so the run's settings are in
    the library's parent folder; the library's own folder is checked too. This
    replaces a per-library .mods.json record: the configuration is already in
    pipeline.settings, and the exact flags DIA-NN received are in its log.
    """
    here = os.path.dirname(os.path.abspath(lib_path))
    for d in (here, os.path.dirname(here)):
        cand = os.path.join(d, "pipeline.settings")
        if os.path.exists(cand):
            return cand
    return None


def check_lib_mods(lib_path: str, expected: List[str], label: str,
                   protein_level: bool = True,
                   enzymes: Optional[List[dict]] = None) -> None:
    """Refuse a library built under another modification configuration.

    The configuration is recovered from the pipeline.settings of the run that
    built the library, and turned into flags exactly as that pass would have
    (protein_level=False for pass 3). A pipeline.settings without modification
    keys comes from before modification support, which applied none. Without a
    pipeline.settings (a library from elsewhere) the configuration cannot be
    recovered: warn loudly, do not assume.
    """
    rec = settings_for_lib(lib_path)
    if rec is None:
        print(f"\n  !! WARNING: {label}: no pipeline.settings found for "
              f"{os.path.basename(lib_path)} (looked in its folder and the one above).\n"
              f"  !! Cannot verify it was built with: {describe_mods(expected)}\n"
              f"  !! If it was built under another configuration the search is "
              f"silently wrong.\n", flush=True)
        return
    with open(rec, encoding="utf-8") as f:
        old = json.load(f)
    built = mod_flags(old if "modifications" in old else None,
                      protein_level=protein_level, include_extra=False)
    if built != expected:
        raise ModConfigError(
            f"{label}: {lib_path} was built with modifications [{describe_mods(built)}] "
            f"but this run uses [{describe_mods(expected)}]. Rebuild the library or "
            f"run with the matching configuration.")
    if enzymes is not None:
        made = _lib_built_specs(lib_path, old)
        if sorted(map(_rule, made), key=str) != sorted(map(_rule, enzymes), key=str):
            raise ModConfigError(
                f"{label}: {lib_path} was built with enzyme definitions "
                f"[{'; '.join(map(describe_enzyme, made)) or 'unknown'}] but this run "
                f"uses [{'; '.join(map(describe_enzyme, enzymes))}]. Rebuild the "
                f"library or run with the matching enzymes and missed cleavages.")


def mod_config_from_args(args) -> dict:
    """Modification config from parsed run_pipeline.py arguments."""
    presets = [p.strip().lower() for p in str(args.mods).split(",") if p.strip()]
    if presets == ["none"]:
        presets = []
    return make_mod_config(presets, args.fixed_mod, args.var_mod, args.var_mods,
                           args.met_excision, args.diann_arg,
                           args.monitor_mod, args.no_cut_after_mod)


_PASS_DIR_RE = re.compile(r"^pass[1-4]_")


# ─── Enzyme definitions: effective missed cleavages, provenance, guards ──────

def parse_enzyme_mc(items: List[str]) -> Dict[str, int]:
    """--enzyme-mc CODE=N (repeatable) -> {CODE: N}. ValueError on a bad entry."""
    out: Dict[str, int] = {}
    for kv in items or []:
        code, sep, n = str(kv).partition("=")
        code = code.strip().upper()
        if not sep or not code:
            raise ValueError(f"--enzyme-mc {kv!r}: expected CODE=N, e.g. EX=3")
        if code not in ENZYME_REGISTRY:
            raise ValueError(f"--enzyme-mc {kv!r}: unknown enzyme code '{code}' "
                             f"(known: {', '.join(ENZYME_REGISTRY)})")
        try:
            val = int(n.strip())
        except ValueError:
            raise ValueError(f"--enzyme-mc {kv!r}: N must be a whole number")
        if val < 0:
            raise ValueError(f"--enzyme-mc {kv!r}: N must be 0 or more")
        out[code] = val
    return out


def resolve_enzyme_mc(codes: List[str], global_mc: int,
                      overrides: Optional[Dict[str, int]] = None) -> Dict[str, int]:
    """Effective missed cleavages per enzyme: the --enzyme-mc override, else the
    enzyme's own default, else the run's global --missed-cleavages."""
    overrides = overrides or {}
    out = {}
    for c in codes:
        own = ENZYME_REGISTRY[c].missed_cleavages
        out[c] = overrides.get(c, own if own is not None else global_mc)
    return out


def enzyme_spec(code: str, missed_cleavages: int) -> dict:
    """The definition one enzyme ran with, as recorded in pipeline.settings."""
    e = ENZYME_REGISTRY[code]
    return {"code": code, "name": e.name, "cut": e.cut, "sense": e.sense,
            "cutResidues": e.cut_res, "blocking": e.noclv,
            "missedCleavages": int(missed_cleavages)}


def _rule(spec: dict) -> tuple:
    """What determines the search space: everything but the code and the name."""
    return (spec.get("cut"), spec.get("sense"), spec.get("cutResidues"),
            spec.get("blocking"), spec.get("missedCleavages"))


def describe_enzyme(spec: dict) -> str:
    return (f"{spec.get('code')} {spec.get('name')} [cut={spec.get('cut')}, "
            f"sense={spec.get('sense')}, residues={spec.get('cutResidues')}, "
            f"blocking={spec.get('blocking') or '-'}, "
            f"missed cleavages={spec.get('missedCleavages')}]")


def recorded_enzyme_specs(settings: dict) -> Dict[str, dict]:
    """{code: definition} from a pipeline.settings dict.

    Settings written before enzymeDefinitions existed recorded only codes and one
    global missedCleavages. No existing code's definition has changed since, so
    the current registry definition with that value is exact for them.
    """
    if isinstance(settings.get("enzymeDefinitions"), list):
        return {d.get("code"): d for d in settings["enzymeDefinitions"]}
    mc = settings.get("missedCleavages", 1)
    return {c: enzyme_spec(c, mc) for c in settings.get("enzymes", [])
            if c in ENZYME_REGISTRY}


def check_resume_enzymes(out_dir: str, specs: List[dict]) -> None:
    """Refuse to resume into out_dir if any enzyme of this run was not recorded
    there with the same definition (rule and effective missed cleavages)."""
    settings = os.path.join(out_dir, "pipeline.settings")
    if not os.path.exists(settings) or not any(
            _PASS_DIR_RE.match(d) for d in os.listdir(out_dir)):
        return
    with open(settings, encoding="utf-8") as f:
        before = recorded_enzyme_specs(json.load(f))
    for spec in specs:
        old = before.get(spec["code"])
        if old is None or _rule(old) != _rule(spec):
            recorded = ("; ".join(describe_enzyme(d) for d in before.values())
                        or "(no enzymes recorded)")
            raise ModConfigError(
                f"{out_dir} holds results made with enzyme definitions [{recorded}]; "
                f"this run uses [{describe_enzyme(spec)}]"
                + ("" if old is None else f" where it recorded [{describe_enzyme(old)}]")
                + ". Steps would be skipped and their outputs mixed. Use a new "
                  "--out-dir, or rerun with the original enzymes and missed cleavages.")


_LIB_ENZYME_DIR_RE = re.compile(r"^pass(?:1|3_semiExt)_(.+)$")


def _lib_built_specs(lib_path: str, settings: dict) -> List[dict]:
    """The definitions a library was built with: its enzyme's for a per-enzyme
    folder (pass1_<CODE>, pass3_semiExt_<CODE>), else the run's whole set."""
    recorded = recorded_enzyme_specs(settings)
    m = _LIB_ENZYME_DIR_RE.match(os.path.basename(os.path.dirname(os.path.abspath(lib_path))))
    if m:
        return [recorded[m.group(1)]] if m.group(1) in recorded else []
    return list(recorded.values())


def check_resume_mods(out_dir: str, mods: dict) -> None:
    """Refuse to resume into out_dir if its outputs used other modifications."""
    settings = os.path.join(out_dir, "pipeline.settings")
    if not os.path.exists(settings) or not any(
            _PASS_DIR_RE.match(d) for d in os.listdir(out_dir)):
        return
    with open(settings, encoding="utf-8") as f:
        old = json.load(f)
    before = mod_flags(old if "modifications" in old else None, include_extra=False)
    now    = mod_flags(mods, include_extra=False)
    if before != now:
        raise ModConfigError(
            f"{out_dir} holds results made with modifications [{describe_mods(before)}]; "
            f"this run uses [{describe_mods(now)}]. Steps would be skipped and their "
            f"outputs mixed. Use a new --out-dir, or rerun with the original "
            f"modifications (results from before modification support: --mods none).")


# ─── Utilities ────────────────────────────────────────────────────────────────

def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def for_windows_exe(cmd: list) -> list:
    """Under WSL, the argv a Windows executable (diann.exe ...) can read.

    The only place paths are converted: POSIX paths on /mnt/<drive>/ become
    Windows paths, a path inside the Linux filesystem is refused (ValueError,
    naming it) before anything runs. Unchanged unless on WSL with a Windows exe.
    """
    if not (wslpaths.is_wsl() and wslpaths.is_windows_exe(cmd[0])):
        return cmd
    return wslpaths.windows_argv(cmd)


def write_tsv(df: "pd.DataFrame", path: str):
    """Write df as TSV under a temporary name, then rename it into place.

    The restart checks accept any output that exists and has some size, so an
    interrupted write must never leave a partial file under the final name.
    """
    tmp = path + ".part"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


# DIA-NN writes its outputs under their final names, so an interrupted call
# (cancel, crash, power loss) can leave a truncated library or report that a
# restart would skip on. run_cmd records the outputs of the call in progress
# here; recover_interrupted() sets them aside on the next start.
INFLIGHT = "pipeline.inflight"


def diann_targets(cmd: list) -> List[str]:
    """The outputs of one DIA-NN call that the restart checks look at."""
    cmd = [str(c) for c in cmd]
    out = []
    for flag in ("--out", "--out-lib"):
        if flag in cmd[:-1]:
            out.append(cmd[cmd.index(flag) + 1])
    lib = cmd[cmd.index("--out-lib") + 1] if "--out-lib" in cmd[:-1] else ""
    if lib.endswith(".tsv"):                 # a --predictor call writes this one
        out.append(lib[:-4] + ".predicted.speclib")
    return out


def recover_interrupted(out_dir: str) -> List[str]:
    """Set aside the outputs of a DIA-NN call that did not finish.

    Each is renamed to <name>.interrupted (kept for inspection, never deleted),
    so the step runs again instead of being skipped. Returns the new names.
    """
    marker = os.path.join(out_dir, INFLIGHT)
    if not os.path.exists(marker):
        return []
    try:
        with open(marker, encoding="utf-8") as f:
            info = json.load(f)
    except (OSError, ValueError):
        info = {}
    moved = []
    for t in info.get("targets", []):
        if os.path.exists(t):
            dst = t + ".interrupted"
            if os.path.exists(dst):
                dst = f"{t}.{datetime.now():%Y%m%d-%H%M%S}.interrupted"
            os.replace(t, dst)
            moved.append(dst)
    print(f"  NOTE: the previous run stopped during '{info.get('label', '?')}'. "
          f"Its outputs were set aside and the step runs again: "
          f"{', '.join(os.path.basename(m) for m in moved) or 'none were written'}",
          flush=True)
    os.remove(marker)
    return moved


# DIA-NN's version, from the banner of the first DIA-NN call of this run (#38).
_DIANN_BANNER: Optional[str] = None


def _note_diann_banner(lines: List[str], log_path: str) -> None:
    """Record DIA-NN's banner once per run: print it and write it into the
    runtime block of pipeline.settings next to the log. Never raises: a missing
    or unreadable banner stays 'unknown'."""
    global _DIANN_BANNER
    try:
        parse = _build_info.parse_diann_banner if _build_info else (lambda _l: "unknown")
        banner = parse(lines)
        if banner == "unknown":
            return
        _DIANN_BANNER = banner
        print(f"  [runtime] DIA-NN: {banner}", flush=True)
        settings = os.path.join(os.path.dirname(os.path.abspath(log_path)),
                                "pipeline.settings")
        if not os.path.exists(settings):
            return
        with open(settings, encoding="utf-8") as f:
            data = json.load(f)
        rt = data.get("runtime")
        if not isinstance(rt, dict):
            return
        rt.setdefault("diann", {})["banner"] = banner
        with open(settings + ".part", "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        os.replace(settings + ".part", settings)
    except Exception as e:  # provenance must never fail a run
        print(f"  WARNING: could not record the DIA-NN banner ({e})", flush=True)


def run_cmd(cmd: list, log_path: str, label: str = "",
            diann_exe: str = "diann") -> Tuple[str, float]:
    """Run cmd (first element 'diann' replaced with diann_exe), tee to log.

    While it runs, its outputs are recorded in INFLIGHT next to the log; the
    record is removed only after a clean exit.
    """
    is_diann = bool(cmd) and str(cmd[0]) == "diann"
    cmd = [diann_exe if str(c) == "diann" else c for c in cmd]
    targets = diann_targets(cmd)
    marker = os.path.join(os.path.dirname(os.path.abspath(log_path)), INFLIGHT)
    if targets:
        with open(marker + ".part", "w", encoding="utf-8") as f:
            json.dump({"label": label, "targets": targets}, f)
        os.replace(marker + ".part", marker)
    cmd = for_windows_exe(cmd)
    cmd_str = " ".join(
        f'"{c}"' if (" " in str(c) and not str(c).startswith('"')) else str(c)
        for c in cmd
    )
    print(f"\n[{ts()}] {label}", flush=True)
    print(f"  $ {cmd_str}\n", flush=True)
    t0 = time.time()
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    with open(log_path, "a", encoding="utf-8") as lf:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", env=env,
        )
        head: List[str] = []
        want_banner = is_diann and _DIANN_BANNER is None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            lf.write(line)
            if want_banner and len(head) < 5:
                head.append(line)
                if len(head) == 5:
                    _note_diann_banner(head, log_path); want_banner = False
        proc.wait()
        if want_banner:
            _note_diann_banner(head, log_path)
    elapsed = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit {proc.returncode}): {cmd_str}")
    if targets:
        os.remove(marker)
    print(f"  [{ts()}] Done in {elapsed / 60:.1f} min ({elapsed:.0f} s)", flush=True)
    return cmd_str, elapsed


def parse_diann_log_params(log_path: str) -> dict:
    """Extract scan window, MS1 acc, MS2 acc from a DIA-NN log (ceiling-rounded)."""
    params: dict = {"window": None, "ms1_acc": None, "mass_acc": None}
    if not os.path.exists(log_path):
        return params
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if "Scan window radius set to" in line:
                m = re.search(r"Scan window radius set to\s+(\d+)", line)
                if m:
                    params["window"] = int(m.group(1))
            if "Recommended MS1 mass accuracy setting:" in line:
                m = re.search(r"([\d.]+)\s*ppm", line)
                if m:
                    params["ms1_acc"] = math.ceil(float(m.group(1)))
            if "Optimised mass accuracy:" in line:
                m = re.search(r"([\d.]+)\s*ppm", line)
                if m:
                    params["mass_acc"] = math.ceil(float(m.group(1)))
    return params


def find_speclib(
    per_enzyme_speclibs: Dict[str, str],
    enzyme_code: str,
    enzyme_name: str,
    speclib_dir: Optional[str] = None,
) -> Optional[str]:
    """
    Locate a pre-predicted .predicted.speclib for this enzyme.
    Priority:
      1. Exact path in per_enzyme_speclibs[enzyme_code]
      2. Directory scan of speclib_dir for enzyme_name prefix match
    """
    # 1. Explicit per-enzyme path
    path = per_enzyme_speclibs.get(enzyme_code, "").strip()
    if path and os.path.exists(path) and os.path.getsize(path) > 100:
        return path

    # 2. Fallback: scan speclib_dir
    if speclib_dir and os.path.isdir(speclib_dir):
        stem = file_stem(enzyme_name)
        # another enzyme whose name extends this one ('elastase_extended' for
        # 'elastase') must not have its libraries taken for this enzyme's
        longer = [e.stem for e in ENZYME_REGISTRY.values()
                  if e.stem != stem and e.stem.startswith(stem)]
        for root, _, files in os.walk(speclib_dir):
            for f in sorted(files):
                fl = f.lower()
                if ((fl.startswith(enzyme_name.lower()) or fl.startswith(stem))
                        and not any(fl.startswith(x) for x in longer)
                        and f.endswith(".predicted.speclib")
                        and os.path.getsize(os.path.join(root, f)) > 100):
                    return os.path.join(root, f)
    return None


def combine_libraries(libs_dict: Dict[str, str], out_path: str) -> pd.DataFrame:
    """Concatenate per-enzyme empirical TSV libraries."""
    dfs = []
    for name, path in libs_dict.items():
        if not (os.path.exists(path) and os.path.getsize(path) > 100):
            print(f"  WARNING: library missing/empty for {name}: {path}")
            continue
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if len(df) == 0:
            continue
        if "Enzyme" not in df.columns:
            df.insert(0, "Enzyme", name)
        n = df["transition_group_id"].nunique() if "transition_group_id" in df.columns else "?"
        print(f"  {name}: {len(df):,} rows, {n} precursors")
        dfs.append(df)
    if not dfs:
        raise RuntimeError("No empirical libraries to combine.")
    combined = pd.concat(dfs, ignore_index=True)
    write_tsv(combined, out_path)
    n = combined["transition_group_id"].nunique() if "transition_group_id" in combined.columns else "?"
    print(f"  Combined: {len(combined):,} rows, {n} unique precursors → {os.path.basename(out_path)}")
    return combined


def _emit(tag: str, data: dict):
    print(f"##{tag}## {json.dumps(data)}", flush=True)


def _find_plotter(diann_exe: str) -> Optional[str]:
    """Return path to dia-nn-plotter: check same dir as diann_exe first, then PATH."""
    diann_path = os.path.abspath(diann_exe) if os.sep in diann_exe else None
    if diann_path and os.path.isfile(diann_path):
        parent = os.path.dirname(diann_path)
        for name in ("dia-nn-plotter.exe", "dia-nn-plotter"):
            candidate = os.path.join(parent, name)
            if os.path.isfile(candidate):
                return candidate
    for name in ("dia-nn-plotter", "dia-nn-plotter.exe"):
        found = shutil.which(name)
        if found:
            return found
    return None


def run_plotter(report_dir: str, diann_exe: str) -> None:
    """
    Run dia-nn-plotter on report.stats.tsv + report.tsv → report.pdf.
    Non-fatal: warnings are printed but the pipeline continues.
    Skipped silently if report.pdf already exists (idempotent on resume).
    """
    stats_tsv  = os.path.join(report_dir, "report.stats.tsv")
    report_tsv = os.path.join(report_dir, "report.tsv")
    report_pdf = os.path.join(report_dir, "report.pdf")

    if not (os.path.exists(stats_tsv) and os.path.exists(report_tsv)):
        return
    if os.path.exists(report_pdf):
        print(f"  [plotter] {os.path.basename(report_pdf)} already exists — skipped")
        return

    plotter = _find_plotter(diann_exe)
    if not plotter:
        print(f"  [plotter] dia-nn-plotter not found — skipping PDF report")
        return

    print(f"\n  [plotter] Generating {os.path.basename(report_pdf)} …", flush=True)
    try:
        proc = subprocess.run(
            for_windows_exe([plotter, stats_tsv, report_tsv, report_pdf]),
            capture_output=True, text=True, timeout=180,
        )
        if proc.returncode == 0:
            print(f"  [plotter] Done → {os.path.relpath(report_pdf, report_dir)}", flush=True)
        else:
            msg = (proc.stderr or proc.stdout or "").strip()[:300]
            print(f"  [plotter] WARNING: exit {proc.returncode}: {msg}", flush=True)
    except subprocess.TimeoutExpired:
        print(f"  [plotter] WARNING: timed out after 180 s", flush=True)
    except Exception as e:
        print(f"  [plotter] WARNING: {e}", flush=True)


def emit_chart_data(report_path: str, combined_lib_path: Optional[str],
                    pass_label: str) -> None:
    """
    Read a DIA-NN report and emit ##CHARTDATA## JSON lines for:
      - per_run : precursors + peptides per raw file
      - enzyme_breakdown : unique precursors per enzyme (needs Enzyme col in lib)
    """
    if not os.path.exists(report_path):
        return
    try:
        df = pd.read_csv(report_path, sep="\t", low_memory=False)
    except Exception as e:
        print(f"  WARNING: chart_data read failed ({e})")
        return

    q_col = "Q.Value" if "Q.Value" in df.columns else None
    fdr   = df[df[q_col] <= 0.01].copy() if q_col else df.copy()

    run_col  = next((c for c in ("Run", "File.Name") if c in fdr.columns), None)
    prec_col = "Precursor.Id"       if "Precursor.Id"       in fdr.columns else None
    pep_col  = ("Modified.Sequence" if "Modified.Sequence"  in fdr.columns else
                "Stripped.Sequence" if "Stripped.Sequence"  in fdr.columns else None)
    pg_col   = ("Protein.Group"     if "Protein.Group"      in fdr.columns else
                "ProteinGroup"      if "ProteinGroup"        in fdr.columns else None)

    # ── Per-run counts ────────────────────────────────────────────────────────
    if run_col and (prec_col or pep_col):
        runs        = sorted(fdr[run_col].unique())
        run_labels  = [os.path.splitext(os.path.basename(str(r)))[0] for r in runs]
        per_prec    = [int(fdr[fdr[run_col] == r][prec_col].nunique())
                       if prec_col else 0 for r in runs]
        per_pep     = [int(fdr[fdr[run_col] == r][pep_col].nunique())
                       if pep_col else 0 for r in runs]
        per_prot    = [int(fdr[fdr[run_col] == r][pg_col].nunique())
                       if pg_col else 0 for r in runs]
        _emit("CHARTDATA", {
            "pass":       pass_label,
            "type":       "per_run",
            "runs":       run_labels,
            "precursors": per_prec,
            "peptides":   per_pep,   # Modified.Sequence (falls back to Stripped.Sequence)
            "proteins":   per_prot,  # unique Protein.Group per run
        })

    # ── Enzyme breakdown ──────────────────────────────────────────────────────
    if combined_lib_path and os.path.exists(combined_lib_path) and prec_col:
        try:
            lib_cols = lambda c: c in ("transition_group_id", "Enzyme")
            lib = pd.read_csv(combined_lib_path, sep="\t", low_memory=False,
                              usecols=lib_cols)
            if "transition_group_id" in lib.columns and "Enzyme" in lib.columns:
                enz_map    = (lib.drop_duplicates("transition_group_id")
                                 .set_index("transition_group_id")["Enzyme"].to_dict())
                enz_series = fdr[prec_col].map(enz_map)
                enz_counts = enz_series.value_counts().dropna()
                if len(enz_counts):
                    _emit("CHARTDATA", {
                        "pass":    pass_label,
                        "type":    "enzyme_breakdown",
                        "enzymes": list(enz_counts.index),
                        "counts":  [int(v) for v in enz_counts.values],
                    })
        except Exception as e:
            print(f"  WARNING: enzyme chart data failed ({e})")


# ─── Virtual FASTA helpers ────────────────────────────────────────────────────

def _parse_fasta(fasta_path: str) -> Dict[str, str]:
    proteins: Dict[str, str] = {}
    header = seq_parts = None
    with open(fasta_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    proteins[header] = "".join(seq_parts)
                fw = line[1:].split()[0]
                parts = fw.split("|")
                header = parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else fw
                seq_parts = []
            elif header is not None:
                seq_parts.append(line.upper())
    if header is not None:
        proteins[header] = "".join(seq_parts)
    return proteins


def _find_cleavage_sites(seq: str, enz: EnzymeDef) -> List[int]:
    """
    Return peptide-boundary positions for the given enzyme sense.

    sense=1 (C-terminal cutter): site is placed AFTER the cut residue.
      Cleavage blocked if the NEXT residue is in enz.noclv.
    sense=0 (N-terminal cutter): site is placed BEFORE the cut residue.
      The cut residue begins the downstream peptide.
    """
    cut_set   = set(enz.cut_res)
    noclv_set = set(enz.noclv)
    sites     = [0]

    if enz.sense == 1:
        for i, aa in enumerate(seq):
            if aa in cut_set:
                nxt = seq[i + 1] if i + 1 < len(seq) else None
                if nxt is None or nxt not in noclv_set:
                    sites.append(i + 1)
    else:  # sense == 0: N-terminal cutter
        for i, aa in enumerate(seq):
            if aa in cut_set and i > 0:
                sites.append(i)

    if sites[-1] != len(seq):
        sites.append(len(seq))
    return sites


def _digest_semi(
    seq: str,
    prot_id: str,
    enz: EnzymeDef,
    min_len: int = 7,
    max_len: int = 30,
    missed: int = 1,
) -> List[Tuple[str, str, str]]:
    """
    Generate all semi-specific peptides for one protein + enzyme.
    Returns list of (peptide_seq, prot_id, specificity).
    Works for both C-terminal (sense=1) and N-terminal (sense=0) cutters.
    """
    sites = _find_cleavage_sites(seq, enz)
    site_set = set(sites)
    results: List[Tuple[str, str, str]] = []
    n = len(sites)

    def add(pep: str, spec: str, start: Optional[int] = None):
        # Specificity is decided by the peptide's actual termini, not by the
        # loop that produced it: a truncation of a missed-cleavage span can land
        # on a cleavage site and is then fully specific (bug before 2026-09-19).
        if start is not None:
            end = start + len(pep)
            spec = ("fully_specific" if start in site_set and end in site_set
                    else "semi_specific")
        if min_len <= len(pep) <= max_len:
            results.append((pep, prot_id, spec))

    for i in range(n - 1):
        for j in range(i + 1, min(i + missed + 2, n)):
            start, end = sites[i], sites[j]
            # Fully specific (both termini at cleavage sites)
            add(seq[start:end], "fully_specific", start)
            # N-truncated (N-terminus moved; semi unless it lands on a site)
            for new_start in range(start + 1, end - min_len + 1):
                add(seq[new_start:end], "semi_specific", new_start)
        if i + 1 < n:
            next_site = sites[min(i + missed + 1, n - 1)]
            # C-truncated (C-terminus moved; semi unless it lands on a site)
            for new_end in range(sites[i] + min_len, next_site):
                if new_end != sites[i + 1]:
                    add(seq[sites[i]:new_end], "semi_specific", sites[i])
    return results


def generate_semi_vfasta(
    fasta_path: str,
    enzyme_code: str,
    enz: EnzymeDef,
    out_vfasta: str,
    out_map: str,
    protein_whitelist: Optional[set] = None,
    min_len: int = 7,
    max_len: int = 30,
    missed: int = 1,
) -> int:
    """Generate semi-specific virtual FASTA + mapping TSV. Returns peptide count."""
    proteins = _parse_fasta(fasta_path)
    if protein_whitelist is not None:
        before  = len(proteins)
        after_d = {k: v for k, v in proteins.items() if k in protein_whitelist}
        print(f"  [{enz.name}] Protein filter: {before:,} → {len(after_d):,}")
        if len(after_d) == 0 and before > 0:
            _fasta_sample = sorted(proteins.keys())[:5]
            _wl_sample    = sorted(protein_whitelist)[:5]
            print(f"  [{enz.name}] ERROR: FASTA/analysis ID mismatch — 0 proteins survived filter.")
            print(f"  [{enz.name}]   FASTA IDs (first 5) : {_fasta_sample}")
            print(f"  [{enz.name}]   Whitelist IDs (first 5): {_wl_sample}")
            print(f"  [{enz.name}]   Check that the FASTA matches the organism of the speclib/analysis.")
        proteins = after_d

    pep_to_prots: Dict[str, set] = defaultdict(set)
    pep_to_spec:  Dict[str, str] = {}

    for prot_id, seq in proteins.items():
        for pep, pid, spec in _digest_semi(seq, prot_id, enz, min_len, max_len, missed):
            pep_to_prots[pep].add(pid)
            # fully specific in ANY protein context -> fully specific
            if pep_to_spec.get(pep) != "fully_specific":
                pep_to_spec[pep] = spec

    map_rows = []
    with open(out_vfasta + ".part", "w") as ffa:
        for idx, (pep, prots) in enumerate(sorted(pep_to_prots.items())):
            vhdr = f"vpep_{enzyme_code}_{idx:07d}"
            ffa.write(f">{vhdr}\n{pep}\n")
            map_rows.append({
                "VirtualHeader":   vhdr,
                "PeptideSequence": pep,
                "ProteinGroup":    ";".join(sorted(prots)),
                "Enzyme":          enz.name,
                "Specificity":     pep_to_spec[pep],
            })

    _MAP_COLS = ["VirtualHeader", "PeptideSequence", "ProteinGroup", "Enzyme", "Specificity"]
    map_df    = pd.DataFrame(map_rows, columns=_MAP_COLS) if map_rows else pd.DataFrame(columns=_MAP_COLS)
    write_tsv(map_df, out_map)
    os.replace(out_vfasta + ".part", out_vfasta)   # last: the restart check needs both
    n_semi = int(map_df["Specificity"].eq("semi_specific").sum())
    n_full = int(map_df["Specificity"].eq("fully_specific").sum())
    print(f"  [{enz.name}] {len(map_rows):,} peptides "
          f"({n_semi:,} semi, {n_full:,} fully-specific)")
    return len(map_rows)


def _reannotate_protein_groups(df: pd.DataFrame,
                                vmap_df: Optional[pd.DataFrame]) -> pd.DataFrame:
    if vmap_df is None or "PeptideSequence" not in df.columns:
        return df
    lookup = vmap_df.set_index("PeptideSequence")["ProteinGroup"].to_dict()
    for col in ("ProteinGroup", "UniprotID"):
        if col not in df.columns:
            continue
        mask = df[col].astype(str).str.startswith("vpep_")
        if not mask.any():
            continue
        remapped = df.loc[mask, "PeptideSequence"].map(lookup)
        df.loc[mask, col] = remapped.where(remapped.notna(), df.loc[mask, col])
    return df


# ─── Summary helpers ──────────────────────────────────────────────────────────

def print_id_summary(report_path: str, label: str):
    if not os.path.exists(report_path):
        print(f"  (report not found: {report_path})")
        return
    try:
        df     = pd.read_csv(report_path, sep="\t", low_memory=False)
        q_col  = "Q.Value"         if "Q.Value"         in df.columns else None
        pq_col = ("Protein.Q.Value" if "Protein.Q.Value" in df.columns
                  else "PG.Q.Value" if "PG.Q.Value"      in df.columns else None)
        flt    = df
        if q_col and pq_col:
            flt = df[(df[q_col] <= 0.01) & (df[pq_col] <= 0.01)]
        elif q_col:
            flt = df[df[q_col] <= 0.01]
        pg_col  = next((c for c in ("Protein.Group", "ProteinGroup")
                        if c in flt.columns), None)
        n_prot  = flt[pg_col].nunique()              if pg_col else "?"
        n_pep   = (flt["Stripped.Sequence"].nunique()
                   if "Stripped.Sequence" in flt.columns else "?")
        n_prec  = (flt["Precursor.Id"].nunique()     if "Precursor.Id"     in flt.columns else
                   flt["transition_group_id"].nunique()
                   if "transition_group_id" in flt.columns else "?")
        print(f"\n  ── {label} (q ≤ 0.01) ──────────────")
        print(f"     Proteins  : {n_prot}")
        print(f"     Peptides  : {n_pep}")
        print(f"     Precursors: {n_prec}")
        _emit("METRIC", {"label": label, "proteins": n_prot,
                          "peptides": n_pep, "precursors": n_prec})
    except Exception as e:
        print(f"  (could not parse report: {e})")


# ─── Pass 1 ───────────────────────────────────────────────────────────────────

def run_pass1(
    raw_files: List[str],
    out_dir: str,
    enzyme_codes: List[str],
    fasta: str,
    per_enzyme_speclibs: Dict[str, str],
    speclib_dir: Optional[str],
    threads: int,
    diann_exe: str,
    log_path: str,
    timing: dict,
    min_pep_len: int = 7,
    max_pep_len: int = 30,
    missed_cleavages: int = 1,
    mods: Optional[dict] = None,
    enzyme_mc: Optional[Dict[str, int]] = None,
) -> Dict[str, Tuple[str, dict]]:
    """
    Predict + search per enzyme. enzyme_mc: effective missed cleavages per code
    (resolve_enzyme_mc); an enzyme not in it uses missed_cleavages.
    Per-enzyme speclib looked up first from per_enzyme_speclibs, then speclib_dir scan.
    Returns {code: (emp_lib_path, log_params_dict)}.
    """
    results: Dict[str, Tuple[str, dict]] = {}
    f_args: List[str] = []
    for r in raw_files:
        f_args += ["--f", r]
    mod_args  = mod_flags(mods)
    lib_flags = mod_flags(mods, include_extra=False)

    for code in enzyme_codes:
        enz    = ENZYME_REGISTRY[code]
        mc     = (enzyme_mc or {}).get(code, missed_cleavages)
        spec   = [enzyme_spec(code, mc)]
        p1_dir = os.path.join(out_dir, f"pass1_{code}")
        os.makedirs(p1_dir, exist_ok=True)

        emp_lib = os.path.join(p1_dir, f"{enz.stem}_empirical.tsv")
        report1 = os.path.join(p1_dir, "report.tsv")
        log1    = report1.replace(".tsv", ".log.txt")

        _emit("STATUS", {"pass": f"pass1_{code}", "status": "running"})

        if os.path.exists(emp_lib) and os.path.getsize(emp_lib) > 100:
            check_lib_mods(emp_lib, lib_flags, f"Pass1 {enz.name} empirical library",
                           enzymes=spec)
            print(f"\n  [SKIP pass1_{code}] empirical lib exists")
            params = parse_diann_log_params(log1)
            results[code] = (emp_lib, params)
            _emit("STATUS", {"pass": f"pass1_{code}", "status": "done", "skipped": True})
            continue

        # ── Locate or predict speclib ─────────────────────────────────────────
        speclib = find_speclib(per_enzyme_speclibs, code, enz.name, speclib_dir)
        if speclib:
            check_lib_mods(speclib, lib_flags, f"Pass1 {enz.name} prebuilt speclib",
                           enzymes=spec)
            print(f"  [{enz.name}] Reusing speclib: {os.path.basename(speclib)}")
        else:
            out_lib_pred = os.path.join(p1_dir, f"{enz.stem}_predicted.tsv")
            speclib_auto = out_lib_pred.replace(".tsv", ".predicted.speclib")
            if not (os.path.exists(speclib_auto) and os.path.getsize(speclib_auto) > 100):
                _, elapsed = run_cmd([
                    "diann",
                    "--fasta",            fasta,
                    "--fasta-search",
                    "--cut",              enz.cut,
                    "--missed-cleavages", str(mc),
                    "--min-pep-len",      str(min_pep_len),
                    "--max-pep-len",      str(max_pep_len),
                    *mod_args,
                    "--predictor",
                    "--gen-spec-lib",
                    "--out-lib",          out_lib_pred,
                    "--out",              os.path.join(p1_dir, "report_predict.tsv"),
                    "--threads",          str(threads),
                ], log_path, f"Pass1 predict — {enz.name} [sense={enz.sense}, cut={enz.cut}]",
                             diann_exe)
                timing[f"pass1_predict_{code}"] = elapsed
                speclib = speclib_auto if os.path.exists(speclib_auto) else out_lib_pred
            else:
                check_lib_mods(speclib_auto, lib_flags, f"Pass1 {enz.name} predicted speclib",
                               enzymes=spec)
                print(f"  [SKIP predict] {enz.name}: speclib exists")
            speclib = speclib_auto if os.path.exists(speclib_auto) else out_lib_pred

        # ── Search pass1 ─────────────────────────────────────────────────────
        _, elapsed = run_cmd([
            "diann",
            *f_args,
            "--lib",            speclib,
            *mod_args,
            "--gen-spec-lib",
            "--full-profiling",
            "--smart-profiling",
            "--out-lib",        emp_lib,
            "--out",            report1,
            "--threads",        str(threads),
        ], log_path, f"Pass1 search — {enz.name}", diann_exe)
        timing[f"pass1_search_{code}"] = elapsed
        run_plotter(p1_dir, diann_exe)

        params = parse_diann_log_params(log1)
        print(f"  [{enz.name}] Optimised: window={params['window']}, "
              f"ms1={params['ms1_acc']} ppm, ms2={params['mass_acc']} ppm")
        results[code] = (emp_lib, params)
        _emit("STATUS", {"pass": f"pass1_{code}", "status": "done",
                          "elapsed_s": round(elapsed, 1), **params})
    return results


# ─── Pass 2 ───────────────────────────────────────────────────────────────────

def run_pass2(
    raw_files: List[str],
    out_dir: str,
    enzyme_codes: List[str],
    combined_lib: str,
    fasta: str,
    threads: int,
    diann_exe: str,
    log_path: str,
    timing: dict,
    mods: Optional[dict] = None,
) -> Tuple[str, dict]:
    """
    Combined library search + matrices.
    --fasta enables Protein.Names / Genes annotation in matrices (no --reannotate).
    Returns (report_path, log_params).
    """
    codes_str = "".join(enzyme_codes)
    p2_dir    = os.path.join(out_dir, f"pass2_{codes_str}")
    os.makedirs(p2_dir, exist_ok=True)
    report2   = os.path.join(p2_dir, "report.tsv")
    log2      = report2.replace(".tsv", ".log.txt")

    _emit("STATUS", {"pass": f"pass2_{codes_str}", "status": "running"})

    if os.path.exists(report2) and os.path.getsize(report2) > 100:
        print(f"\n  [SKIP pass2] report exists")
        params = parse_diann_log_params(log2)
        _emit("STATUS", {"pass": f"pass2_{codes_str}", "status": "done", "skipped": True})
        return report2, params

    f_args: List[str] = []
    for r in raw_files:
        f_args += ["--f", r]

    _, elapsed = run_cmd([
        "diann",
        *f_args,
        "--lib",            combined_lib,
        "--fasta",          fasta,    # populates Protein.Names / Genes; NO --reannotate
        "--cut",            "",       # permissive: accept all peptide types in lib
        *mod_flags(mods),
        "--matrices",
        "--matrix-qvalue",  "0.01",
        "--matrix-spec-q",  "0.01",
        "--out",            report2,
        "--threads",        str(threads),
    ], log_path, f"Pass2 — combined {codes_str} library + matrices", diann_exe)
    timing[f"pass2_{codes_str}"] = elapsed
    run_plotter(p2_dir, diann_exe)

    params = parse_diann_log_params(log2)
    _emit("STATUS", {"pass": f"pass2_{codes_str}", "status": "done",
                      "elapsed_s": round(elapsed, 1), **params})
    return report2, params


# ─── Pass 3 ───────────────────────────────────────────────────────────────────

def run_pass3(
    raw_files: List[str],
    out_dir: str,
    enzyme_codes: List[str],
    pass1_results: Dict[str, Tuple[str, dict]],
    fasta: str,
    threads: int,
    diann_exe: str,
    log_path: str,
    timing: dict,
    min_pep_len: int = 7,
    max_pep_len: int = 30,
    missed_cleavages: int = 1,
    protein_whitelist: Optional[set] = None,
    mods: Optional[dict] = None,
    enzyme_mc: Optional[Dict[str, int]] = None,
) -> Dict[str, Tuple[str, Optional[pd.DataFrame]]]:
    """
    Per-enzyme semi-specific extension via virtual FASTA.
    Inline FASTA reannotation (--reannotate) during pass3 search ensures empirical
    libraries have correct ProteinGroup / Genes / ProteinName before combining.
    Fixed mass accuracy params propagated from pass1 per enzyme.
    Returns {code: (emp_lib_path, vmap_df)}.
    """
    results: Dict[str, Tuple[str, Optional[pd.DataFrame]]] = {}
    f_args: List[str] = []
    for r in raw_files:
        f_args += ["--f", r]
    # Each virtual-FASTA entry is one candidate peptide, so protein-level options
    # (protein N-term mods, --met-excision) would hit every peptide: left out here.
    mod_args  = mod_flags(mods, protein_level=False)
    lib_flags = mod_flags(mods, protein_level=False, include_extra=False)

    for code in enzyme_codes:
        enz  = ENZYME_REGISTRY[code]
        # The limit takes effect here, in the Python digest that builds the
        # virtual FASTA; the pass-3 prediction itself runs with 0 (one peptide
        # per entry).
        mc   = (enzyme_mc or {}).get(code, missed_cleavages)
        spec = [enzyme_spec(code, mc)]

        p3_dir = os.path.join(out_dir, f"pass3_semiExt_{code}")
        os.makedirs(p3_dir, exist_ok=True)

        # Fixed params from pass1 per enzyme (None → DIA-NN auto-optimises)
        _, p1_params = pass1_results.get(code, (None, {}))
        window   = p1_params.get("window")   if p1_params else None
        ms1_acc  = p1_params.get("ms1_acc")  if p1_params else None
        mass_acc = p1_params.get("mass_acc") if p1_params else None

        _emit("STATUS", {"pass": f"pass3_semiExt_{code}", "status": "running"})

        # This enzyme's own empirical library, defined before any early exit:
        # an empty virtual FASTA must not report another enzyme's library.
        emp_lib = os.path.join(p3_dir, f"{enz.stem}_semi_empirical.tsv")

        # Virtual FASTA
        vfasta = os.path.join(p3_dir, f"{code}_semi_vfasta.fasta")
        vmap   = os.path.join(p3_dir, f"{code}_semi_vfasta_map.tsv")
        if not (os.path.exists(vfasta) and os.path.exists(vmap)):
            print(f"\n[{ts()}] Generating virtual FASTA for {enz.name} "
                  f"(sense={'C-terminal' if enz.sense else 'N-terminal'})...")
            n_pep = generate_semi_vfasta(fasta, code, enz, vfasta, vmap,
                                         protein_whitelist=protein_whitelist,
                                         min_len=min_pep_len, max_len=max_pep_len,
                                         missed=mc)
        else:
            print(f"  [SKIP vfasta] {enz.name}: {os.path.basename(vfasta)} exists")
            vmap_df_check = pd.read_csv(vmap, sep="\t")
            n_pep = len(vmap_df_check)

        # If the virtual FASTA is empty (e.g. FASTA/speclib organism mismatch) skip pass3
        if n_pep == 0:
            print(f"  [SKIP pass3_{code}] virtual FASTA is empty — "
                  f"verify the FASTA matches the analysis organism.")
            results[code] = (emp_lib, None)
            _emit("STATUS", {"pass": f"pass3_semiExt_{code}", "status": "skipped",
                              "skipped": True})
            continue

        vmap_df = pd.read_csv(vmap, sep="\t") if os.path.exists(vmap) else None

        # Predict virtual speclib
        out_tsv      = os.path.join(p3_dir, f"{enz.stem}_semi_predicted.tsv")
        semi_speclib = out_tsv.replace(".tsv", ".predicted.speclib")
        if not (os.path.exists(semi_speclib) and os.path.getsize(semi_speclib) > 100):
            _, elapsed = run_cmd([
                "diann",
                "--fasta",            vfasta,
                "--fasta-search",
                "--cut",              "",       # no cut: one peptide per virtual entry
                "--missed-cleavages", "0",
                "--min-pep-len",      str(min_pep_len),
                "--max-pep-len",      str(max_pep_len),
                *mod_args,
                "--predictor",
                "--gen-spec-lib",
                "--out-lib",          out_tsv,
                "--out",              os.path.join(p3_dir, "report_predict.tsv"),
                "--threads",          str(threads),
            ], log_path, f"Pass3 predict virtual FASTA — {enz.name}", diann_exe)
            timing[f"pass3_predict_{code}"] = elapsed
        else:
            check_lib_mods(semi_speclib, lib_flags, f"Pass3 {enz.name} predicted speclib",
                           protein_level=False, enzymes=spec)
            print(f"  [SKIP predict] {enz.name}: {os.path.basename(semi_speclib)} exists")

        # Pass3 search — inline reannotation maps vpep_ headers → real accessions
        report3 = os.path.join(p3_dir, "report.tsv")

        if os.path.exists(emp_lib) and os.path.getsize(emp_lib) > 100:
            check_lib_mods(emp_lib, lib_flags, f"Pass3 {enz.name} empirical library",
                           protein_level=False, enzymes=spec)
            print(f"  [SKIP pass3_{code}] empirical lib exists")
        else:
            cmd = [
                "diann",
                *f_args,
                "--lib",            semi_speclib,
                "--fasta",          fasta,
                "--reannotate",     # maps vpep_ headers → real UniProt IDs (pass3 only)
                "--cut",            enz.cut,
                "--semi",
                *mod_args,
                "--gen-spec-lib",
                "--full-profiling",
                "--smart-profiling",
                "--no-prot-inf",
                "--out-lib",        emp_lib,
                "--out",            report3,
                "--threads",        str(threads),
            ]
            if window   is not None: cmd += ["--window",       str(window)]
            if ms1_acc  is not None: cmd += ["--mass-acc-ms1", str(ms1_acc)]
            if mass_acc is not None: cmd += ["--mass-acc",      str(mass_acc)]
            pstr = (f"[w={window}, ms1={ms1_acc}, ms2={mass_acc} ppm]"
                    if window is not None else "[auto-optimise]")
            _, elapsed = run_cmd(
                cmd, log_path,
                f"Pass3 semi search (inline reannotate) — {enz.name} {pstr}",
                diann_exe,
            )
            timing[f"pass3_search_{code}"] = elapsed
            run_plotter(p3_dir, diann_exe)

        results[code] = (emp_lib, vmap_df)
        _emit("STATUS", {"pass": f"pass3_semiExt_{code}", "status": "done"})

    return results


# ─── Pass 4 ───────────────────────────────────────────────────────────────────

def run_pass4(
    raw_files: List[str],
    out_dir: str,
    enzyme_codes: List[str],
    ext_combined: str,
    pass2_params: dict,
    fasta: str,
    threads: int,
    diann_exe: str,
    log_path: str,
    timing: dict,
    mods: Optional[dict] = None,
) -> str:
    """
    Combined semi-ext search + matrices.
    --fasta populates Protein.Names / Genes in matrix output (no --reannotate).
    Fixed params propagated from pass2 log.
    Returns report path.
    """
    codes_str = "".join(enzyme_codes)
    p4_dir    = os.path.join(out_dir, f"pass4_semiExt_{codes_str}")
    os.makedirs(p4_dir, exist_ok=True)
    report4   = os.path.join(p4_dir, "report.tsv")

    window   = pass2_params.get("window")
    ms1_acc  = pass2_params.get("ms1_acc")
    mass_acc = pass2_params.get("mass_acc")

    _emit("STATUS", {"pass": f"pass4_semiExt_{codes_str}", "status": "running"})

    if os.path.exists(report4) and os.path.getsize(report4) > 100:
        print(f"\n  [SKIP pass4] report exists")
        _emit("STATUS", {"pass": f"pass4_semiExt_{codes_str}",
                          "status": "done", "skipped": True})
        return report4

    f_args: List[str] = []
    for r in raw_files:
        f_args += ["--f", r]

    cmd = [
        "diann",
        *f_args,
        "--lib",            ext_combined,
        "--fasta",          fasta,    # populates annotation in matrices; NO --reannotate
        "--cut",            "",       # permissive: lib already annotated upstream (pass3)
        *mod_flags(mods),
        "--matrices",
        "--matrix-qvalue",  "0.01",
        "--matrix-spec-q",  "0.01",
        "--out",            report4,
        "--threads",        str(threads),
    ]
    if window   is not None: cmd += ["--window",       str(window)]
    if ms1_acc  is not None: cmd += ["--mass-acc-ms1", str(ms1_acc)]
    if mass_acc is not None: cmd += ["--mass-acc",      str(mass_acc)]

    pstr = (f"[w={window}, ms1={ms1_acc}, ms2={mass_acc} ppm]"
            if window is not None else "[auto-optimise]")
    _, elapsed = run_cmd(
        cmd, log_path,
        f"Pass4 — semi-ext {codes_str} + matrices {pstr}",
        diann_exe,
    )
    timing[f"pass4_{codes_str}"] = elapsed
    run_plotter(p4_dir, diann_exe)
    _emit("STATUS", {"pass": f"pass4_semiExt_{codes_str}",
                      "status": "done", "elapsed_s": round(elapsed, 1)})
    return report4


def _copy_fasta_to_out(fasta_path: Optional[str], out_dir: str) -> Optional[str]:
    """Copy the analysis FASTA into ``out_dir`` as report.used.fasta.

    Idempotent: if destination exists with same size, skip the copy.
    Returns the destination path (or None if nothing was copied).
    """
    if not fasta_path or not os.path.exists(fasta_path):
        return None
    dst = os.path.join(out_dir, "report.used.fasta")
    try:
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(fasta_path):
            return dst
        shutil.copyfile(fasta_path, dst)
        print(f"  FASTA copy: {os.path.basename(dst)}")
        return dst
    except OSError as e:
        print(f"  WARNING: failed to copy FASTA to {dst}: {e}")
        return None


# In-source-decay filter of the proteoform step (issue #42); set from the CLI in main()
PF_ISD = {"enabled": True, "rtWindowMin": 0.05}


def _run_proteoform(report_annotated: str, fasta_used: Optional[str],
                    out_dir: str, min_support: int = 2) -> None:
    """Run proteoform_analysis on the annotated report. No-op if the module
    failed to import, inputs are missing, or the call raises."""
    if _pf is None:
        return
    if not (os.path.exists(report_annotated) and fasta_used and os.path.exists(fasta_used)):
        return
    print(f"\n  Proteoform analysis (min_support={min_support}, ISD filter "
          f"{'on, ' + str(PF_ISD['rtWindowMin']) + ' min' if PF_ISD['enabled'] else 'off'})...")
    try:
        _pf.run(report_annotated, fasta_used, out_dir, min_support=min_support,
                isd_filter=PF_ISD["enabled"], isd_rt_window=PF_ISD["rtWindowMin"])
    except Exception as e:
        print(f"  WARNING: proteoform analysis failed: {e}")


def _copy_peptide_mapper(out_dir: str, report_annotated: str,
                         fasta_used: Optional[str]) -> None:
    """Copy bin/viewer/viewer.html to ``out_dir`` as PeptideMapper.html.

    The viewer auto-fetches ``./report_annotated.tsv`` and ``./report.used.fasta``
    as sibling files when opened over HTTP (e.g. served by mEDIA.exe).
    Reports can easily reach multi-GB scale, so we do NOT inline the data;
    we ship the small template (~90 KB) and let the browser stream the
    siblings on demand. report_annotated / fasta_used args are kept for
    signature compatibility but unused.
    """
    del report_annotated, fasta_used  # signature compatibility
    viewer_src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "viewer", "viewer.html")
    if not os.path.exists(viewer_src):
        return
    try:
        dst = os.path.join(out_dir, "PeptideMapper.html")
        shutil.copyfile(viewer_src, dst)
        size_kb = os.path.getsize(dst) / 1024
        print(f"  Peptide mapper: {os.path.basename(dst)} ({size_kb:.0f} KB; fetches siblings over HTTP)")
    except OSError as e:
        print(f"  WARNING: failed to write PeptideMapper.html: {e}")


def _load_fasta_dict(fasta_path: str) -> Dict[str, str]:
    """Parse FASTA → {accession: uppercase_sequence}.

    Accession extraction matches proteoform_analysis convention: for
    UniProt-style ``sp|P12345|NAME_HUMAN`` use the middle token; otherwise
    the first whitespace-delimited token of the header.
    """
    out: Dict[str, str] = {}
    cur_acc, cur = None, []
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_acc is not None:
                    out[cur_acc] = "".join(cur).upper()
                header = line[1:]
                if "|" in header:
                    parts = header.split("|")
                    cur_acc = parts[1] if len(parts) >= 2 else parts[0]
                else:
                    cur_acc = header.split()[0]
                cur = []
            else:
                cur.append(line.strip())
        if cur_acc is not None:
            out[cur_acc] = "".join(cur).upper()
    return out


def _add_sequence_positions(rpt: pd.DataFrame, fasta_seqs: Dict[str, str]) -> None:
    """Add Sequence.Start / Sequence.Stop columns in-place (1-based).

    For multi-protein groups in ``Protein.Ids`` (semicolon-separated), emits
    a semicolon-aligned list of positions in the same order. Missing parent
    or peptide-not-found yields an empty token in that slot.
    """
    pid_col = ("Protein.Ids"   if "Protein.Ids"   in rpt.columns else
               "Protein.Group" if "Protein.Group" in rpt.columns else None)
    if pid_col is None or "Stripped.Sequence" not in rpt.columns:
        print("  WARNING: cannot add Sequence.Start/Stop — missing Protein.Ids or Stripped.Sequence")
        return

    cache: Dict[Tuple[str, str], Tuple[Optional[int], Optional[int]]] = {}
    def lookup(acc: str, pep: str) -> Tuple[Optional[int], Optional[int]]:
        key = (acc, pep)
        v = cache.get(key)
        if v is not None:
            return v
        seq = fasta_seqs.get(acc)
        if seq is None:
            v = (None, None)
        else:
            i = seq.find(pep)
            v = (None, None) if i < 0 else (i + 1, i + len(pep))
        cache[key] = v
        return v

    n_rows, n_unmapped = 0, 0
    starts, stops = [], []
    for pids, pep in zip(rpt[pid_col].astype(str), rpt["Stripped.Sequence"].astype(str)):
        n_rows += 1
        accs = [a.strip() for a in pids.split(";")]
        ss, tt, any_hit = [], [], False
        for a in accs:
            s, t = lookup(a, pep)
            if s is None:
                ss.append(""); tt.append("")
            else:
                ss.append(str(s)); tt.append(str(t)); any_hit = True
        if not any_hit:
            n_unmapped += 1
        starts.append(";".join(ss))
        stops.append(";".join(tt))
    rpt["Sequence.Start"] = starts
    rpt["Sequence.Stop"]  = stops
    if n_unmapped:
        print(f"  WARNING: {n_unmapped:,}/{n_rows:,} rows had no peptide-in-protein hit")


def annotate_report(report_path: str, lib_path: str, out_dir: str,
                    default_specificity: Optional[str] = None,
                    fasta_path: Optional[str] = None):
    """Extend a DIA-NN report with Enzyme + Specificity from the combined library.

    When ``fasta_path`` is given, also adds Sequence.Start / Sequence.Stop
    (1-based positions of Stripped.Sequence within each parent protein in
    Protein.Ids, semicolon-aligned for multi-protein groups).
    """
    if not (os.path.exists(report_path) and os.path.exists(lib_path)):
        return
    lib    = pd.read_csv(lib_path, sep="\t", low_memory=False,
                         usecols=lambda c: c in ("transition_group_id", "Enzyme", "Specificity"))
    if "Specificity" not in lib.columns:
        lib["Specificity"] = default_specificity
    lookup = (lib.drop_duplicates("transition_group_id")
                 .set_index("transition_group_id")[["Enzyme", "Specificity"]])
    rpt    = pd.read_csv(report_path, sep="\t", low_memory=False)
    id_col = "Precursor.Id" if "Precursor.Id" in rpt.columns else "transition_group_id"
    rpt["Enzyme"]      = rpt[id_col].map(lookup["Enzyme"])
    rpt["Specificity"] = rpt[id_col].map(lookup["Specificity"])
    unmapped = rpt["Enzyme"].isna().sum()
    if unmapped:
        print(f"  WARNING: {unmapped:,} precursors unmapped in annotation")
    if fasta_path and os.path.exists(fasta_path):
        fasta_seqs = _load_fasta_dict(fasta_path)
        print(f"  FASTA: {len(fasta_seqs):,} entries — annotating Sequence.Start/Stop")
        _add_sequence_positions(rpt, fasta_seqs)
    out_path = os.path.join(out_dir, "report_annotated.tsv")
    write_tsv(rpt, out_path)
    print(f"  Annotated report: {os.path.basename(out_path)}")

    q_col  = "Q.Value"         if "Q.Value"         in rpt.columns else None
    pq_col = ("Protein.Q.Value" if "Protein.Q.Value" in rpt.columns
              else "PG.Q.Value" if "PG.Q.Value"      in rpt.columns else None)
    if q_col:
        flt = rpt[rpt[q_col] <= 0.01]
        if pq_col:
            flt = flt[flt[pq_col] <= 0.01]
        if "Enzyme" in flt.columns and "Specificity" in flt.columns:
            prec_col = "Precursor.Id" if "Precursor.Id" in flt.columns else "transition_group_id"
            summary  = (flt.drop_duplicates(prec_col)
                           .groupby(["Enzyme", "Specificity"])[prec_col].count()
                           .rename("N").reset_index()
                           .sort_values(["Enzyme", "Specificity"]))
            print(f"\n  Precursors at 1% FDR by Enzyme × Specificity:")
            for _, row in summary.iterrows():
                print(f"    {row['Enzyme']:<15}  {row['Specificity']:<18}  {row['N']:>6,}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Multi-enzyme DIA Proteomics Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--raw",              action="append", nargs="+",
                    required=True,
                    help="Input RAW file(s). "
                         "Accepts both --raw f1 f2 f3 and --raw f1 --raw f2 forms.")
    ap.add_argument("--fasta",            default=None,
                    help="FASTA protein database (required for pass1 prediction and annotation)")
    ap.add_argument("--out-dir",          required=True)
    ap.add_argument("--enzymes",          default="TCGE",
                    help=f"Enzyme codes, e.g. TCGE or T,C,G,E. "
                         f"Known: {','.join(ENZYME_REGISTRY)}")
    ap.add_argument("--threads",          type=int, default=14)
    ap.add_argument("--semi-ext",         action="store_true",
                    help="Enable semi-specific extension (pass3 + pass4)")
    ap.add_argument("--no-isd-filter",    action="store_true",
                    help="Proteoform step: keep in-source-decay products (issue #42)")
    ap.add_argument("--isd-rt-window",    type=float, default=0.05,
                    help="Proteoform step: max apex RT difference (min) to a longer, more "
                         "intense parent for an in-source-decay flag (default 0.05)")
    ap.add_argument("--speclib",          action="append", default=[],
                    metavar="CODE=PATH",
                    help="Per-enzyme predicted speclib, e.g. --speclib T=/path/T.speclib "
                         "(can be repeated for multiple enzymes)")
    ap.add_argument("--speclib-dir",      default=None,
                    help="Fallback: scan this directory for .predicted.speclib files")
    ap.add_argument("--lib",              default=None,
                    help="Skip pass1; use this combined library TSV for pass2")
    ap.add_argument("--skip-pass1",       action="store_true")
    ap.add_argument("--diann",            default=os.environ.get("MEDIA_DIANN_EXE", "diann"),
                    help="DIA-NN executable path (default: $MEDIA_DIANN_EXE, "
                         "else 'diann' from PATH)")
    ap.add_argument("--min-pep-len",      type=int, default=7)
    ap.add_argument("--max-pep-len",      type=int, default=30)
    ap.add_argument("--missed-cleavages", type=int, default=1,
                    help="missed cleavages for enzymes without their own value "
                         "(default 1); see --enzyme-mc")
    ap.add_argument("--enzyme-mc",        action="append", default=[], metavar="CODE=N",
                    help="missed cleavages for one enzyme, e.g. --enzyme-mc EX=2 "
                         "(repeatable). Precedence: this, then the enzyme's own "
                         "default (EX: 3), then --missed-cleavages")
    ap.add_argument("--enzyme-config",    default=None,
                    help="JSON file extending the enzyme registry")
    mg = ap.add_argument_group(
        "modifications",
        "Applied identically to every DIA-NN call and recorded in pipeline.settings. "
        "Default: carbamidomethyl (C) fixed, for alkylated samples.")
    mg.add_argument("--mods",          default=",".join(DEFAULT_MOD_PRESETS),
                    help=f"comma-separated presets, or 'none' "
                         f"(known: {', '.join(MOD_PRESETS)}; default: %(default)s)")
    mg.add_argument("--fixed-mod",     action="append", default=[], metavar="NAME,MASS,SITES",
                    help="additional fixed modification, DIA-NN syntax, "
                         "e.g. UniMod:4,57.021464,C (repeatable)")
    mg.add_argument("--var-mod",       action="append", default=[], metavar="NAME,MASS,SITES",
                    help="additional variable modification; sites may include "
                         "n (peptide N-term) and *n (protein N-term) (repeatable)")
    mg.add_argument("--var-mods",      type=int, default=1, metavar="N",
                    help="maximum variable modifications per peptide (default 1)")
    mg.add_argument("--monitor-mod",   action="append", default=[], metavar="NAME",
                    help="treat this variable modification as a PTM of interest; adds "
                         "DIA-NN --monitor-mod, which is what produces the PTM.Q.Value "
                         "and PTM.Site.Confidence columns (repeatable)")
    mg.add_argument("--no-cut-after-mod", action="append", default=[], metavar="NAME",
                    help="the modified residue blocks cleavage, e.g. GlyGly on K against "
                         "trypsin; adds DIA-NN --no-cut-after-mod (repeatable)")
    mg.add_argument("--met-excision",  action="store_true",
                    help="variable protein N-term methionine excision")
    mg.add_argument("--diann-arg",     action="append", default=[], metavar="ARG",
                    help="escape hatch: raw argument appended to every DIA-NN call, "
                         "write as --diann-arg=--flag (repeatable)")
    return ap


_CLI_PATH_ARGS = ("fasta", "out_dir", "speclib_dir", "lib", "enzyme_config", "diann")


def normalise_cli_paths(args, to_posix=None):
    """Windows-form path arguments -> POSIX (call under WSL only).

    Everything downstream then sees one convention; in particular the library
    guard (settings_for_lib) walks a library's folders with os.path, which a
    Windows-form path would defeat.
    """
    conv = to_posix or wslpaths.to_posix
    args.raw = [[conv(r) for r in group] for group in args.raw]
    for name in _CLI_PATH_ARGS:
        if getattr(args, name, None):
            setattr(args, name, conv(getattr(args, name)))
    args.speclib = [f"{code}={conv(path)}" if sep else kv
                    for kv in args.speclib
                    for code, sep, path in [kv.partition("=")]]
    return args


def parse_enzyme_codes(text: str) -> Tuple[List[str], List[str]]:
    """'TCGEX' or 'T,C,G,EX' -> (['T','C','G','EX'], unknown). Case-insensitive;
    a 2-character code is tried before a 1-character one, so EX is not E + X."""
    raw = text.replace(",", "").upper()
    codes, unknown, i = [], [], 0
    while i < len(raw):
        if raw[i:i+2] in ENZYME_REGISTRY:
            codes.append(raw[i:i+2]); i += 2
        elif raw[i] in ENZYME_REGISTRY:
            codes.append(raw[i]);     i += 1
        else:
            unknown.append(raw[i]);   i += 1
    return codes, unknown


def main(argv: Optional[List[str]] = None):
    ap   = build_arg_parser()
    args = ap.parse_args(argv)
    PF_ISD["enabled"] = not args.no_isd_filter
    PF_ISD["rtWindowMin"] = args.isd_rt_window
    if wslpaths.is_wsl():
        try:
            normalise_cli_paths(args)
        except ValueError as e:
            ap.error(str(e))
    try:
        mods = mod_config_from_args(args)
    except ModConfigError as e:
        ap.error(str(e))

    # Flatten --raw (action=append,nargs="+") → each --raw appends a sublist
    raw_files = [os.path.abspath(r)
                 for group in args.raw for r in group]
    fasta       = os.path.abspath(args.fasta) if args.fasta else None
    out_dir     = os.path.abspath(args.out_dir)
    speclib_dir = os.path.abspath(args.speclib_dir) if args.speclib_dir else None
    os.makedirs(out_dir, exist_ok=True)

    # Per-enzyme speclib dict: parse CODE=PATH pairs
    per_enzyme_speclibs: Dict[str, str] = {}
    for kv in args.speclib:
        if "=" in kv:
            code, _, path = kv.partition("=")
            per_enzyme_speclibs[code.strip().upper()] = path.strip()

    # Extend enzyme registry if requested
    if args.enzyme_config:
        with open(args.enzyme_config) as f:
            for code, d in json.load(f).items():
                ENZYME_REGISTRY[code] = EnzymeDef(**d)

    # Parse enzyme codes
    enzyme_codes, unknown = parse_enzyme_codes(args.enzymes)
    for u in unknown:
        print(f"  WARNING: unknown enzyme code '{u}' — skipped")

    if not enzyme_codes:
        ap.error(f"No valid enzyme codes. Known: {','.join(ENZYME_REGISTRY)}")
    codes_str = "".join(enzyme_codes)

    # Effective missed cleavages per enzyme, resolved once and passed down.
    try:
        mc_overrides = parse_enzyme_mc(args.enzyme_mc)
    except ValueError as e:
        ap.error(str(e))
    for code in mc_overrides:
        if code not in enzyme_codes:
            print(f"  WARNING: --enzyme-mc {code}=… is for an enzyme not in this run — ignored")
    mc_overrides = {c: n for c, n in mc_overrides.items() if c in enzyme_codes}
    enzyme_mc    = resolve_enzyme_mc(enzyme_codes, args.missed_cleavages, mc_overrides)
    enzyme_specs = [enzyme_spec(c, enzyme_mc[c]) for c in enzyme_codes]

    log_path      = os.path.join(out_dir, "pipeline.log")
    settings_path = os.path.join(out_dir, "pipeline.settings")
    timing:  Dict[str, float] = {}

    # Resuming into a folder whose outputs were made with other modifications
    # would mix configurations. Settings without the key predate modification
    # support, when no modification flag was ever emitted.
    try:
        check_resume_mods(out_dir, mods)
        check_resume_enzymes(out_dir, enzyme_specs)
    except ModConfigError as e:
        ap.error(str(e))
    recover_interrupted(out_dir)

    # ── Save pipeline settings for potential resume ───────────────────────────
    settings_data = {
        "timestamp":       datetime.now().isoformat(),
        "rawFiles":        raw_files,
        "fasta":           fasta or "",
        "outDir":          out_dir,
        "enzymes":         enzyme_codes,
        "speclibPaths":    per_enzyme_speclibs,
        "semiExt":         args.semi_ext,
        "proteoformIsdFilter": {"enabled": not args.no_isd_filter,
                                "rtWindowMin": args.isd_rt_window},
        "skipPass1":       args.skip_pass1,
        "combinedLib":     args.lib or "",
        "minPepLen":       args.min_pep_len,
        "maxPepLen":       args.max_pep_len,
        "missedCleavages": args.missed_cleavages,
        "enzymeMissedCleavages": mc_overrides,
        "enzymeDefinitions":     enzyme_specs,
        **mods,
        "diannModFlags":   mod_flags(mods),
        "buildInfo":       BUILD_INFO,
        "runtime":         RUNTIME_INFO(args.diann),
        "diannExe":        args.diann,
        "threads":         args.threads,
    }
    with open(settings_path, "w", encoding="utf-8") as _sf:
        json.dump(settings_data, _sf, indent=2)
    print(f"  Settings: {settings_path}", flush=True)

    print(f"\n{'='*70}")
    print(f"  Multi-enzyme DIA Pipeline — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  RAW     : {', '.join(os.path.basename(r) for r in raw_files)}")
    if fasta:
        print(f"  FASTA   : {os.path.basename(fasta)}")
    print(f"  OUT     : {out_dir}")
    print(f"  Enzymes : {codes_str} "
          f"({', '.join(ENZYME_REGISTRY[c].name for c in enzyme_codes)})")
    print(f"  Missed  : "
          + ", ".join(f"{c}={enzyme_mc[c]}" for c in enzyme_codes)
          + f"  (default {args.missed_cleavages})")
    print(f"  Threads : {args.threads}")
    print(f"  SemiExt : {'yes (pass3+pass4)' if args.semi_ext else 'no'}")
    print(f"  Build   : {_build_info.one_line(BUILD_INFO) if _build_info else 'unknown'}")
    print(f"  Runtime : {_build_info.runtime_one_line(settings_data['runtime']) if _build_info else 'unknown'}")
    print(f"  DIA-NN  : {args.diann}  (version: from its banner at the first DIA-NN call)")
    print(f"  Mods    : {describe_mods(mod_flags(mods))}")
    if not any(m["sites"] == "C" and m["type"] == "fixed" for m in mods["modifications"]):
        print("  !! WARNING: no fixed cysteine modification. Correct only for "
              "non-alkylated samples; alkylated cysteine peptides will be missed.")
    if per_enzyme_speclibs:
        for code, p in per_enzyme_speclibs.items():
            print(f"  Speclib {code}: {os.path.basename(p)}")
    if speclib_dir:
        print(f"  SpeclibDir: {speclib_dir}")
    print(f"{'='*70}\n")

    _emit("STATUS", {"pass": "pipeline", "status": "running",
                      "codes": codes_str, "semi_ext": args.semi_ext})

    # ── Pass 1 ────────────────────────────────────────────────────────────────
    combined_lib:    Optional[str]                  = None
    pass2_params:    dict                           = {}
    pass1_results:   Dict[str, Tuple[str, dict]]   = {}
    pass1_params_map: Dict[str, dict]               = {}

    if args.lib:
        combined_lib = os.path.abspath(args.lib)
        print(f"  Using provided combined library: {combined_lib}")
        try:
            check_lib_mods(combined_lib, mod_flags(mods, include_extra=False),
                           "Provided combined library (--lib)", enzymes=enzyme_specs)
        except ModConfigError as e:
            ap.error(str(e))
    elif not args.skip_pass1:
        if not fasta:
            ap.error("--fasta is required for pass1 prediction")
        print("─" * 70)
        print("  PASS 1 — per-enzyme predict + search")
        print("─" * 70)
        pass1_results = run_pass1(
            raw_files, out_dir, enzyme_codes, fasta,
            per_enzyme_speclibs, speclib_dir,
            args.threads, args.diann, log_path, timing,
            args.min_pep_len, args.max_pep_len, args.missed_cleavages,
            mods=mods, enzyme_mc=enzyme_mc,
        )
        pass1_params_map = {code: res[1] for code, res in pass1_results.items()}

        p2_dir       = os.path.join(out_dir, f"pass2_{codes_str}")
        os.makedirs(p2_dir, exist_ok=True)
        combined_lib = os.path.join(p2_dir, f"combined_{codes_str}.tsv")
        if not (os.path.exists(combined_lib) and os.path.getsize(combined_lib) > 100):
            print(f"\n[{ts()}] Combining empirical libraries...")
            emp_libs = {ENZYME_REGISTRY[c].name: res[0]
                        for c, res in pass1_results.items()
                        if os.path.exists(res[0]) and os.path.getsize(res[0]) > 100}
            combine_libraries(emp_libs, combined_lib)
        else:
            print(f"  [SKIP] Combined lib exists")

    # ── Pass 2 ────────────────────────────────────────────────────────────────
    if not combined_lib or not os.path.exists(combined_lib):
        print("WARNING: No combined library — skipping pass2+.")
        return
    if not fasta:
        ap.error("--fasta is required for pass2 matrix annotation")

    print("\n" + "─" * 70)
    print(f"  PASS 2 — combined {codes_str} library + matrices")
    print("─" * 70)
    report2, pass2_params = run_pass2(
        raw_files, out_dir, enzyme_codes, combined_lib, fasta,
        args.threads, args.diann, log_path, timing, mods=mods,
    )
    print_id_summary(report2, f"Pass2 {codes_str}")
    emit_chart_data(report2, combined_lib, f"pass2_{codes_str}")
    p2_dir = os.path.dirname(report2)
    annotate_report(report2, combined_lib,
                    p2_dir, default_specificity="fully_specific",
                    fasta_path=fasta)
    fasta_used_p2 = _copy_fasta_to_out(fasta, p2_dir)
    report_annotated_p2 = os.path.join(p2_dir, "report_annotated.tsv")
    # Proteoform analysis + viewer copy run only at the canonical final pass.
    # When --semi-ext is set, that's pass4; otherwise pass2 is final.
    if not args.semi_ext:
        _run_proteoform(report_annotated_p2, fasta_used_p2, p2_dir)
        _copy_peptide_mapper(p2_dir, report_annotated_p2, fasta_used_p2)

    # ── Pass 3 + Pass 4 (optional semi-ext) ──────────────────────────────────
    if not args.semi_ext:
        print(f"\n  Semi-extension not requested.  Done.")
        _print_timing(timing)
        print(f"\n  Results: {out_dir}")
        _emit("STATUS", {"pass": "pipeline", "status": "done"})
        return

    # Protein whitelist from pass2 report
    protein_whitelist = None
    if report2 and os.path.exists(report2):
        try:
            peek   = pd.read_csv(report2, sep="\t", nrows=0)
            pg_col = ("Protein.Group" if "Protein.Group" in peek.columns else "ProteinGroup")
            pg_df  = pd.read_csv(report2, sep="\t", usecols=[pg_col], low_memory=False)
            protein_whitelist = set()
            for pg in pg_df[pg_col].dropna().astype(str):
                for acc in pg.split(";"):
                    acc = acc.strip()
                    if acc:
                        protein_whitelist.add(acc)
            print(f"\n  Protein whitelist: {len(protein_whitelist):,} accessions from pass2")
        except Exception as e:
            print(f"  WARNING: could not build whitelist: {e}")

    print("\n" + "─" * 70)
    print("  PASS 3 — per-enzyme semi-specific extension")
    print("           (inline reannotation; fixed params from pass1)")
    print("─" * 70)
    pass3_results = run_pass3(
        raw_files, out_dir, enzyme_codes,
        {code: (res[0], pass1_params_map.get(code, {}))
         for code, res in pass1_results.items()},
        fasta, args.threads, args.diann, log_path, timing,
        args.min_pep_len, args.max_pep_len, args.missed_cleavages,
        protein_whitelist, mods=mods, enzyme_mc=enzyme_mc,
    )

    # Build extended combined library (pass2 fully-specific + semi)
    p4_dir       = os.path.join(out_dir, f"pass4_semiExt_{codes_str}")
    os.makedirs(p4_dir, exist_ok=True)
    ext_combined = os.path.join(p4_dir, f"{codes_str}_empirical_extSemi.tsv")

    if not (os.path.exists(ext_combined) and os.path.getsize(ext_combined) > 100):
        print(f"\n[{ts()}] Building extended combined library (pass2 + semi)...")
        p1_df     = pd.read_csv(combined_lib, sep="\t", low_memory=False)
        if "Specificity" not in p1_df.columns:
            p1_df["Specificity"] = "fully_specific"
        pass1_ids = set(p1_df["transition_group_id"])

        semi_dfs = []
        for code, (emp_lib, vmap_df) in pass3_results.items():
            enz = ENZYME_REGISTRY[code]
            if not (os.path.exists(emp_lib) and os.path.getsize(emp_lib) > 100):
                continue
            df = pd.read_csv(emp_lib, sep="\t", low_memory=False)
            if len(df) == 0:
                continue
            df = _reannotate_protein_groups(df, vmap_df)
            if "Enzyme" not in df.columns:
                df.insert(0, "Enzyme", enz.name)
            if "Specificity" not in df.columns:
                if vmap_df is not None and "PeptideSequence" in df.columns:
                    spec_lup    = vmap_df.set_index("PeptideSequence")["Specificity"].to_dict()
                    df["Specificity"] = df["PeptideSequence"].map(spec_lup).fillna("semi_specific")
                else:
                    df["Specificity"] = "semi_specific"
            semi_dfs.append(df)

        if not semi_dfs:
            raise RuntimeError("No semi empirical libraries — pass3 may have failed.")
        semi_df  = pd.concat(semi_dfs, ignore_index=True)
        # Multiplex: a peptide semi-specific for the enzyme whose library found
        # it can be fully specific for another enzyme of the same experiment.
        # Only peptides fully specific under NO enzyme are semi_specific.
        fully_any = set()
        for code, (_lib, vmap_df) in pass3_results.items():
            if vmap_df is not None:
                fully_any |= set(vmap_df.loc[vmap_df["Specificity"] == "fully_specific",
                                             "PeptideSequence"])
        if "PeptideSequence" in semi_df.columns and fully_any:
            semi_df.loc[semi_df["PeptideSequence"].isin(fully_any), "Specificity"] = "fully_specific"
        sm_new   = semi_df[~semi_df["transition_group_id"].isin(pass1_ids)].copy()
        n_new    = sm_new["transition_group_id"].nunique()
        n_dup    = semi_df["transition_group_id"].nunique() - n_new
        print(f"  Pass2 precursors  : {len(pass1_ids):,}")
        print(f"  Semi new (unique) : {n_new:,}  (duplicates excluded: {n_dup:,})")
        merged   = pd.concat([p1_df, sm_new], ignore_index=True)
        write_tsv(merged, ext_combined)
        print(f"  Merged total      : {merged['transition_group_id'].nunique():,}"
              f"  →  {os.path.basename(ext_combined)}")
    else:
        print(f"  [SKIP] Extended combined lib exists")

    print("\n" + "─" * 70)
    print(f"  PASS 4 — semi-ext {codes_str} + matrices (fixed params from pass2)")
    print("─" * 70)
    report4 = run_pass4(
        raw_files, out_dir, enzyme_codes, ext_combined, pass2_params,
        fasta, args.threads, args.diann, log_path, timing, mods=mods,
    )
    print_id_summary(report4, f"Pass4 semi-ext {codes_str}")
    emit_chart_data(report4, ext_combined, f"pass4_semiExt_{codes_str}")
    annotate_report(report4, ext_combined, p4_dir, fasta_path=fasta)
    fasta_used_p4 = _copy_fasta_to_out(fasta, p4_dir)
    report_annotated_p4 = os.path.join(p4_dir, "report_annotated.tsv")
    _run_proteoform(report_annotated_p4, fasta_used_p4, p4_dir)
    _copy_peptide_mapper(p4_dir, report_annotated_p4, fasta_used_p4)

    _print_timing(timing)
    print(f"\n  Results: {out_dir}")
    print(f"  All done — {ts()}")
    _emit("STATUS", {"pass": "pipeline", "status": "done"})


def _print_timing(timing: dict):
    if not timing:
        return
    total = sum(timing.values())
    print(f"\n{'─'*55}")
    print("  TIMING SUMMARY")
    print(f"{'─'*55}")
    for step, elapsed in sorted(timing.items()):
        print(f"  {step:<38} {elapsed / 60:>5.1f} min")
    print(f"  {'TOTAL':<38} {total / 60:>5.1f} min")


if __name__ == "__main__":
    main()
