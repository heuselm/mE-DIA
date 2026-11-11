export interface BuildInfo {
  commit:      string
  branch:      string
  dirty?:      boolean
  built_utc?:  string
  describe?:   string
  line?:       string
}

export interface Config {
  diannExe: string
  threads:  number
}

// One DIA-NN modification: --fixed-mod / --var-mod NAME,MASS,SITES (DIA-NN 1.8.1).
// sites: amino-acid letters, 'n' (peptide N-term), '*n' (protein N-term).
export interface Modification {
  name:  string
  mass:  number
  sites: string
  type:  'fixed' | 'variable'
  // PTM of interest: adds DIA-NN --monitor-mod, which is what produces the
  // PTM.Q.Value / PTM.Site.Confidence columns. Variable modifications only.
  monitor?:    boolean
  // the modified residue blocks cleavage (DIA-NN --no-cut-after-mod), e.g. GlyGly on K
  noCutAfter?: boolean
}

// Mirrors MOD_PRESETS in bin/run_pipeline.py (monoisotopic, UniMod compositions).
export const MOD_PRESETS: { key: string; label: string; mod: Modification }[] = [
  { key: 'carbamidomethyl', label: 'Carbamidomethyl (C), fixed',
    mod: { name: 'UniMod:4',  mass: 57.021464, sites: 'C',  type: 'fixed' } },
  { key: 'oxidation',       label: 'Oxidation (M), variable',
    mod: { name: 'UniMod:35', mass: 15.994915, sites: 'M',  type: 'variable' } },
  { key: 'acetyl-nterm',    label: 'Acetyl (protein N-term), variable',
    mod: { name: 'UniMod:1',  mass: 42.010565, sites: '*n', type: 'variable',
           monitor: true } },
  { key: 'phospho',         label: 'Phospho (STY), variable, localised',
    mod: { name: 'UniMod:21', mass: 79.966331, sites: 'STY', type: 'variable',
           monitor: true } },
  { key: 'gg',              label: 'GlyGly (K), variable, localised, blocks cleavage',
    mod: { name: 'UniMod:121', mass: 114.042927, sites: 'K', type: 'variable',
           monitor: true, noCutAfter: true } },
]

export interface JobConfig {
  rawFiles:        string[]
  fasta:           string
  outDir:          string
  enzymes:         string[]              // 1- or 2-letter codes
  speclibPaths:    Record<string, string> // CODE → path to predicted speclib
  semiExt:         boolean
  skipPass1:       boolean
  combinedLib:     string
  minPepLen:       number
  maxPepLen:       number
  missedCleavages: number              // default for enzymes without their own
  // Per-enzyme overrides only (CODE → N), sent as --enzyme-mc CODE=N. An enzyme
  // without one uses its own default (ENZYME_DEFS mc) or missedCleavages.
  enzymeMissedCleavages: Record<string, number>
  modifications:   Modification[]   // applied identically to every DIA-NN pass
  maxVarMods:      number
  metExcision:     boolean
  extraDiannArgs:  string[]         // raw escape hatch, appended to every DIA-NN call
}

export type PassStatus = 'pending' | 'running' | 'done' | 'failed' | 'skipped' | 'cancelled'

export interface PassState {
  pass:      string
  status:    PassStatus
  elapsedS?: number
  skipped?:  boolean
  window?:   number
  ms1_acc?:  number
  mass_acc?: number
}

export interface Metric {
  label:      string
  proteins:   number | string
  peptides:   number | string
  precursors: number | string
}

export interface ChartDataPerRun {
  pass:       string
  type:       'per_run'
  runs:       string[]
  precursors: number[]
  peptides:   number[]   // Modified.Sequence (falls back to Stripped.Sequence)
  proteins?:  number[]   // unique Protein.Group per run
}

export interface ChartDataEnzyme {
  pass:    string
  type:    'enzyme_breakdown'
  enzymes: string[]
  counts:  number[]
}

export type ChartData = ChartDataPerRun | ChartDataEnzyme

export interface RunState {
  runId:        string
  status:       'running' | 'done' | 'failed' | 'cancelled'
  outputLines:  string[]
  passStatuses: Record<string, PassState>
  metrics:      Metric[]
  chartData:    ChartData[]
}

// Mirrors ENZYME_REGISTRY in bin/run_pipeline.py. mc: the enzyme's own
// missed-cleavage default, where it has one.
export const ENZYME_DEFS: Record<string, { name: string; cut: string; mc?: number }> = {
  T:  { name: 'Trypsin',       cut: 'K*,R*,!*P'       },
  T2: { name: 'Trypsin/P',     cut: 'K*,R*'            },
  C:  { name: 'Chymotrypsin',  cut: 'F*,W*,Y*,L*,!*P' },
  G:  { name: 'GluC',          cut: 'D*,E*,!*P'        },
  E:  { name: 'Elastase (AV)', cut: 'A*,V*'            },
  EX: { name: 'Elastase ext. (AVSGLI)', cut: 'A*,V*,S*,G*,L*,I*', mc: 3 },
  L:  { name: 'LysC',          cut: 'K*,!*P'           },
  LN: { name: 'LysN',          cut: '*K'               },
  R:  { name: 'ArgC',          cut: 'R*,!*P'           },
  A:  { name: 'Asp-N',         cut: '*D'               },
  M:  { name: 'CNBr',          cut: 'M*'               },
  F:  { name: 'PepsinA',       cut: 'F*,L*,!*P'        },
}

export const DEFAULT_CONFIG: Config = {
  diannExe: 'diann',
  threads:  14,
}

export const DEFAULT_JOB: JobConfig = {
  rawFiles:        [],
  fasta:           '',
  outDir:          '',
  enzymes:         ['T', 'C', 'G', 'E'],
  speclibPaths:    {},
  semiExt:         false,
  skipPass1:       false,
  combinedLib:     '',
  minPepLen:       7,
  maxPepLen:       30,
  missedCleavages: 1,
  enzymeMissedCleavages: {},
  // Alkylated samples are the norm: carbamidomethyl (C) on by default, visible
  // in the Modifications panel and recorded in pipeline.settings.
  modifications:   [MOD_PRESETS[0].mod],
  maxVarMods:      1,
  metExcision:     false,
  extraDiannArgs:  [],
}

// The missed cleavages an enzyme runs with: same precedence as the pipeline's
// resolve_enzyme_mc (override, then the enzyme's default, then the global value).
export function effectiveMissedCleavages(
  code: string, job: Pick<JobConfig, 'missedCleavages' | 'enzymeMissedCleavages'>,
): number {
  return job.enzymeMissedCleavages?.[code] ?? ENZYME_DEFS[code]?.mc ?? job.missedCleavages
}
