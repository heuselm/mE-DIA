import { useState, useEffect, useRef, useCallback } from 'react'
import {
  Config, JobConfig, RunState, Metric, PassState, ChartData, ChartDataPerRun,
  BuildInfo, ENZYME_DEFS, DEFAULT_CONFIG, DEFAULT_JOB, effectiveMissedCleavages,
  Modification, MOD_PRESETS,
} from './types'

// ─── API helpers ─────────────────────────────────────────────────────────────

const API = (path: string) => `/api${path}`

async function apiFetch(path: string, method = 'GET', body?: object) {
  const r = await fetch(API(path), {
    method,
    headers: body ? { 'Content-Type': 'application/json' } : undefined,
    body: body ? JSON.stringify(body) : undefined,
  })
  if (!r.ok) throw new Error(await r.text())
  return r.json()
}

// A dialog that cannot open (no display, WSL without interop, timeout) comes
// back as an HTTP error with a reason; show it rather than silently doing nothing.
async function browse(endpoint: string): Promise<string | string[]> {
  try {
    const d = await apiFetch(endpoint)
    return d.path ?? d.paths ?? ''
  } catch (e: unknown) {
    const raw = e instanceof Error ? e.message : String(e)
    let msg = raw
    try { msg = JSON.parse(raw).detail ?? raw } catch { /* plain text */ }
    alert(`Browse: ${msg}`)
    return ''
  }
}

// ─── Sub-components ──────────────────────────────────────────────────────────

function PathInput({
  label, value, onChange, browseEp, placeholder, mono = false,
}: {
  label: string; value: string; onChange: (v: string) => void
  browseEp: string; placeholder?: string; mono?: boolean
}) {
  const [busy, setBusy] = useState(false)
  const doBrowse = async () => {
    setBusy(true)
    try {
      const r = await browse(browseEp)
      if (typeof r === 'string' && r) onChange(r)
    } finally { setBusy(false) }
  }
  return (
    <div className="mb-3">
      <label className="block text-xs font-medium text-gray-400 mb-1">{label}</label>
      <div className="flex gap-1">
        <input
          value={value}
          onChange={e => onChange(e.target.value)}
          placeholder={placeholder}
          className={`flex-1 bg-gray-900 border border-gray-700 rounded px-2 py-1.5 text-sm
            text-gray-100 placeholder-gray-600 focus:outline-none focus:border-blue-500
            ${mono ? 'font-mono text-xs' : ''}`}
        />
        <button
          onClick={doBrowse}
          disabled={busy}
          className="px-2 py-1 bg-gray-700 hover:bg-gray-600 rounded text-xs text-gray-300
            disabled:opacity-50 whitespace-nowrap"
        >
          {busy ? '…' : 'Browse'}
        </button>
      </div>
    </div>
  )
}

function RawFilesInput({
  files, onChange,
}: { files: string[]; onChange: (v: string[]) => void }) {
  const [busy, setBusy] = useState(false)

  const addFiles = async () => {
    setBusy(true)
    try {
      const r = await browse('/browse/files')
      const added = Array.isArray(r) ? r : [r]
      const merged = [...files, ...added.filter(p => p && !files.includes(p))]
      onChange(merged)
    } finally { setBusy(false) }
  }

  const remove = (idx: number) => onChange(files.filter((_, i) => i !== idx))

  return (
    <div className="mb-3">
      <div className="flex items-center justify-between mb-1">
        <label className="text-xs font-medium text-gray-400">RAW Files</label>
        <button
          onClick={addFiles}
          disabled={busy}
          className="text-xs px-2 py-0.5 bg-blue-700 hover:bg-blue-600 rounded text-white
            disabled:opacity-50"
        >
          {busy ? '…' : '+ Add files'}
        </button>
      </div>
      {files.length === 0 && (
        <p className="text-xs text-gray-600 italic">No files selected</p>
      )}
      <div className="space-y-1 max-h-36 overflow-y-auto">
        {files.map((f, i) => (
          <div key={i} className="flex items-center gap-1 bg-gray-900 rounded px-2 py-1">
            <span className="flex-1 text-xs font-mono text-gray-300 truncate" title={f}>
              {f.split(/[\\/]/).pop()}
            </span>
            <button
              onClick={() => remove(i)}
              className="text-gray-600 hover:text-red-400 text-xs ml-1 flex-shrink-0"
            >✕</button>
          </div>
        ))}
      </div>
    </div>
  )
}

function EnzymeSelector({
  selected, onChange,
}: { selected: string[]; onChange: (v: string[]) => void }) {
  const toggle = (code: string) => {
    if (selected.includes(code)) onChange(selected.filter(c => c !== code))
    else onChange([...selected, code])
  }
  return (
    <div className="mb-3">
      <label className="block text-xs font-medium text-gray-400 mb-1">Enzymes</label>
      <div className="flex flex-wrap gap-1">
        {Object.entries(ENZYME_DEFS).map(([code, def]) => {
          const active = selected.includes(code)
          return (
            <button
              key={code}
              onClick={() => toggle(code)}
              title={`${def.name}  (${def.cut})`}
              className={`px-2 py-1 rounded text-xs font-medium transition-colors
                ${active
                  ? 'bg-blue-600 text-white'
                  : 'bg-gray-800 text-gray-400 hover:bg-gray-700'}`}
            >
              {code} – {def.name}
            </button>
          )
        })}
      </div>
    </div>
  )
}

// Missed cleavages per selected enzyme: the effective value, editable as an
// override; clearing a field returns the enzyme to its default. The line below
// shows what each pass-1 prediction will receive.
function MissedCleavagesPanel({ job, setJob }: {
  job: JobConfig; setJob: React.Dispatch<React.SetStateAction<JobConfig>>
}) {
  if (job.enzymes.length === 0) return null
  const setOverride = (code: string, raw: string) =>
    setJob(j => {
      const next = { ...(j.enzymeMissedCleavages ?? {}) }
      const n = parseInt(raw)
      if (raw.trim() === '' || Number.isNaN(n) || n < 0) delete next[code]
      else next[code] = n
      return { ...j, enzymeMissedCleavages: next }
    })
  return (
    <div className="mb-3">
      <label className="block text-xs font-medium text-gray-400 mb-1">
        Missed cleavages per enzyme
        <span className="text-gray-600 font-normal"> (default {job.missedCleavages}, set under Advanced)</span>
      </label>
      <div className="flex flex-wrap gap-2">
        {job.enzymes.map(code => {
          const overridden = job.enzymeMissedCleavages?.[code] !== undefined
          const own        = ENZYME_DEFS[code]?.mc
          return (
            <label key={code} className="flex items-center gap-1 text-xs text-gray-400"
              title={overridden ? 'override; clear to use the default'
                     : own !== undefined ? `${code}'s own default` : 'global default'}>
              <span className="font-mono">{code}</span>
              <input
                type="number" min={0} max={9}
                value={effectiveMissedCleavages(code, job)}
                onChange={e => setOverride(code, e.target.value)}
                className={`w-12 bg-gray-900 border rounded px-1 py-0.5 text-xs text-gray-100
                  focus:outline-none focus:border-blue-500
                  ${overridden ? 'border-blue-500' : 'border-gray-700'}`}
              />
              {overridden && (
                <button onClick={() => setOverride(code, '')}
                  className="text-gray-500 hover:text-gray-300" title="reset to default">↺</button>
              )}
            </label>
          )
        })}
      </div>
      <div className="mt-1 text-xs text-gray-500 font-mono break-all">
        {job.enzymes.map(code =>
          `${code}: --cut ${ENZYME_DEFS[code]?.cut ?? '?'} --missed-cleavages ${effectiveMissedCleavages(code, job)}`
        ).join('  ·  ')}
      </div>
    </div>
  )
}

function SpeclibInputs({
  enzymes, paths, onChange,
}: {
  enzymes: string[]
  paths: Record<string, string>
  onChange: (v: Record<string, string>) => void
}) {
  const [busy, setBusy] = useState<string | null>(null)

  const doBrowse = async (code: string) => {
    setBusy(code)
    try {
      const p = await browse(`/browse/file?title=Select+speclib+for+${code}`)
      if (typeof p === 'string' && p) onChange({ ...paths, [code]: p })
    } finally { setBusy(null) }
  }

  if (enzymes.length === 0) return null

  return (
    <div className="mb-3">
      <label className="block text-xs font-medium text-gray-400 mb-1">
        Per-enzyme predicted speclib <span className="text-gray-600">(optional)</span>
      </label>
      <div className="space-y-1.5">
        {enzymes.map(code => (
          <div key={code} className="flex gap-1 items-center">
            <span className="text-xs font-mono text-blue-400 w-6 flex-shrink-0">{code}</span>
            <input
              value={paths[code] ?? ''}
              onChange={e => onChange({ ...paths, [code]: e.target.value })}
              placeholder={`${ENZYME_DEFS[code]?.name ?? code}.predicted.speclib`}
              className="flex-1 bg-gray-900 border border-gray-700 rounded px-2 py-1
                text-xs font-mono text-gray-100 placeholder-gray-700
                focus:outline-none focus:border-blue-500"
            />
            <button
              onClick={() => doBrowse(code)}
              disabled={busy === code}
              className="px-2 py-1 bg-gray-700 hover:bg-gray-600 rounded text-xs
                text-gray-300 disabled:opacity-50 whitespace-nowrap"
            >
              {busy === code ? '…' : 'Browse'}
            </button>
          </div>
        ))}
      </div>
    </div>
  )
}

// ─── Simple CSS bar chart ─────────────────────────────────────────────────────

function SimpleBarChart({
  items, color = 'bg-blue-500', maxOverride,
}: { items: { label: string; value: number }[]; color?: string; maxOverride?: number }) {
  const max = maxOverride ?? Math.max(...items.map(d => d.value), 1)
  return (
    <div className="space-y-1.5">
      {items.map(d => (
        <div key={d.label} className="flex items-center gap-2 text-xs">
          <span className="w-40 text-right text-gray-400 truncate shrink-0" title={d.label}>
            {d.label}
          </span>
          <div className="flex-1 bg-gray-800 rounded overflow-hidden h-4 min-w-0">
            <div
              className={`h-full ${color} rounded transition-all duration-500`}
              style={{ width: `${(d.value / max) * 100}%` }}
            />
          </div>
          <span className="w-16 text-gray-300 text-right shrink-0">
            {d.value.toLocaleString()}
          </span>
        </div>
      ))}
    </div>
  )
}

// ─── Pipeline flow overview ────────────────────────────────────────────────────

type StepStatus = 'pending' | 'running' | 'done' | 'failed' | 'skipped' | 'cancelled'

function StepNode({
  label, status, elapsed,
}: { label: string; status: StepStatus; elapsed?: number }) {
  const ring: Record<StepStatus, string> = {
    pending:   'border-gray-700 text-gray-600',
    running:   'border-blue-500 text-blue-400',
    done:      'border-green-500 text-green-400',
    failed:    'border-red-500  text-red-400',
    skipped:   'border-gray-600 text-gray-500',
    cancelled: 'border-yellow-600 text-yellow-500',
  }
  const icon: Record<StepStatus, string> = {
    pending: '○', running: '↻', done: '✓', failed: '✗', skipped: '⏭', cancelled: '⊘',
  }
  const bg: Record<StepStatus, string> = {
    pending: 'bg-gray-900', running: 'bg-blue-950', done: 'bg-green-950',
    failed: 'bg-red-950', skipped: 'bg-gray-900', cancelled: 'bg-yellow-950',
  }
  return (
    <div className={`flex flex-col items-center gap-0.5 border rounded px-2 py-1.5
      min-w-[72px] ${ring[status]} ${bg[status]}`}>
      <span className={`text-base leading-none ${status === 'running' ? 'spinning' : ''}`}>
        {icon[status]}
      </span>
      <span className="text-xs font-mono text-center leading-tight">{label}</span>
      {elapsed !== undefined && (
        <span className="text-xs text-gray-500">{(elapsed / 60).toFixed(1)}m</span>
      )}
    </div>
  )
}

function PipelineFlow({
  enzymes, semiExt, passStatuses,
}: {
  enzymes:      string[]
  semiExt:      boolean
  passStatuses: Record<string, PassState>
}) {
  if (enzymes.length === 0) return null

  const codes = enzymes.join('')
  const getStatus = (key: string): StepStatus =>
    (passStatuses[key]?.status as StepStatus) ?? 'pending'
  const getElapsed = (key: string) => passStatuses[key]?.elapsedS

  return (
    <div className="px-4 py-3 border-b border-gray-800 bg-gray-900/50">
      <div className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">
        Pipeline Flow
      </div>
      <div className="flex flex-col gap-2">

        {/* Pass 1 */}
        <div className="flex items-center gap-1 flex-wrap">
          <span className="text-xs text-gray-500 w-14 shrink-0">Pass 1</span>
          <div className="flex gap-1 flex-wrap">
            {enzymes.map(code => (
              <StepNode
                key={code}
                label={`p1_${code}`}
                status={getStatus(`pass1_${code}`)}
                elapsed={getElapsed(`pass1_${code}`)}
              />
            ))}
          </div>
        </div>

        {/* Arrow */}
        <div className="ml-16 text-gray-700 text-xs">↓ combine libs</div>

        {/* Pass 2 */}
        <div className="flex items-center gap-1">
          <span className="text-xs text-gray-500 w-14 shrink-0">Pass 2</span>
          <StepNode
            label={`p2_${codes}`}
            status={getStatus(`pass2_${codes}`)}
            elapsed={getElapsed(`pass2_${codes}`)}
          />
        </div>

        {semiExt && (<>
          {/* Arrow */}
          <div className="ml-16 text-gray-700 text-xs">↓ semi-ext (virtual FASTA)</div>

          {/* Pass 3 */}
          <div className="flex items-center gap-1 flex-wrap">
            <span className="text-xs text-gray-500 w-14 shrink-0">Pass 3</span>
            <div className="flex gap-1 flex-wrap">
              {enzymes.map(code => (
                <StepNode
                  key={code}
                  label={`p3_${code}`}
                  status={getStatus(`pass3_semiExt_${code}`)}
                  elapsed={getElapsed(`pass3_semiExt_${code}`)}
                />
              ))}
            </div>
          </div>

          {/* Arrow */}
          <div className="ml-16 text-gray-700 text-xs">↓ combine semi libs</div>

          {/* Pass 4 */}
          <div className="flex items-center gap-1">
            <span className="text-xs text-gray-500 w-14 shrink-0">Pass 4</span>
            <StepNode
              label={`p4_${codes}`}
              status={getStatus(`pass4_semiExt_${codes}`)}
              elapsed={getElapsed(`pass4_semiExt_${codes}`)}
            />
          </div>
        </>)}
      </div>
    </div>
  )
}

// ─── Run charts (results tab) ─────────────────────────────────────────────────

function RunCharts({ chartData }: { chartData: ChartData[] }) {
  if (chartData.length === 0) return null

  // Shared scales across all per_run passes
  const perRunData = chartData.filter(cd => cd.type === 'per_run') as ChartDataPerRun[]
  const globalMaxPrecPep = Math.max(
    ...perRunData.flatMap(cd => [...cd.precursors, ...cd.peptides]),
    1,
  )
  const globalMaxProt = Math.max(
    ...perRunData.flatMap(cd => cd.proteins ?? []),
    1,
  )

  const ENZYME_COLORS = ['bg-blue-500','bg-emerald-500','bg-orange-500','bg-purple-500',
                         'bg-yellow-500','bg-pink-500','bg-cyan-500']

  return (
    <div className="space-y-6">
      {chartData.map((cd, i) => {
        if (cd.type === 'per_run') {
          const precItems = cd.runs.map((r, j) => ({ label: r, value: cd.precursors[j] ?? 0 }))
          const pepItems  = cd.runs.map((r, j) => ({ label: r, value: cd.peptides[j]   ?? 0 }))
          const protItems = cd.proteins
            ? cd.runs.map((r, j) => ({ label: r, value: (cd.proteins as number[])[j] ?? 0 }))
            : null

          return (
            <div key={i} className="bg-gray-900 rounded-lg p-4 space-y-4">
              <div className="text-sm font-medium text-blue-300">
                {cd.pass} — Per-run identifications (q ≤ 0.01)
              </div>

              {/* Precursors + Modified Peptides — shared scale */}
              <div>
                <div className="text-xs text-gray-500 mb-1 flex items-center gap-1">
                  Precursors &amp; Modified Peptides
                  <span className="text-gray-700 text-xs">(shared scale, max {globalMaxPrecPep.toLocaleString()})</span>
                </div>
                <div className="grid grid-cols-1 xl:grid-cols-2 gap-4">
                  <div>
                    <div className="text-xs text-gray-600 mb-1.5">Precursors</div>
                    <SimpleBarChart items={precItems} color="bg-blue-500" maxOverride={globalMaxPrecPep} />
                  </div>
                  <div>
                    <div className="text-xs text-gray-600 mb-1.5">Modified Peptides</div>
                    <SimpleBarChart items={pepItems} color="bg-teal-500" maxOverride={globalMaxPrecPep} />
                  </div>
                </div>
              </div>

              {/* Protein Groups — separate scale */}
              {protItems && (
                <div className="border-t border-gray-800 pt-3">
                  <div className="text-xs text-gray-500 mb-1 flex items-center gap-1">
                    Protein Groups
                    <span className="text-gray-700 text-xs">(shared scale, max {globalMaxProt.toLocaleString()})</span>
                  </div>
                  <SimpleBarChart items={protItems} color="bg-violet-500" maxOverride={globalMaxProt} />
                </div>
              )}
            </div>
          )
        }

        if (cd.type === 'enzyme_breakdown') {
          const items = cd.enzymes.map((e, j) => ({ label: e, value: cd.counts[j] ?? 0 }))
          const localMax = Math.max(...items.map(x => x.value), 1)
          return (
            <div key={i} className="bg-gray-900 rounded-lg p-4">
              <div className="text-sm font-medium text-blue-300 mb-3">
                {cd.pass} — Precursors by enzyme (q ≤ 0.01)
              </div>
              <div className="space-y-1.5">
                {items.map((d, j) => (
                  <div key={d.label} className="flex items-center gap-2 text-xs">
                    <span className="w-28 text-right text-gray-400 shrink-0">{d.label}</span>
                    <div className="flex-1 bg-gray-800 rounded overflow-hidden h-5 min-w-0">
                      <div
                        className={`h-full ${ENZYME_COLORS[j % ENZYME_COLORS.length]} rounded transition-all duration-500`}
                        style={{ width: `${(d.value / localMax) * 100}%` }}
                      />
                    </div>
                    <span className="w-16 text-gray-300 text-right shrink-0">
                      {d.value.toLocaleString()}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          )
        }

        return null
      })}
    </div>
  )
}

function PassBadge({ state }: { state?: PassState }) {
  if (!state) return null
  const icons: Record<string, string> = {
    pending:   '○',
    running:   '↻',
    done:      '✓',
    failed:    '✗',
    skipped:   '⏭',
    cancelled: '⊘',
  }
  const colors: Record<string, string> = {
    pending:   'text-gray-500',
    running:   'text-blue-400 spinning',
    done:      'text-green-400',
    failed:    'text-red-400',
    skipped:   'text-gray-500',
    cancelled: 'text-yellow-500',
  }
  return (
    <div className="flex items-center gap-1.5 bg-gray-800 rounded px-2 py-1 text-xs">
      <span className={colors[state.status] ?? 'text-gray-400'}>
        {icons[state.status] ?? '?'}
      </span>
      <span className="text-gray-300 font-mono">{state.pass}</span>
      {state.elapsedS !== undefined && (
        <span className="text-gray-500">{(state.elapsedS / 60).toFixed(1)}m</span>
      )}
    </div>
  )
}

function MetricCard({ m }: { m: Metric }) {
  return (
    <div className="bg-gray-800 rounded p-3 text-sm">
      <div className="text-blue-300 font-medium text-xs mb-2 truncate">{m.label}</div>
      <div className="grid grid-cols-3 gap-1 text-center">
        {([['Proteins', m.proteins], ['Peptides', m.peptides], ['Precursors', m.precursors]] as const).map(
          ([key, val]) => (
            <div key={key}>
              <div className="text-gray-400 text-xs">{key}</div>
              <div className="text-white font-bold text-sm">
                {typeof val === 'number' ? val.toLocaleString() : val}
              </div>
            </div>
          )
        )}
      </div>
    </div>
  )
}

// ─── Main App ─────────────────────────────────────────────────────────────────

// ─── Modifications ───────────────────────────────────────────────────────────

const MOD_SITES_RE = /^(\*n|n|[ACDEFGHIKLMNPQRSTUVWY])+$/

function fmtMass(m: number): string {
  return m.toFixed(6).replace(/0+$/, '').replace(/\.$/, '')
}

const sameMod = (a: Modification, b: Modification) =>
  a.name === b.name && a.sites === b.sites && a.type === b.type
  && !!a.monitor === !!b.monitor && !!a.noCutAfter === !!b.noCutAfter

// The DIA-NN flags the pipeline will emit (mirrors mod_flags in run_pipeline.py).
export function modFlags(job: JobConfig): string[] {
  const out: string[] = []
  for (const m of job.modifications) {
    out.push(m.type === 'fixed' ? '--fixed-mod' : '--var-mod', `${m.name},${fmtMass(m.mass)},${m.sites}`)
    // emitted per modification, immediately after it, exactly as the pipeline does
    if (m.monitor) out.push('--monitor-mod', m.name)
    if (m.noCutAfter) out.push('--no-cut-after-mod', m.name)
  }
  if (job.modifications.some(m => m.type === 'variable')) out.push('--var-mods', String(job.maxVarMods))
  if (job.metExcision) out.push('--met-excision')
  return [...out, ...job.extraDiannArgs]
}

// null when valid, else the reason (same rules as make_mod in run_pipeline.py).
function modError(name: string, mass: string, sites: string): string | null {
  if (!name.trim() || /[,\s]/.test(name.trim())) return 'name: required, no commas or spaces'
  const m = Number(mass)
  if (!mass.trim() || !isFinite(m) || m === 0 || Math.abs(m) > 2000) return 'mass: a non-zero delta mass in Da'
  if (!MOD_SITES_RE.test(sites.trim())) return "sites: amino-acid letters, n or *n"
  return null
}

function ModificationsPanel({ job, setJob }: {
  job: JobConfig
  setJob: React.Dispatch<React.SetStateAction<JobConfig>>
}) {
  const [name,  setName]  = useState('')
  const [mass,  setMass]  = useState('')
  const [sites, setSites] = useState('')
  const [type,  setType]  = useState<'fixed' | 'variable'>('variable')
  const [extra, setExtra] = useState(job.extraDiannArgs.join(' '))
  const err = (name || mass || sites) ? modError(name, mass, sites) : null

  const toggle = (mod: Modification, on: boolean) => setJob(j => ({
    ...j,
    modifications: on ? [...j.modifications.filter(m => !sameMod(m, mod)), mod]
                      : j.modifications.filter(m => !sameMod(m, mod)),
  }))
  const custom = job.modifications.filter(m => !MOD_PRESETS.some(p => sameMod(p.mod, m)))
  const flags  = modFlags(job)
  const cysFixed = job.modifications.some(m => m.type === 'fixed' && m.sites === 'C')

  return (
    <div className="mt-3">
      <label className="block text-xs text-gray-400 mb-1">Modifications (all passes)</label>
      {MOD_PRESETS.map(p => (
        <label key={p.key} className="flex items-center gap-2 text-xs text-gray-300 cursor-pointer">
          <input type="checkbox" className="accent-blue-500"
            checked={job.modifications.some(m => sameMod(m, p.mod))}
            onChange={e => toggle(p.mod, e.target.checked)} />
          {p.label}
        </label>
      ))}
      {custom.map(m => (
        <div key={`${m.name}-${m.sites}-${m.type}`} className="flex items-center gap-2 text-xs text-gray-300">
          <span className="font-mono">{m.name},{fmtMass(m.mass)},{m.sites}</span>
          <span className="text-gray-500">{m.type}</span>
          <button onClick={() => toggle(m, false)} className="ml-auto text-gray-500 hover:text-red-400">remove</button>
        </div>
      ))}
      <details className="mt-1">
        <summary className="text-xs text-gray-500 cursor-pointer hover:text-gray-400">Advanced modifications</summary>
        <div className="mt-1 grid grid-cols-4 gap-1">
          <input placeholder="UniMod:21" value={name} onChange={e => setName(e.target.value)}
            className="bg-gray-900 border border-gray-700 rounded px-1 py-1 text-xs text-gray-100" />
          <input placeholder="79.966331" value={mass} onChange={e => setMass(e.target.value)}
            className="bg-gray-900 border border-gray-700 rounded px-1 py-1 text-xs text-gray-100" />
          <input placeholder="STY" value={sites} onChange={e => setSites(e.target.value)}
            className="bg-gray-900 border border-gray-700 rounded px-1 py-1 text-xs text-gray-100" />
          <select value={type} onChange={e => setType(e.target.value as 'fixed' | 'variable')}
            className="bg-gray-900 border border-gray-700 rounded px-1 py-1 text-xs text-gray-100">
            <option value="variable">variable</option>
            <option value="fixed">fixed</option>
          </select>
        </div>
        <div className="flex items-center gap-2 mt-1">
          <button disabled={!name || !!err}
            onClick={() => {
              toggle({ name: name.trim(), mass: Number(mass), sites: sites.trim(), type }, true)
              setName(''); setMass(''); setSites('')
            }}
            className="text-xs px-2 py-0.5 rounded bg-gray-700 text-gray-200 disabled:opacity-40">Add</button>
          {err && <span className="text-xs text-red-400">{err}</span>}
        </div>
        <div className="flex items-center gap-3 mt-2 text-xs text-gray-400">
          <label className="flex items-center gap-1">Max variable
            <input type="number" min={0} max={5} value={job.maxVarMods}
              onChange={e => setJob(j => ({ ...j, maxVarMods: Math.max(0, parseInt(e.target.value) || 0) }))}
              className="w-12 bg-gray-900 border border-gray-700 rounded px-1 py-0.5 text-xs text-gray-100" />
          </label>
          <label className="flex items-center gap-1 cursor-pointer">
            <input type="checkbox" className="accent-blue-500" checked={job.metExcision}
              onChange={e => setJob(j => ({ ...j, metExcision: e.target.checked }))} />
            N-term Met excision
          </label>
        </div>
        <label className="block text-xs text-gray-500 mt-2 mb-0.5">Extra DIA-NN arguments (escape hatch, every pass)</label>
        <input value={extra} placeholder="e.g. --peak-center"
          onChange={e => setExtra(e.target.value)}
          onBlur={() => setJob(j => ({ ...j, extraDiannArgs: extra.split(/\s+/).filter(Boolean) }))}
          className="w-full bg-gray-900 border border-gray-700 rounded px-2 py-1 text-xs text-gray-100 font-mono" />
      </details>
      <div className="mt-2 text-xs">
        <span className="text-gray-500">Will be applied: </span>
        <span className="font-mono text-gray-300 break-all">{flags.length ? flags.join(' ') : 'no modifications'}</span>
      </div>
      {!cysFixed && (
        <div className="mt-1 text-xs text-amber-400">
          No fixed cysteine modification: correct only for non-alkylated samples.
        </div>
      )}
    </div>
  )
}

export default function App() {
  const [config,      setConfigState] = useState<Config>(DEFAULT_CONFIG)
  const [job,         setJob]         = useState<JobConfig>(DEFAULT_JOB)
  const [run,         setRun]         = useState<RunState | null>(null)
  const [autoScroll,  setAutoScroll]  = useState(true)
  const [tab,         setTab]         = useState<'terminal' | 'results'>('terminal')
  const [configDirty, setConfigDirty] = useState(false)
  const [resumeModal, setResumeModal] = useState<{ settingsPath: string } | null>(null)
  const [savedMsg,    setSavedMsg]    = useState('')

  const termRef = useRef<HTMLDivElement>(null)
  const wsRef   = useRef<WebSocket | null>(null)

  // Load config on mount
  useEffect(() => {
    apiFetch('/config').then(cfg => {
      setConfigState({ ...DEFAULT_CONFIG, ...cfg })
    }).catch(() => {})
  }, [])

  // Build provenance: which branch/commit this binary came from. Shown in the
  // header because mEDIA.exe is not tracked by git and can be renamed.
  const [version, setVersion] = useState<BuildInfo | null>(null)
  useEffect(() => {
    apiFetch('/version').then(setVersion).catch(() => {})
  }, [])

  // Auto-scroll terminal
  useEffect(() => {
    if (autoScroll && termRef.current) {
      termRef.current.scrollTop = termRef.current.scrollHeight
    }
  }, [run?.outputLines.length, autoScroll])

  const setConfig = (update: Partial<Config>) => {
    setConfigState(prev => ({ ...prev, ...update }))
    setConfigDirty(true)
  }

  const saveConfig = async () => {
    await apiFetch('/config', 'POST', config)
    setConfigDirty(false)
  }

  // WebSocket connection
  const connectWs = useCallback((runId: string) => {
    const proto = location.protocol === 'https:' ? 'wss' : 'ws'
    const wsUrl = `${proto}://${location.host}/ws/${runId}`
    const ws    = new WebSocket(wsUrl)
    wsRef.current = ws

    ws.onmessage = (ev) => {
      try {
        const msg = JSON.parse(ev.data)
        if (msg.type === 'line') {
          const line: string = msg.data
          setRun(prev => {
            if (!prev) return prev
            // Parse ##STATUS##, ##METRIC##, ##CHARTDATA##
            let passStatuses = prev.passStatuses
            let metrics      = prev.metrics
            let chartData    = prev.chartData
            const statusMatch    = line.match(/^##STATUS## (.+)$/)
            const metricMatch    = line.match(/^##METRIC## (.+)$/)
            const chartDataMatch = line.match(/^##CHARTDATA## (.+)$/)
            if (statusMatch) {
              try {
                const d = JSON.parse(statusMatch[1])
                passStatuses = { ...passStatuses, [d.pass]: d }
              } catch {}
            }
            if (metricMatch) {
              try { metrics = [...metrics, JSON.parse(metricMatch[1])] } catch {}
            }
            if (chartDataMatch) {
              try { chartData = [...chartData, JSON.parse(chartDataMatch[1])] } catch {}
            }
            return {
              ...prev,
              outputLines: [...prev.outputLines, line],
              passStatuses,
              metrics,
              chartData,
            }
          })
        } else if (msg.type === 'status') {
          setRun(prev => prev ? { ...prev, status: msg.data.status } : prev)
        }
      } catch {}
    }

    ws.onerror  = () => setRun(prev => prev ? { ...prev, status: 'failed' } : prev)
    ws.onclose  = () => {}
  }, [])

  const doStartRun = async () => {
    if (configDirty) await saveConfig()
    setRun({
      runId:        '',
      status:       'running',
      outputLines:  [],
      passStatuses: {},
      metrics:      [],
      chartData:    [],
    })
    setTab('terminal')
    try {
      const { runId } = await apiFetch('/runs', 'POST', { config, job })
      setRun(prev => prev ? { ...prev, runId } : null)
      connectWs(runId)
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : String(e)
      setRun(prev => prev ? {
        ...prev,
        status: 'failed',
        outputLines: [...(prev?.outputLines ?? []), `[error] ${msg}`],
      } : null)
    }
  }

  const handleRun = async () => {
    if (run?.status === 'running') return
    if (job.outDir) {
      try {
        const check = await apiFetch(`/check-outdir?path=${encodeURIComponent(job.outDir)}`)
        if (check.hasExisting) {
          setResumeModal({ settingsPath: check.settingsPath })
          return
        }
      } catch { /* ignore — just proceed */ }
    }
    await doStartRun()
  }

  const loadSettingsFile = async () => {
    try {
      const filePath = await browse('/browse/file?title=Select+pipeline.settings+file')
      if (typeof filePath !== 'string' || !filePath) return
      const s = await apiFetch(`/load-settings?path=${encodeURIComponent(filePath)}`)
      setJob(j => ({
        ...j,
        rawFiles:        Array.isArray(s.rawFiles)    ? s.rawFiles        : j.rawFiles,
        fasta:           s.fasta           ?? j.fasta,
        outDir:          s.outDir          ?? j.outDir,
        enzymes:         Array.isArray(s.enzymes)     ? s.enzymes         : j.enzymes,
        speclibPaths:    s.speclibPaths    ?? j.speclibPaths,
        semiExt:         s.semiExt         ?? j.semiExt,
        skipPass1:       s.skipPass1       ?? j.skipPass1,
        combinedLib:     s.combinedLib     ?? j.combinedLib,
        minPepLen:       s.minPepLen       ?? j.minPepLen,
        maxPepLen:       s.maxPepLen       ?? j.maxPepLen,
        missedCleavages: s.missedCleavages ?? j.missedCleavages,
        // Older settings files have no per-enzyme overrides: none, as before.
        enzymeMissedCleavages: (s.enzymeMissedCleavages && typeof s.enzymeMissedCleavages === 'object')
                                 ? s.enzymeMissedCleavages : {},
        // A settings file without the key predates modification support: that
        // run applied none, so loading it reproduces it rather than the default.
        modifications:   Array.isArray(s.modifications) ? s.modifications : [],
        maxVarMods:      s.maxVarMods      ?? 1,
        metExcision:     s.metExcision     ?? false,
        extraDiannArgs:  Array.isArray(s.extraDiannArgs) ? s.extraDiannArgs : [],
      }))
      if (s.diannExe || s.threads) {
        setConfigState(prev => ({
          ...prev,
          ...(s.diannExe  ? { diannExe: s.diannExe }   : {}),
          ...(s.threads   ? { threads:  s.threads  }   : {}),
        }))
      }
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : String(e)
      alert(`Failed to load settings: ${msg}`)
    }
  }

  const saveSettingsFile = async () => {
    if (!job.outDir) {
      alert('Set an output directory first — settings will be saved there.')
      return
    }
    const sep    = job.outDir.includes('\\') ? '\\' : '/'
    const target = job.outDir.replace(/[\\/]+$/, '') + sep + 'pipeline.settings'
    const data   = {
      timestamp:       new Date().toISOString(),
      rawFiles:        job.rawFiles,
      fasta:           job.fasta,
      outDir:          job.outDir,
      enzymes:         job.enzymes,
      speclibPaths:    job.speclibPaths,
      semiExt:         job.semiExt,
      skipPass1:       job.skipPass1,
      combinedLib:     job.combinedLib,
      minPepLen:       job.minPepLen,
      maxPepLen:       job.maxPepLen,
      missedCleavages: job.missedCleavages,
      enzymeMissedCleavages: job.enzymeMissedCleavages,
      modifications:   job.modifications,
      maxVarMods:      job.maxVarMods,
      metExcision:     job.metExcision,
      extraDiannArgs:  job.extraDiannArgs,
      diannExe:        config.diannExe,
      threads:         config.threads,
    }
    try {
      await apiFetch('/save-settings', 'POST', { path: target, data })
      setSavedMsg('✓ Saved')
      setTimeout(() => setSavedMsg(''), 2500)
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : String(e)
      alert(`Failed to save settings: ${msg}`)
    }
  }

  const handleCancel = async () => {
    if (!run?.runId) return
    try {
      // Answers once the pipeline and DIA-NN have exited (or the wait ran out).
      const r = await apiFetch(`/runs/${run.runId}`, 'DELETE')
      setRun(prev => prev ? { ...prev, status: r.status ?? 'cancelled' } : prev)
      if (r.stillRunning?.length)
        alert(`Run cancelled, but process(es) ${r.stillRunning.join(', ')} did not exit. ` +
              'Stop them (Task Manager) before reusing the output folder.')
      wsRef.current?.close()
    } catch {}
  }

  const [shutdownState, setShutdownState] = useState<'idle' | 'stopping' | 'stopped'>('idle')

  const handleShutdown = async () => {
    if (shutdownState !== 'idle') return
    const msg = isRunning
      ? 'A pipeline run is in progress — shutting down will terminate it. Continue?'
      : 'Shut down the mE-DIA server? You will need to relaunch mEDIA.exe to use it again.'
    if (!window.confirm(msg)) return
    setShutdownState('stopping')
    try { await apiFetch('/shutdown', 'POST') } catch {}
    // The server hard-exits ~0.4 s after responding, so any further requests
    // will fail. Close the websocket and flip into the "stopped" splash.
    try { wsRef.current?.close() } catch {}
    setTimeout(() => setShutdownState('stopped'), 600)
  }

  const isRunning = run?.status === 'running'

  // Classify terminal lines for colour coding
  const classifyLine = (line: string): string => {
    if (line.startsWith('##METRIC##')) return 'metric'
    if (line.startsWith('##STATUS##')) return 'status'
    if (/warning|warn/i.test(line)) return 'warn'
    if (/error|fail|critical/i.test(line)) return 'error'
    return ''
  }

  if (shutdownState === 'stopped') {
    return (
      <div className="flex h-screen bg-gray-950 text-gray-100 items-center justify-center">
        <div className="bg-gray-900 border border-gray-800 rounded-2xl shadow-2xl p-10 max-w-md text-center">
          <div className="text-4xl mb-3">🛑</div>
          <div className="text-lg font-semibold mb-1">mE-DIA server stopped</div>
          <div className="text-sm text-gray-400 mb-6">
            Port 8765 has been freed. You can close this browser tab.
          </div>
          <div className="text-xs text-gray-500">
            To resume, relaunch <span className="font-mono text-gray-300">mEDIA.exe</span>.
          </div>
        </div>
      </div>
    )
  }

  return (
    <div className="flex h-screen bg-gray-950 text-gray-100 overflow-hidden">

      {/* ── Left sidebar ──────────────────────────────────────────────── */}
      <div className="w-80 min-w-72 flex flex-col border-r border-gray-800 overflow-y-auto">

        {/* Header */}
        <div className="p-4 border-b border-gray-800 flex items-center gap-2">
          <span className="text-xl">🧬</span>
          <div className="flex-1">
            <div className="font-semibold text-sm">mE-DIA Pipeline</div>
            <div className="text-xs text-gray-500">
              Multi-enzyme DIA-NN
              {version && (
                <span
                  title={`built ${version.built_utc ?? 'unknown'}` +
                         (version.dirty ? ' from a tree with uncommitted changes' : '')}
                  className={'ml-1 font-mono ' + (version.dirty ? 'text-amber-400' : '')}>
                  {' '}· {version.branch}@{version.commit}{version.dirty ? '-dirty' : ''}
                </span>
              )}
            </div>
          </div>
          <button
            onClick={handleShutdown}
            disabled={shutdownState !== 'idle'}
            title="Stop the mE-DIA server and free port 8765"
            className="text-xs font-semibold px-2.5 py-1 rounded border border-gray-700 text-gray-400
              hover:text-red-300 hover:border-red-700 hover:bg-red-950/30
              disabled:opacity-40 disabled:cursor-not-allowed transition-colors">
            {shutdownState === 'stopping' ? 'Stopping…' : 'Quit'}
          </button>
        </div>

        <div className="p-4 space-y-1 flex-1">

          {/* ── DIA-NN Settings ────────────────────────────────────── */}
          <details open className="mb-1">
            <summary className="text-xs font-semibold text-gray-400 uppercase tracking-wider
              cursor-pointer select-none mb-2 hover:text-gray-300">
              ⚙ DIA-NN Settings
            </summary>
            <PathInput
              label="DIA-NN executable"
              value={config.diannExe}
              onChange={v => setConfig({ diannExe: v })}
              browseEp="/browse/exe"
              placeholder="diann  or  C:/DIA-NN/diann.exe"
              mono
            />
            <div className="mb-3">
              <label className="block text-xs font-medium text-gray-400 mb-1">
                CPU Threads
              </label>
              <input
                type="number"
                min={1}
                max={128}
                value={config.threads}
                onChange={e => setConfig({ threads: parseInt(e.target.value) || 14 })}
                className="w-24 bg-gray-900 border border-gray-700 rounded px-2 py-1.5
                  text-sm text-gray-100 focus:outline-none focus:border-blue-500"
              />
            </div>
            {configDirty && (
              <button
                onClick={saveConfig}
                className="text-xs px-3 py-1 bg-indigo-700 hover:bg-indigo-600 rounded
                  text-white mb-2"
              >
                Save settings
              </button>
            )}
          </details>

          {/* ── Job Setup ──────────────────────────────────────────── */}
          <details open>
            <summary className="text-xs font-semibold text-gray-400 uppercase tracking-wider
              cursor-pointer select-none mb-2 hover:text-gray-300">
              📂 Job Setup
            </summary>

            {/* Settings file buttons */}
            <div className="flex items-center gap-2 mb-3">
              <button
                onClick={loadSettingsFile}
                className="flex-1 py-1 bg-gray-800 hover:bg-gray-700 border border-gray-700
                  rounded text-xs text-gray-300 hover:text-white"
              >
                📂 Load Settings
              </button>
              <button
                onClick={saveSettingsFile}
                className="flex-1 py-1 bg-gray-800 hover:bg-gray-700 border border-gray-700
                  rounded text-xs text-gray-300 hover:text-white"
              >
                💾 Save Settings
              </button>
              {savedMsg && (
                <span className="text-xs text-green-400 whitespace-nowrap">{savedMsg}</span>
              )}
            </div>

            <RawFilesInput
              files={job.rawFiles}
              onChange={v => setJob(j => ({ ...j, rawFiles: v }))}
            />

            <PathInput
              label="FASTA database"
              value={job.fasta}
              onChange={v => setJob(j => ({ ...j, fasta: v }))}
              browseEp="/browse/fasta"
              placeholder="proteome.fasta"
              mono
            />

            <PathInput
              label="Output directory"
              value={job.outDir}
              onChange={v => setJob(j => ({ ...j, outDir: v }))}
              browseEp="/browse/dir"
              placeholder="results/"
              mono
            />

            <EnzymeSelector
              selected={job.enzymes}
              onChange={v => setJob(j => ({ ...j, enzymes: v }))}
            />

            <MissedCleavagesPanel job={job} setJob={setJob} />

            <SpeclibInputs
              enzymes={job.enzymes}
              paths={job.speclibPaths}
              onChange={v => setJob(j => ({ ...j, speclibPaths: v }))}
            />

            {/* Options */}
            <div className="space-y-2 mb-3">
              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={job.semiExt}
                  onChange={e => setJob(j => ({ ...j, semiExt: e.target.checked }))}
                  className="accent-blue-500"
                />
                <span className="text-sm text-gray-300">Semi-specific extension (pass3+4)</span>
              </label>
              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={job.skipPass1}
                  onChange={e => setJob(j => ({ ...j, skipPass1: e.target.checked }))}
                  className="accent-blue-500"
                />
                <span className="text-sm text-gray-300">Skip pass1 (provide library)</span>
              </label>
            </div>

            {job.skipPass1 && (
              <PathInput
                label="Combined library (TSV, replaces pass1)"
                value={job.combinedLib}
                onChange={v => setJob(j => ({ ...j, combinedLib: v }))}
                browseEp="/browse/file"
                placeholder="tcge_combined.tsv"
                mono
              />
            )}

            <ModificationsPanel job={job} setJob={setJob} />

            {/* Advanced */}
            <details className="mt-2">
              <summary className="text-xs text-gray-500 cursor-pointer hover:text-gray-400">
                Advanced peptide params
              </summary>
              <div className="mt-2 grid grid-cols-3 gap-2">
                {([
                  ['Min len', 'minPepLen', 4, 20] as const,
                  ['Max len', 'maxPepLen', 10, 50] as const,
                  ['Missed (default)', 'missedCleavages', 0, 3] as const,
                ]).map(([lbl, key, mn, mx]) => (
                  <div key={key}>
                    <label className="block text-xs text-gray-500 mb-0.5">{lbl}</label>
                    <input
                      type="number"
                      min={mn}
                      max={mx}
                      value={(job as unknown as Record<string, number>)[key]}
                      onChange={e =>
                        setJob(j => ({ ...j, [key]: parseInt(e.target.value) }))
                      }
                      className="w-full bg-gray-900 border border-gray-700 rounded px-2 py-1
                        text-xs text-gray-100 focus:outline-none focus:border-blue-500"
                    />
                  </div>
                ))}
              </div>
            </details>
          </details>
        </div>

        {/* Run / Cancel button */}
        <div className="p-4 border-t border-gray-800">
          {isRunning ? (
            <button
              onClick={handleCancel}
              className="w-full py-2.5 rounded font-semibold text-sm
                bg-red-700 hover:bg-red-600 text-white transition-colors"
            >
              ⏹ Cancel Run
            </button>
          ) : (
            <button
              onClick={handleRun}
              disabled={job.rawFiles.length === 0 || job.enzymes.length === 0}
              className="w-full py-2.5 rounded font-semibold text-sm
                bg-blue-600 hover:bg-blue-500 text-white transition-colors
                disabled:bg-gray-700 disabled:text-gray-500 disabled:cursor-not-allowed"
            >
              ▶ Run Pipeline
            </button>
          )}
          {run && (
            <div className={`mt-2 text-center text-xs font-medium ${
              run.status === 'done'      ? 'text-green-400' :
              run.status === 'failed'    ? 'text-red-400'   :
              run.status === 'cancelled' ? 'text-yellow-400' :
                                           'text-blue-400'
            }`}>
              {run.status === 'running'   ? '⏳ Running...' :
               run.status === 'done'      ? '✓ Done' :
               run.status === 'failed'    ? '✗ Failed' :
               run.status === 'cancelled' ? '⊘ Cancelled' : ''}
            </div>
          )}
        </div>
      </div>

      {/* ── Right panel ──────────────────────────────────────────────── */}
      <div className="flex-1 flex flex-col overflow-hidden">

        {/* Pipeline flow overview */}
        {run && (
          <PipelineFlow
            enzymes={job.enzymes}
            semiExt={job.semiExt}
            passStatuses={run.passStatuses}
          />
        )}

        {/* Tab bar */}
        <div className="flex border-b border-gray-800 bg-gray-900">
          {(['terminal', 'results'] as const).map(t => (
            <button
              key={t}
              onClick={() => setTab(t)}
              className={`px-5 py-2.5 text-sm font-medium border-b-2 transition-colors
                ${tab === t
                  ? 'border-blue-500 text-blue-400'
                  : 'border-transparent text-gray-500 hover:text-gray-300'}`}
            >
              {t === 'terminal'
                ? `Terminal${run ? ` (${run.outputLines.length} lines)` : ''}`
                : `Results${run?.metrics.length
                    ? ` (${run.metrics.length} passes${run.chartData.length
                        ? `, ${run.chartData.length} charts` : ''})` : ''}`}
            </button>
          ))}
          {tab === 'terminal' && (
            <div className="ml-auto flex items-center pr-3">
              <label className="flex items-center gap-1.5 text-xs text-gray-500 cursor-pointer">
                <input
                  type="checkbox"
                  checked={autoScroll}
                  onChange={e => setAutoScroll(e.target.checked)}
                  className="accent-blue-500"
                />
                Auto-scroll
              </label>
            </div>
          )}
        </div>

        {/* Terminal tab */}
        {tab === 'terminal' && (
          <div
            ref={termRef}
            onScroll={() => {
              if (!termRef.current) return
              const { scrollTop, scrollHeight, clientHeight } = termRef.current
              setAutoScroll(scrollTop + clientHeight >= scrollHeight - 40)
            }}
            className="flex-1 overflow-y-auto bg-gray-950 p-4 font-mono text-xs leading-5"
          >
            {!run && (
              <div className="text-gray-600 italic mt-8 text-center text-sm">
                Configure a run on the left and click ▶ Run Pipeline
              </div>
            )}
            {run?.outputLines.map((line, i) => {
              const cls = classifyLine(line)
              // Hide structured lines — already parsed into UI state
              if (line.startsWith('##STATUS##') || line.startsWith('##METRIC##')
                  || line.startsWith('##CHARTDATA##')) return null
              return (
                <div key={i} className={`terminal-line ${cls}`}>
                  {line || '\u00a0'}
                </div>
              )
            })}
            {isRunning && (
              <div className="text-blue-400 animate-pulse mt-1">█</div>
            )}
          </div>
        )}

        {/* Results tab */}
        {tab === 'results' && (
          <div className="flex-1 overflow-y-auto p-4 space-y-6">
            {!run?.metrics.length && !run?.chartData.length && (
              <div className="text-gray-600 italic mt-8 text-center text-sm">
                Results appear here as each pass completes
              </div>
            )}

            {/* ID summary cards */}
            {(run?.metrics.length ?? 0) > 0 && (
              <div>
                <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-2">
                  Identification Summary
                </h3>
                <div className="grid grid-cols-1 sm:grid-cols-2 xl:grid-cols-3 gap-3">
                  {run?.metrics.map((m, i) => <MetricCard key={i} m={m} />)}
                </div>
              </div>
            )}

            {/* Charts */}
            {(run?.chartData.length ?? 0) > 0 && (
              <div>
                <h3 className="text-xs font-semibold text-gray-500 uppercase tracking-wider mb-3">
                  Visualisations
                </h3>
                <RunCharts chartData={run?.chartData ?? []} />
              </div>
            )}
          </div>
        )}
      </div>

      {/* Resume modal */}
      {resumeModal && (
        <div className="fixed inset-0 bg-black/70 flex items-center justify-center z-50">
          <div className="bg-gray-900 border border-gray-700 rounded-lg p-6 max-w-lg w-full mx-4 shadow-2xl">
            <h3 className="text-base font-semibold text-white mb-2">
              ⚠ Existing Run Detected
            </h3>
            <p className="text-sm text-gray-300 mb-1">
              The output folder already contains pipeline results:
            </p>
            <p className="text-xs font-mono text-gray-500 bg-gray-800 rounded px-2 py-1.5 mb-3 break-all">
              {job.outDir}
            </p>
            <p className="text-sm text-gray-300 mb-5">
              The pipeline will skip passes whose outputs already exist and resume
              from the last incomplete step.
            </p>
            {resumeModal.settingsPath && (
              <p className="text-xs text-gray-500 mb-5">
                Settings from previous run:{' '}
                <span className="font-mono">{resumeModal.settingsPath.split(/[\\/]/).pop()}</span>
              </p>
            )}
            <div className="flex gap-3 justify-end">
              <button
                onClick={() => setResumeModal(null)}
                className="px-4 py-2 rounded text-sm text-gray-400 hover:text-gray-200
                  bg-gray-800 hover:bg-gray-700 border border-gray-700"
              >
                Cancel
              </button>
              <button
                onClick={() => { setResumeModal(null); doStartRun() }}
                className="px-4 py-2 rounded text-sm font-semibold text-white
                  bg-blue-600 hover:bg-blue-500"
              >
                ▶ Continue Pipeline
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
