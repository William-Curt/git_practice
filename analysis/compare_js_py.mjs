// Cross-check the JavaScript analysis (src/analysis.js) against the Python reference output.
//
//   python analysis/langmuir_analysis.py data/sweep_log.csv --quantize-float32 --out ref.json --dump-segments 1,10,30
//   node analysis/compare_js_py.mjs ref.json            # add 'none' or 'const' as a 2nd argument for those baselines
//
// The CSV is parsed here the same way src/build.py packs it (float32 columns), so both sides see identical inputs.
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { createRequire } from 'module';
const require = createRequire(import.meta.url);
const here = path.dirname(fileURLToPath(import.meta.url));
const LP = require(path.join(here, '..', 'src', 'analysis.js'));
const [refPath, mode = 'ref', tolArg = '1e-9', csvPath = path.join(here, '..', 'data', 'sweep_log.csv')] = process.argv.slice(2);
const TOL = +tolArg;

function loadCsv(file) {
  const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(l => l.length);
  const head = lines[0].split(','), col = name => head.indexOf(name);
  const iT = col('timestamp'), iD = col('DAC'), iI1 = col('I1'), iI2 = col('I2'), iV2 = col('V2');
  const n = lines.length - 1, ts = new Float32Array(n), V2 = new Float32Array(n), I1 = new Float32Array(n), I2 = new Float32Array(n), DAC = new Uint16Array(n);
  const t0 = Date.parse(lines[1].split(',')[iT] + 'Z');
  const num = s => (/^nan$/i.test(s) || s === '') ? NaN : +s;
  for (let i = 0; i < n; i++) {
    const f = lines[i + 1].split(',');
    ts[i] = (Date.parse(f[iT] + 'Z') - t0) / 1000; DAC[i] = Math.round(+f[iD]); V2[i] = num(f[iV2]); I1[i] = num(f[iI1]); I2[i] = num(f[iI2]);
  }
  return { bounds: LP.segmentBounds(DAC), cols: { ts, V2, I1, I2, DAC } };
}
const data = loadCsv(csvPath);
const ref = JSON.parse(fs.readFileSync(refPath));
const p = Object.assign({}, LP.DEFAULTS);
if (mode === 'none') p.baseline = { mode: 'none' };
if (mode === 'const') p.baseline = { mode: 'const', value: 1.48e-7 };
if (ref.params) {
  for (const [k, v] of Object.entries({ f_ion: 'fIon', window: 'window', te_lo: 'teLo', te_hi: 'teHi', area: 'area', mass: 'mass' })) if (ref.params[k] !== undefined && ref.params[k] !== null) p[v] = ref.params[k];
  if (ref.params.channel) p.channel = ref.params.channel;
  if (ref.params.baseline_value !== undefined && ref.params.baseline_value !== null && mode === 'const') p.baseline.value = ref.params.baseline_value;
}
const out = LP.analyzeAll(data, p);
const KEYS = ['Vmin', 'Vmax', 'sigma_I', 'sigma_b', 'a_i', 'b_i', 'Vf', 'dVf', 'Vp', 'Te', 'Te_r2', 'Te_lo', 'Te_hi', 'Ies', 'Iis', 'Iis_vf', 'ne', 'ni', 'lambdaD', 'Vp_expected', 'I_lo', 'I_hi', 'dIdV_max', 'hysteresis'];
const FLOOR = { a_i: 1e-11, b_i: 1e-12, Iis: 1e-11, Iis_vf: 1e-11, Ies: 1e-11, I_lo: 1e-11, I_hi: 1e-11, dIdV_max: 1e-11, sigma_I: 1e-13, sigma_b: 1e-13, Is: 1e-11, dIdV: 1e-11, Ie_s: 1e-11, Iion: 1e-11, Te_local: 1e-4 };
const rel = (a, b, key) => {
  if (a === null || a === undefined || !Number.isFinite(a)) a = null; if (b === null || b === undefined || !Number.isFinite(b)) b = null;
  if (a === null && b === null) return 0; if (a === null || b === null) return Infinity;
  return Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), FLOOR[key] || 1e-6);
};
const worst = {}, mismatches = [];
for (const s of ref.segments) {
  const r = out.results[s.k]; if (!r) { mismatches.push(`segment ${s.k} missing in JS`); continue; }
  for (const key of KEYS) { const e = rel(r[key], s[key], key); if (!(worst[key] >= e)) worst[key] = e; if (e > TOL) mismatches.push(`k=${s.k} ${key}: js=${r[key]} py=${s[key]}`); }
  const fj = r.flags.slice().sort().join(','), fp = (s.flags || []).slice().sort().join(',');
  if (fj !== fp) mismatches.push(`k=${s.k} flags: js=[${fj}] py=[${fp}]`);
  if (s.Te_window) for (let i = 0; i < 2; i++) if (rel(r.Te_window[i], s.Te_window[i]) > TOL) mismatches.push(`k=${s.k} Te_window[${i}]: js=${r.Te_window[i]} py=${s.Te_window[i]}`);
  if ((s.Te_npts || 0) !== r.Te_npts) mismatches.push(`k=${s.k} Te_npts: js=${r.Te_npts} py=${s.Te_npts}`);
  if (s.n !== r.n) mismatches.push(`k=${s.k} n: js=${r.n} py=${s.n}`);
  if (s.n_glitch !== r.nGlitch) mismatches.push(`k=${s.k} n_glitch: js=${r.nGlitch} py=${s.n_glitch}`);
  if (s.dump) {
    for (const jk of ['Is', 'dIdV', 'Ie_s', 'Iion', 'Te_local']) {
      const A = r[jk], B = s.dump[jk]; if (!A || !B) { if ((A && A.length) || (B && B.length)) mismatches.push(`k=${s.k} dump ${jk} presence differs`); continue; }
      let e = 0, nn = 0; for (let i = 0; i < B.length; i++) { const x = rel(A[i], B[i], jk); if (x === Infinity) nn++; else if (x > e) e = x; }
      if (e > TOL || nn) mismatches.push(`k=${s.k} dump ${jk}: max rel ${e.toExponential(2)} null-mismatch ${nn}`);
    }
    if (s.dump.te_mask) { let diff = 0; for (let i = 0; i < s.dump.te_mask.length; i++) if (!!s.dump.te_mask[i] !== !!r.teMask[i]) diff++; if (diff) mismatches.push(`k=${s.k} te_mask differs at ${diff} points`); }
  }
}
if (out.refLine && ref.a_r != null) for (const [a, b] of [['a_r', 'a'], ['b_r', 'b']]) if (rel(out.refLine[b], ref[a]) > TOL) mismatches.push(`baseline ${a}: js=${out.refLine[b]} py=${ref[a]}`);
console.log(`compared ${ref.segments.length} segments (baseline ${mode}, tolerance ${TOL}); worst relative differences:`);
console.log(Object.entries(worst).map(([k, v]) => `${k}=${v === Infinity ? 'NULL-MISMATCH' : v.toExponential(1)}`).join('  '));
console.log(mismatches.length ? `MISMATCHES (${mismatches.length}):\n` + mismatches.slice(0, 40).join('\n') : 'ALL MATCH');
process.exit(mismatches.length ? 1 : 0);
