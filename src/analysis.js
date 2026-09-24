/* Langmuir sweep analysis — pure functions, no DOM. Implements ANALYSIS_SPEC v2, mirroring langmuir_analysis.py step for step. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.LP = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const E = 1.602176634e-19, ME = 9.1093837015e-31, AMU = 1.66053906660e-27, EPS0 = 8.8541878128e-12;

  const DEFAULTS = {
    channel: 'I2', baseline: { mode: 'ref', refFrom: 74, refTo: 79, value: 1.48e-7 },
    fIon: 0.20, window: 31, teMode: 'auto', teLo: 0.02, teHi: 0.30, teVlo: -4, teVhi: 2,
    area: 1.0e-5, mass: 39.948,
  };

  /* ---- segmentation: half-sweeps between DAC ramp reversals ---- */
  function segmentBounds(dac) {
    // d[i] = sign(dac[i+1] - dac[i]); a new half-sweep starts at row i where d[i] != d[i-1] (the extremum row opens the next segment)
    const n = dac.length, bounds = [0];
    let prev = 0;
    for (let i = 0; i + 1 < n; i++) {
      const d = Math.sign(dac[i + 1] - dac[i]);
      if (i > 0 && d !== prev) bounds.push(i);
      prev = d;
    }
    bounds.push(n);
    return bounds;
  }

  /* ---- ordinary least squares y = a + b x over an index list (centred sums) ---- */
  function ols(x, y, idx) {
    const n = idx.length;
    let sx = 0, sy = 0;
    for (let k = 0; k < n; k++) { sx += x[idx[k]]; sy += y[idx[k]]; }
    const mx = sx / n, my = sy / n;
    let sxx = 0, sxy = 0, syy = 0;
    for (let k = 0; k < n; k++) {
      const dx = x[idx[k]] - mx, dy = y[idx[k]] - my;
      sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    if (!(sxx > 0)) return { a: NaN, b: NaN, r2: NaN, n };
    const b = sxy / sxx, a = my - b * mx;
    let ssRes = 0;
    for (let k = 0; k < n; k++) { const r = y[idx[k]] - (a + b * x[idx[k]]); ssRes += r * r; }
    const r2 = syy > 0 ? 1 - ssRes / syy : NaN;
    return { a, b, r2, n };
  }
  const range = n => { const a = new Array(n); for (let i = 0; i < n; i++) a[i] = i; return a; };

  /* ---- centred moving average, window clipped at the ends ---- */
  function smooth(I, w) {
    const n = I.length, h = (w - 1) >> 1, out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const lo = Math.max(0, i - h), hi = Math.min(n - 1, i + h);
      let s = 0; for (let j = lo; j <= hi; j++) s += I[j];
      out[i] = s / (hi - lo + 1);
    }
    return out;
  }

  /* ---- windowed OLS slope dI/dV, same clipped window ---- */
  function derivative(V, I, w) {
    const n = V.length, h = (w - 1) >> 1, out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const lo = Math.max(0, i - h), hi = Math.min(n - 1, i + h), m = hi - lo + 1;
      let sx = 0, sy = 0;
      for (let j = lo; j <= hi; j++) { sx += V[j]; sy += I[j]; }
      const mx = sx / m, my = sy / m;
      let sxx = 0, sxy = 0;
      for (let j = lo; j <= hi; j++) { const dx = V[j] - mx; sxx += dx * dx; sxy += dx * (I[j] - my); }
      out[i] = sxx > 0 ? sxy / sxx : NaN;
    }
    return out;
  }

  /* linear interpolation of (xs ascending, ys) at x, clamped at the ends */
  function interp(xs, ys, x) {
    const n = xs.length;
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= x) lo = mid; else hi = mid; }
    const d = xs[hi] - xs[lo];
    if (!(d > 0)) return ys[lo];
    const t = (x - xs[lo]) / d;
    return ys[lo] + t * (ys[hi] - ys[lo]);
  }

  /* ---- Step 0a: overrange screening in time order (always on I2) ---- */
  function badRows(data, k) {
    const b0 = data.bounds[k], b1 = data.bounds[k + 1], I2 = data.cols.I2, bad = new Uint8Array(b1 - b0);
    let n = 0;
    for (let i = b0; i < b1; i++) {
      let isBad = I2[i] < -5e-8;
      if (!isBad && i > b0 && Math.abs(I2[i] - I2[i - 1]) > 1e-6) isBad = true;
      if (!isBad && i < b1 - 1 && Math.abs(I2[i + 1] - I2[i]) > 1e-6) isBad = true;
      if (isBad) { bad[i - b0] = 1; n++; }
    }
    return { bad, n };
  }

  /* ---- Step 0b: extract a segment, screen, drop NaN, sort by V ---- */
  function prepSegment(data, k, channel) {
    const b0 = data.bounds[k], b1 = data.bounds[k + 1];
    const Vcol = data.cols.V2, Icol = data.cols[channel];
    const { bad, n: nGlitch } = badRows(data, k);
    const idx = [], excluded = [];
    for (let i = b0; i < b1; i++) {
      if (bad[i - b0]) { excluded.push(i); continue; }
      if (Number.isFinite(Icol[i])) idx.push(i);
    }
    idx.sort((p, q) => (Vcol[p] - Vcol[q]) || (p - q)); // stable by V
    const n = idx.length, V = new Float64Array(n), I = new Float64Array(n), rows = new Int32Array(n);
    for (let j = 0; j < n; j++) { V[j] = Vcol[idx[j]]; I[j] = Icol[idx[j]]; rows[j] = idx[j]; }
    const ts = data.cols.ts;
    return {
      k, V, I, rows, n, nDropped: (b1 - b0) - n, nGlitch, excluded,
      dir: (b1 - b0) > 1 && Vcol[b1 - 1] > Vcol[b0] ? 'up' : 'down',
      tStart: ts[b0], duration: ts[b1 - 1] - ts[b0], nRows: b1 - b0,
    };
  }

  /* ---- Step 0c: baseline ---- */
  function referenceLine(data, refFrom, refTo, channel) {
    const Vcol = data.cols.V2, Icol = data.cols[channel];
    const xs = [], ys = [];
    const lo = Math.max(0, Math.min(refFrom, refTo)), hi = Math.min(data.bounds.length - 2, Math.max(refFrom, refTo));
    for (let k = lo; k <= hi; k++) {
      const { bad } = badRows(data, k), b0 = data.bounds[k], b1 = data.bounds[k + 1];
      for (let i = b0; i < b1; i++) if (!bad[i - b0] && Number.isFinite(Icol[i])) { xs.push(Vcol[i]); ys.push(Icol[i]); }
    }
    if (xs.length < 50) return null;
    const f = ols(xs, ys, range(xs.length));
    return { a: f.a, b: f.b, n: xs.length, from: lo, to: hi };
  }
  function applyBaseline(seg, baseline, refLine) {
    const I = new Float64Array(seg.I);
    if (baseline.mode === 'const') { for (let i = 0; i < I.length; i++) I[i] -= baseline.value; }
    else if (baseline.mode === 'ref' && refLine) { for (let i = 0; i < I.length; i++) I[i] -= refLine.a + refLine.b * seg.V[i]; }
    return I;
  }

  /* ---- Steps 1–10 on a prepared (V, I) ---- */
  function analyze(V, I, p) {
    const n = V.length, flags = [];
    const res = { n, flags, Vmin: NaN, Vmax: NaN, sigma_I: NaN, sigma_b: NaN, a_i: NaN, b_i: NaN, Vf: null, dVf: null, Vp: null, Te: null, Te_r2: NaN,
      Te_window: [NaN, NaN], Te_npts: 0, Te_lo: null, Te_hi: null, te_alpha: NaN, te_beta: NaN, Ies: null, Iis: null, Iis_vf: null, ne: null, ni: null, lambdaD: null, Vp_expected: null,
      I_lo: NaN, I_hi: NaN, dIdV_max: NaN, hysteresis: null, Is: null, dIdV: null, Ie: null, Ie_s: null, Iion: null, Te_local: null, teMask: null, ionMask: null, ip: -1, lo: 0, hi: 0 };
    if (p.nGlitch > 0) flags.push('overrange_event');
    if (n < 50) { flags.push('insufficient'); return res; }
    const w = p.window | 1, h = (w - 1) >> 1;
    const Vmin = V[0], Vmax = V[n - 1];
    res.Vmin = Vmin; res.Vmax = Vmax;

    // Steps 1–2
    const Is = smooth(I, w), dIdV = derivative(V, I, w);
    const lo = h, hi = n - 1 - h;
    let ss = 0; for (let i = 0; i < n; i++) { const r = I[i] - Is[i]; ss += r * r; }
    const sigma_I = Math.sqrt(ss / n);
    const dV = (Vmax - Vmin) / (n - 1);
    const sigma_b = sigma_I / Math.sqrt(dV * dV * (w * w * w - w) / 12);
    Object.assign(res, { Is, dIdV, lo, hi, sigma_I, sigma_b, I_lo: Is[lo], I_hi: Is[hi] });

    // Step 2b coverage
    if (Vmin > -8 || Vmax < 8) { flags.push('partial_sweep'); return res; }

    // Step 3 ion-saturation fit
    const vIon = Vmin + p.fIon * (Vmax - Vmin);
    const ionIdx = [], ionMask = new Uint8Array(n);
    for (let i = 0; i < n; i++) if (V[i] <= vIon) { ionIdx.push(i); ionMask[i] = 1; }
    res.ionMask = ionMask;
    if (ionIdx.length < 5) { flags.push('ion_fit_insufficient'); return res; }
    const fi = ols(V, I, ionIdx); const a_i = fi.a, b_i = fi.b;
    res.a_i = a_i; res.b_i = b_i;

    // Step 4
    const Iion = new Float64Array(n), Ie = new Float64Array(n), Ie_s = new Float64Array(n);
    for (let i = 0; i < n; i++) { Iion[i] = a_i + b_i * V[i]; Ie[i] = I[i] - Iion[i]; Ie_s[i] = Is[i] - Iion[i]; }
    Object.assign(res, { Iion, Ie, Ie_s });

    // Step 5 floating potential
    const cross = [];
    for (let i = 0; i + 1 < n; i++) if (Is[i] < 0 && Is[i + 1] >= 0) cross.push(V[i] + (0 - Is[i]) * (V[i + 1] - V[i]) / (Is[i + 1] - Is[i]));
    let Vf = null;
    if (cross.length === 0) flags.push('no_zero_crossing');
    else if (cross.length === 1) Vf = cross[0];
    else {
      flags.push('vf_multiple_crossings');
      const Vc = cross[Math.floor((cross.length - 1) / 2)];
      const near = []; for (let i = 0; i < n; i++) if (Math.abs(V[i] - Vc) <= 0.5) near.push(i);
      const f = near.length >= 2 ? ols(V, Is, near) : { b: NaN };
      Vf = f.b > 0 ? -f.a / f.b : Vc;
    }
    let dVf = null;
    if (Vf !== null) {
      let iVf = 0, bd = Infinity; for (let i = 0; i < n; i++) { const d = Math.abs(V[i] - Vf); if (d < bd) { bd = d; iVf = i; } }
      dVf = dIdV[iVf] > 0 ? sigma_I / dIdV[iVf] : null;
    }
    res.Vf = Vf; res.dVf = dVf;

    // Step 6 signal gate + plasma potential
    let ip = lo; for (let i = lo; i <= hi; i++) if (dIdV[i] > dIdV[ip]) ip = i;
    res.ip = ip; res.dIdV_max = dIdV[ip];
    const span = Is[hi] - Is[lo];
    if (span < 10 * sigma_I || dIdV[ip] < 5 * sigma_b) {
      for (const f of ['no_zero_crossing', 'vf_multiple_crossings']) { const i = flags.indexOf(f); if (i >= 0) flags.splice(i, 1); }
      flags.push('no_plasma_signal'); res.Vf = null; res.dVf = null; return res;
    }
    let Vp = null;
    const topUnreliable = p.nGlitch > 0;
    if (topUnreliable) flags.push('top_of_sweep_unreliable');
    else {
      if (V[ip] <= Vmax - 1.0 && Is[hi] >= Is[ip] - 3 * sigma_I) { for (let j = ip + 1; j <= hi; j++) if (dIdV[j] <= 0.7 * dIdV[ip]) { Vp = V[ip]; break; } }
      if (Vp === null) { flags.push('vp_beyond_range'); flags.push('not_saturated'); }
    }
    res.Vp = Vp;

    // Step 7 electron temperature
    let iRef = ip;
    if (Vp === null) { iRef = lo; for (let i = lo; i <= hi; i++) if (Ie_s[i] > Ie_s[iRef]) iRef = i; }
    const IeRef = Ie_s[iRef];
    const teIdx = [], teMask = new Uint8Array(n), noise3 = 3 * sigma_I;
    if (p.teMode === 'manual') {
      for (let j = 0; j < n; j++) if (V[j] >= p.teVlo && V[j] <= p.teVhi && Ie_s[j] > noise3) { teIdx.push(j); teMask[j] = 1; }
    } else {
      let j = iRef;
      while (j >= 0 && Ie_s[j] > p.teHi * IeRef) j--;
      while (j >= 0 && Ie_s[j] >= p.teLo * IeRef && Ie_s[j] > noise3 && (Vf === null || V[j] > Vf)) { teIdx.push(j); teMask[j] = 1; j--; }
      teIdx.reverse();
    }
    res.teMask = teMask; res.Te_npts = teIdx.length;
    if (teIdx.length) { let vlo = Infinity, vhi = -Infinity; for (const j of teIdx) { if (V[j] < vlo) vlo = V[j]; if (V[j] > vhi) vhi = V[j]; } res.Te_window = [vlo, vhi]; }
    let Te = null;
    if (teIdx.length < 8) flags.push('te_insufficient');
    else {
      const lnIe = new Float64Array(n); for (const j of teIdx) lnIe[j] = Math.log(Ie_s[j]);
      const f = ols(V, lnIe, teIdx);
      res.te_alpha = f.a; res.te_beta = f.b; res.Te_r2 = f.r2;
      if (f.b > 0) Te = 1 / f.b; else flags.push('te_negative_slope');
      if (!(f.r2 >= 0.9)) flags.push('te_poor_fit');
      const sorted = teIdx.slice().sort((a, b) => V[a] - V[b]), m = sorted.length, half = Math.floor(m / 2);
      const A = sorted.slice(0, half), B = sorted.slice(half);
      if (A.length >= 4 && B.length >= 4) {
        const fa = ols(V, lnIe, A), fb = ols(V, lnIe, B);
        res.Te_lo = fa.b > 0 ? 1 / fa.b : null; res.Te_hi = fb.b > 0 ? 1 / fb.b : null;
        if (res.Te_lo !== null && res.Te_hi !== null && res.Te_hi / res.Te_lo > 1.3) flags.push('te_nonexponential');
      }
    }
    res.Te = Te;

    // Step 8
    res.Ies = topUnreliable ? null : Ie_s[iRef];

    // Step 9 ion saturation current
    let s9 = 0, c9 = 0; for (let i = 0; i < n; i++) if (V[i] <= Vmin + 1.0) { s9 += Is[i]; c9++; }
    const IisVmin = c9 ? s9 / c9 : NaN;
    res.Iis = IisVmin; res.Iis_vf = Vf !== null ? a_i + b_i * Vf : null;
    const iisNoise = Math.abs(IisVmin) < noise3;
    if (iisNoise) flags.push('iis_below_noise');
    if (IisVmin >= 0) flags.push('ion_current_positive');

    // Step 10 densities
    if (Te !== null && res.Ies !== null && res.Ies > 0) {
      res.ne = res.Ies / (E * p.area * Math.sqrt(E * Te / (2 * Math.PI * ME)));
      res.lambdaD = Math.sqrt(EPS0 * Te / (res.ne * E));
    }
    if (Te !== null && IisVmin < 0 && !iisNoise) res.ni = Math.abs(IisVmin) / (0.61 * E * p.area * Math.sqrt(E * Te / (p.mass * AMU)));
    if (Te !== null && Vf !== null) res.Vp_expected = Vf + Te * Math.log(Math.sqrt(p.mass * AMU / (2 * Math.PI * ME)) / 0.61);

    // Step 11 local slope temperature
    const Te_local = new Float64Array(n).fill(NaN);
    for (let i = 0; i < n; i++) { const d = dIdV[i] - b_i; if (Ie_s[i] > noise3 && d > 0) Te_local[i] = Ie_s[i] / d; }
    res.Te_local = Te_local;
    return res;
  }

  /* ---- Step 11 hysteresis: h minimising rms of Is_up(V) − Is_down(V + h); positive h = down-sweep lags ---- */
  function hysteresis(up, down) {
    const VU = up.V, IU = up.Is, VD = down.V, ID = down.Is;
    const idx = []; for (let i = 0; i < VU.length; i++) if (VU[i] >= -5 && VU[i] <= 9) idx.push(i);
    if (idx.length < 10) return null;
    const rms = h => { let acc = 0; for (const i of idx) { const d = IU[i] - interp(VD, ID, VU[i] + h); acc += d * d; } return Math.sqrt(acc / idx.length); };
    let h0 = 0, best = Infinity;
    for (let i = 0; i <= 40; i++) { const h = Math.round((-2 + 0.1 * i) * 100) / 100, r = rms(h); if (r < best) { best = r; h0 = h; } }
    let h1 = h0; best = Infinity;
    for (let i = 0; i <= 10; i++) { const h = Math.round((h0 - 0.1 + 0.02 * i) * 100) / 100, r = rms(h); if (r < best) { best = r; h1 = h; } }
    return Math.abs(h1) >= 1.98 ? null : h1;
  }
  function partnerOf(k, nSeg) { let p = (k % 2 === 0) ? k + 1 : k - 1; if (p > nSeg - 1) p = nSeg - 2; return p; }
  const fullSignal = r => r && r.Is && !r.flags.includes('partial_sweep') && !r.flags.includes('insufficient') && !r.flags.includes('no_plasma_signal') && !r.flags.includes('ion_fit_insufficient');
  const pairable = r => fullSignal(r) && !r.flags.includes('overrange_event');

  /* ---- full pipeline over all segments ---- */
  function analyzeSegments(data, p) {
    const nSeg = data.bounds.length - 1, flagsAll = [];
    let refLine = null;
    if (p.baseline.mode === 'ref') { refLine = referenceLine(data, p.baseline.refFrom, p.baseline.refTo, p.channel); if (!refLine) flagsAll.push('baseline_unavailable'); }
    const results = new Array(nSeg);
    for (let k = 0; k < nSeg; k++) {
      const seg = prepSegment(data, k, p.channel);
      const I = applyBaseline(seg, refLine ? p.baseline : (p.baseline.mode === 'ref' ? { mode: 'none' } : p.baseline), refLine);
      const r = analyze(seg.V, I, Object.assign({}, p, { nGlitch: seg.nGlitch }));
      for (const f of flagsAll) r.flags.push(f);
      results[k] = Object.assign(r, { k, dir: seg.dir, tStart: seg.tStart, duration: seg.duration, nDropped: seg.nDropped, nGlitch: seg.nGlitch, excluded: seg.excluded, nRows: seg.nRows, V: seg.V, I, rows: seg.rows });
    }
    return { results, refLine };
  }
  /* pair step; returns an array of per-segment hysteresis values (null where undefined) and flags pairs it could not resolve */
  function hysteresisAll(results) {
    const nSeg = results.length, out = new Array(nSeg).fill(null), unresolved = new Uint8Array(nSeg);
    for (let k = 0; k < nSeg; k++) {
      const q = partnerOf(k, nSeg), a = results[k], b = results[q];
      if (q < k && out[k] !== null) continue;
      if (q < k && (out[q] !== null || unresolved[q])) { out[k] = out[q]; unresolved[k] = unresolved[q]; continue; }
      if (!pairable(a) || !pairable(b) || a.dir === b.dir) continue;
      const hv = a.dir === 'up' ? hysteresis(a, b) : hysteresis(b, a);
      out[k] = hv; out[q] = hv;
      if (hv === null) { unresolved[k] = 1; unresolved[q] = 1; }
    }
    return { values: out, unresolved };
  }
  function applyHysteresis(results, hys) {
    for (let k = 0; k < results.length; k++) {
      const r = results[k]; r.hysteresis = hys.values[k];
      const i = r.flags.indexOf('hysteresis_unresolved'); if (i >= 0) r.flags.splice(i, 1);
      if (hys.unresolved[k]) r.flags.push('hysteresis_unresolved');
    }
  }
  function analyzeAll(data, p) {
    const out = analyzeSegments(data, p);
    applyHysteresis(out.results, hysteresisAll(out.results));
    return out;
  }

  return { DEFAULTS, CONST: { E, ME, AMU, EPS0 }, segmentBounds, ols, smooth, derivative, interp, badRows, prepSegment, referenceLine, applyBaseline, analyze, hysteresis, partnerOf, analyzeSegments, hysteresisAll, applyHysteresis, analyzeAll, fullSignal, pairable };
});
