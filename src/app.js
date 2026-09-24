/* Langmuir Sweep Bench — UI, charts, playback. Depends on d3 (global) and LP (analysis.js). */
(function () {
  'use strict';
  const $ = (s, el) => (el || document).querySelector(s);
  const $$ = (s, el) => Array.from((el || document).querySelectorAll(s));
  const reduceMotion = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
  const NG = 1457, DAC_STEP = 45;

  /* ======================= data ======================= */
  function b64ToBytes(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < u8.length; i++) u8[i] = bin.charCodeAt(i);
    return u8;
  }
  function decode(raw) {
    const enc = raw.meta.encodings || {}, cols = {};
    for (const name of raw.meta.columns) {
      const bytes = b64ToBytes(raw.cols[name]);
      const e = enc[name] || {};
      const dtype = (e.dtype || e.type || (name === 'DAC' || name === 'GSE_I' ? 'uint16' : 'float32')).toLowerCase();
      let arr = dtype.indexOf('16') >= 0 ? new Uint16Array(bytes.buffer, 0, bytes.length >> 1) : new Float32Array(bytes.buffer, 0, bytes.length >> 2);
      if (e.scale !== undefined || e.offset !== undefined) {
        const f = new Float32Array(arr.length), sc = e.scale === undefined ? 1 : e.scale, off = e.offset || 0;
        for (let i = 0; i < f.length; i++) f[i] = arr[i] * sc + off;
        arr = f;
      }
      cols[name] = arr;
    }
    return { meta: raw.meta, bounds: raw.bounds, cols, nSeg: raw.bounds.length - 1 };
  }
  const DATA = decode(window.SWEEP_DATA);
  const NSEG = DATA.nSeg, LAST = NSEG - 1;
  const GAS_LABEL = { '39.948': 'Ar', '4.0026': 'He', '20.18': 'Ne', '83.798': 'Kr', '131.29': 'Xe', '28.014': 'N₂', '31.998': 'O₂', '2.016': 'H₂' };

  /* ======================= state ======================= */
  const state = {
    sweep: 1, playing: false, speed: 1, loop: true,
    channel: 'I2', overlay: 'single', yscale: 'lin', quantity: 'I', lower: 'deriv',
    layers: { points: true, smooth: true, fits: true, deriv: true },
    baseline: { mode: 'ref', refFrom: Math.max(0, NSEG - 7), refTo: Math.max(0, NSEG - 2), value: 1.48e-7 },
    params: { fIon: 0.2, window: 31, teMode: 'auto', teLo: 0.02, teHi: 0.30, teVlo: -4, teVhi: 2, area: 1.0e-5, mass: 39.948 },
    tableOpen: false,
  };
  try {
    const saved = JSON.parse(localStorage.getItem('lsb-state-v2') || 'null');
    if (saved && typeof saved === 'object') {
      for (const k of ['sweep', 'speed', 'loop', 'channel', 'overlay', 'yscale', 'quantity', 'lower', 'tableOpen']) if (k in saved) state[k] = saved[k];
      if (saved.layers) Object.assign(state.layers, saved.layers);
      if (saved.baseline) Object.assign(state.baseline, saved.baseline);
      if (saved.params) Object.assign(state.params, saved.params);
      state.sweep = Math.min(LAST, Math.max(0, state.sweep | 0));
      state.baseline.refFrom = Math.min(LAST, Math.max(0, state.baseline.refFrom | 0)); state.baseline.refTo = Math.min(LAST, Math.max(0, state.baseline.refTo | 0));
    }
  } catch (e) { /* storage unavailable: defaults */ }
  function persist() {
    try { const { playing, ...rest } = state; localStorage.setItem('lsb-state-v2', JSON.stringify(rest)); } catch (e) { /* ignore */ }
  }

  /* ======================= analysis ======================= */
  let results = [], refLine = null, gridCache = [], hysCache = { key: '', value: null };
  function analysisParams() { return Object.assign({ channel: state.channel, baseline: state.baseline }, state.params); }
  function recomputeAll() {
    const p = analysisParams();
    const out = LP.analyzeSegments(DATA, p);
    results = out.results; refLine = out.refLine; gridCache = new Array(NSEG);
    // the pair step only depends on the smoothed curves: cache it across Te/ion/probe parameter changes
    const key = [state.channel, JSON.stringify(state.baseline), p.window, p.fIon].join('|');
    if (hysCache.key !== key) hysCache = { key, value: LP.hysteresisAll(results) };
    LP.applyHysteresis(results, hysCache.value);
  }
  /* per-sweep arrays on the common DAC grid, ascending V (index g = 1456 - DAC/45) */
  function grid(k) {
    if (gridCache[k]) return gridCache[k];
    const r = results[k], V = new Float64Array(NG).fill(NaN), R = new Float64Array(NG).fill(NaN), S = new Float64Array(NG).fill(NaN),
      D = new Float64Array(NG).fill(NaN), E = new Float64Array(NG).fill(NaN), ES = new Float64Array(NG).fill(NaN), TL = new Float64Array(NG).fill(NaN);
    const dac = DATA.cols.DAC;
    for (let i = 0; i < r.n; i++) {
      const g = NG - 1 - Math.round(dac[r.rows[i]] / DAC_STEP);
      if (g < 0 || g >= NG) continue;
      V[g] = r.V[i]; R[g] = r.I[i];
      if (r.Is) { S[g] = r.Is[i]; D[g] = r.dIdV[i]; }
      if (r.Ie) { E[g] = r.Ie[i]; ES[g] = r.Ie_s[i]; }
      if (r.Te_local) TL[g] = r.Te_local[i];
    }
    // excluded (overrange) samples for plotting, baseline-corrected like the rest
    const ex = r.excluded.map(i => [DATA.cols.V2[i], DATA.cols[state.channel][i] - (state.baseline.mode === 'const' ? state.baseline.value : (state.baseline.mode === 'ref' && refLine ? refLine.a + refLine.b * DATA.cols.V2[i] : 0))]).filter(d => Number.isFinite(d[1]));
    return (gridCache[k] = { V, R, S, D, E, ES, TL, ex });
  }
  const partnerOf = k => LP.partnerOf(k, NSEG);

  /* ======================= formatting ======================= */
  const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  const sup = s => String(s).split('').map(c => SUP[c] || c).join('');
  function unitFor(maxAbs) {
    if (!(maxAbs > 0)) return { s: 1e9, u: 'nA' };
    if (maxAbs >= 1e-6) return { s: 1e6, u: 'µA' };
    if (maxAbs >= 1e-9) return { s: 1e9, u: 'nA' };
    return { s: 1e12, u: 'pA' };
  }
  function fmtA(a, digits) {
    if (a === null || a === undefined || !Number.isFinite(a)) return null;
    const u = unitFor(Math.abs(a)); const v = a * u.s;
    const d = digits === undefined ? (Math.abs(v) >= 100 ? 0 : Math.abs(v) >= 10 ? 1 : 2) : digits;
    return { v: (v < 0 ? '−' : '') + Math.abs(v).toFixed(d), u: u.u };
  }
  const fmtAstr = (a, d) => { const f = fmtA(a, d); return f ? f.v + ' ' + f.u : '—'; };
  const fmtV = (v, d) => (v === null || v === undefined || !Number.isFinite(v)) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(d === undefined ? 2 : d) + ' V';
  function fmtSci(n) {
    if (n === null || n === undefined || !Number.isFinite(n) || n <= 0) return '—';
    const ex = Math.floor(Math.log10(n)), m = n / Math.pow(10, ex);
    return m.toFixed(2) + '×10' + sup(ex);
  }
  function fmtLen(m) {
    if (m === null || !Number.isFinite(m)) return '—';
    if (m >= 1e-3) return (m * 1e3).toFixed(2) + ' mm';
    if (m >= 1e-6) return (m * 1e6).toFixed(1) + ' µm';
    return (m * 1e9).toFixed(0) + ' nm';
  }
  const fmtT = s => (s < 60 ? s.toFixed(1) + ' s' : Math.floor(s / 60) + ' min ' + (s % 60).toFixed(0).padStart(2, '0') + ' s');
  const T0 = new Date(DATA.meta.t0_iso + (/Z|[+-]\d\d:\d\d$/.test(DATA.meta.t0_iso) ? '' : 'Z')).getTime();
  const clockAt = s => new Date(T0 + s * 1000).toISOString().slice(11, 19);
  const FLAG_TEXT = {
    partial_sweep: ['Partial sweep (does not reach ±8 V): shown but not analysed', 'warn'],
    insufficient: ['Too few samples in this sweep', 'crit'],
    overrange_event: ['Over-range transient: wrapped samples excluded (grey)', 'crit'],
    top_of_sweep_unreliable: ['Top of sweep lost to the over-range event: no V_p, I_es or n_e', 'warn'],
    hysteresis_unresolved: ['Up/down pair could not be aligned within ±2 V', 'warn'],
    no_plasma_signal: ['No plasma signal above the instrument noise', 'warn'],
    no_zero_crossing: ['Current never crosses zero: V_f outside the sweep (or baseline needed)', 'warn'],
    vf_multiple_crossings: ['Several zero crossings: V_f from a local line fit', 'warn'],
    vp_beyond_range: ['No knee inside the sweep: V_p is beyond +10 V', 'warn'],
    not_saturated: ['Electron current not saturated: I_es and n_e are lower bounds (at +10 V)', 'warn'],
    te_insufficient: ['Too few points in the T_e window', 'crit'],
    te_negative_slope: ['ln I_e slope ≤ 0: no exponential transition found', 'crit'],
    te_poor_fit: ['T_e fit R² below 0.9', 'warn'],
    te_nonexponential: ['Electron branch is not a single exponential: slope T_e differs by >30 % between the window halves', 'info'],
    iis_below_noise: ['Ion current below 3 σ noise: n_i not computed', 'warn'],
    ion_current_positive: ['Ion-side current is positive: offset not removed or no ion current', 'warn'],
    ion_fit_insufficient: ['Too few points for the ion-saturation fit', 'crit'],
    baseline_unavailable: ['Plasma-off reference has too few samples: no baseline applied', 'crit'],
  };

  /* ======================= theme colours for canvas ======================= */
  const tok = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  function hexA(hex, a) {
    const h = hex.replace('#', ''); const n = parseInt(h.length === 3 ? h.split('').map(c => c + c).join('') : h, 16);
    return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
  }

  /* ======================= DOM refs ======================= */
  const el = {
    play: $('#play'), sweep: $('#sweep'), sweepLabel: $('#sweep-label'), sweepTime: $('#sweep-time'), loop: $('#loop'),
    strip: $('#strip'), stripRight: $('#strip-right'), stripLeft: $('#strip-left'),
    ivChart: $('#iv-chart'), ivSvg: $('#iv-svg'), ivCanvas: $('#iv-canvas'), ivTip: $('#iv-tip'), ivTitle: $('#iv-title'), ivCaption: $('#iv-caption'),
    lgPair: $('#lg-pair'), lgGhost: $('#lg-ghost'), lgExcl: $('#lg-excl'),
    tiles: $('#tiles'), flags: $('#flags'), anSub: $('#an-sub'),
    baseRefWrap: $('#base-ref-wrap'), baseFrom: $('#base-from'), baseTo: $('#base-to'), baseFit: $('#base-fit'), baseConstWrap: $('#base-const-wrap'), baseConst: $('#base-const'),
    pWindow: $('#p-window'), vWindow: $('#v-window'), pFion: $('#p-fion'), vFion: $('#v-fion'),
    pTelo: $('#p-telo'), vTelo: $('#v-telo'), pTehi: $('#p-tehi'), vTehi: $('#v-tehi'),
    pVlo: $('#p-vlo'), vVlo: $('#v-vlo'), pVhi: $('#p-vhi'), vVhi: $('#v-vhi'),
    tewinAuto: $('#tewin-auto'), tewinManual: $('#tewin-manual'),
    pArea: $('#p-area'), pGas: $('#p-gas'), pMass: $('#p-mass'), pMassWrap: $('#p-mass-wrap'),
    multiples: $('#multiples'), tblToggle: $('#tbl-toggle'), tblCopy: $('#tbl-copy'), tblWrap: $('#tbl-wrap'), tbl: $('#tbl'),
    facts: $('#facts'),
  };

  /* ======================= header facts ======================= */
  (function facts() {
    const n = DATA.cols.ts.length, ts = DATA.cols.ts;
    $('#src-name').textContent = DATA.meta.source || 'sweep_log.csv';
    $('#src-when').textContent = new Date(T0).toISOString().slice(0, 16).replace('T', ' ') + ' → ' + clockAt(ts[n - 1]).slice(0, 5) + ' UTC (' + fmtT(ts[n - 1]) + ')';
    $('#src-rows').textContent = n.toLocaleString();
    $('#src-sweeps').textContent = String(NSEG);
    const ext = a => { let lo = Infinity, hi = -Infinity; for (let i = 0; i < a.length; i++) { const v = a[i]; if (!Number.isFinite(v)) continue; if (v < lo) lo = v; if (v > hi) hi = v; } return [lo, hi]; };
    const [vmin, vmax] = ext(DATA.cols.V2), [imin, imax] = ext(DATA.cols.I2), [gmin, gmax] = ext(DATA.cols.GSE_I);
    const chips = [
      ['Bias', `${fmtV(vmin, 1)} … ${fmtV(vmax, 1)}`],
      ['Ramp', `${NG} steps · ${DAC_STEP} DAC counts`],
      ['Half-sweep', `${(ts[DATA.bounds[2] - 1] - ts[DATA.bounds[1]]).toFixed(1)} s`],
      ['I2', `${fmtAstr(imin)} … ${fmtAstr(imax)}`],
      ['GSE supply', `${(gmin * 1e3).toFixed(1)} … ${(gmax * 1e3).toFixed(1)} mA`],
    ];
    el.facts.innerHTML = '';
    for (const [k, v] of chips) { const c = document.createElement('span'); c.className = 'chip mono'; const b = document.createElement('b'); b.textContent = k; c.appendChild(b); c.appendChild(document.createTextNode(' ' + v)); el.facts.appendChild(c); }
  })();

  /* ======================= segmented controls ======================= */
  function segControl(id, get, set) {
    const root = $('#' + id); if (!root) return () => {};
    const sync = () => $$('button', root).forEach(b => b.setAttribute('aria-pressed', String(b.dataset.v === String(get()))));
    root.addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; set(b.dataset.v); sync(); });
    sync();
    return sync;
  }

  /* ======================= I–V chart ======================= */
  const iv = { w: 0, mainH: 340, derivH: 120, gap: 22, m: { l: 64, r: 18, t: 14, b: 30 }, tween: null, cur: null };
  const svg = d3.select(el.ivSvg);
  const gRoot = svg.append('g');
  const gGrid = gRoot.append('g'), gShade = gRoot.append('g'), gAxes = gRoot.append('g'), gLines = gRoot.append('g'), gLower = gRoot.append('g'), gFits = gRoot.append('g'), gMark = gRoot.append('g'), gRead = gRoot.append('g'), gHover = gRoot.append('g'), gHit = gRoot.append('g');
  const defs = svg.append('defs');
  defs.append('clipPath').attr('id', 'clip-main').append('rect');
  defs.append('clipPath').attr('id', 'clip-lower').append('rect');
  gLines.attr('clip-path', 'url(#clip-main)'); gFits.attr('clip-path', 'url(#clip-main)'); gLower.attr('clip-path', 'url(#clip-lower)');
  const pathMain = gLines.append('path').attr('class', 'line-main');
  const pathPair = gLines.append('path').attr('class', 'line-pair');
  const pathLower = gLower.append('path').attr('class', 'line-deriv');
  const pathIon = gFits.append('path').attr('class', 'line-ion');
  const pathTe = gFits.append('path').attr('class', 'line-te');
  const shadeIon = gShade.append('rect').attr('class', 'ion-shade');
  const shadeTe = gShade.append('rect').attr('class', 'win-shade');
  const shadeIonD = gShade.append('rect').attr('class', 'ion-shade');
  const shadeTeD = gShade.append('rect').attr('class', 'win-shade');
  const crossX = gHover.append('line').attr('class', 'crosshair');
  const hoverDot = gHover.append('circle').attr('class', 'hover-dot').attr('r', 4.5);
  const hoverDotD = gHover.append('circle').attr('class', 'hover-dot').attr('r', 4);
  const hit = gHit.append('rect').attr('class', 'hit');
  const ctx = el.ivCanvas.getContext('2d');
  let ghostCanvas = null, ghostKey = '';

  function layout() {
    const w = Math.max(320, el.ivChart.clientWidth || 640);
    iv.w = w;
    iv.mainH = w < 520 ? 260 : 340;
    iv.derivH = w < 520 ? 100 : 120;
    iv.m.l = w < 520 ? 52 : 64;
    const showD = state.layers.deriv;
    iv.h = iv.m.t + iv.mainH + (showD ? iv.gap + iv.derivH : 0) + iv.m.b;
    svg.attr('viewBox', `0 0 ${w} ${iv.h}`).attr('width', w).attr('height', iv.h);
    const dpr = window.devicePixelRatio || 1;
    el.ivCanvas.width = Math.round(w * dpr); el.ivCanvas.height = Math.round(iv.h * dpr);
    el.ivCanvas.style.width = w + 'px'; el.ivCanvas.style.height = iv.h + 'px';
    iv.x = d3.scaleLinear().range([iv.m.l, w - iv.m.r]);
    iv.yMainRange = [iv.m.t + iv.mainH, iv.m.t];
    iv.yLowerRange = [iv.m.t + iv.mainH + iv.gap + iv.derivH, iv.m.t + iv.mainH + iv.gap];
    d3.select('#clip-main rect').attr('x', iv.m.l).attr('y', iv.m.t - 2).attr('width', w - iv.m.l - iv.m.r).attr('height', iv.mainH + 4);
    d3.select('#clip-lower rect').attr('x', iv.m.l).attr('y', iv.yLowerRange[1] - 2).attr('width', w - iv.m.l - iv.m.r).attr('height', iv.derivH + 4);
    hit.attr('x', iv.m.l).attr('y', iv.m.t).attr('width', w - iv.m.l - iv.m.r).attr('height', iv.h - iv.m.t - iv.m.b);
    ghostKey = '';
  }

  /* what is plotted on y for a grid: measured I or electron current */
  function yValues(g, kind) {
    const q = state.quantity;
    if (kind === 'raw') return q === 'Ie' ? g.E : g.R;
    return q === 'Ie' ? g.ES : g.S;
  }
  const lowerValues = g => state.lower === 'telocal' ? g.TL : g.D;
  function yExtent(arr, log) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < arr.length; i++) { const v = log ? Math.abs(arr[i]) : arr[i]; if (!Number.isFinite(v)) continue; if (log && v <= 0) continue; if (v < lo) lo = v; if (v > hi) hi = v; }
    return [lo, hi];
  }
  function globalDomain(kind) {
    const log = state.yscale === 'log';
    let lo = Infinity, hi = -Infinity;
    for (let k = 0; k < NSEG; k++) { const e = yExtent(yValues(grid(k), kind), log); if (e[0] < lo) lo = e[0]; if (e[1] > hi) hi = e[1]; }
    return [lo, hi];
  }
  function padDomain(d, log) {
    if (!(d[0] < d[1])) return log ? [1e-9, 1e-6] : [0, 1];
    if (log) return [Math.max(d[0], 1e-10) / 1.5, d[1] * 1.5];
    const span = d[1] - d[0] || Math.abs(d[1]) || 1;
    return [d[0] - 0.06 * span, d[1] + 0.08 * span];
  }
  function targetFrame(k) {
    const g = grid(k), log = state.yscale === 'log';
    const S = yValues(g, 'smooth'), R = yValues(g, 'raw'), L = lowerValues(g);
    let dom;
    if (state.overlay === 'all') dom = globalDomain('smooth');
    else {
      dom = yExtent(S, log);
      if (state.layers.points) { const rr = yExtent(R, log); dom = [Math.min(dom[0], rr[0]), Math.max(dom[1], rr[1])]; }
      if (state.overlay === 'pair') { const pe = yExtent(yValues(grid(partnerOf(k)), 'smooth'), log); dom = [Math.min(dom[0], pe[0]), Math.max(dom[1], pe[1])]; }
      if (!log) dom = [Math.min(dom[0], 0), dom[1]];
    }
    dom = padDomain(dom, log);
    let dd;
    if (state.lower === 'telocal') { const e = yExtent(L, false); dd = [0, Number.isFinite(e[1]) ? Math.min(e[1], 12) * 1.1 : 8]; }
    else { const e = padDomain(yExtent(L, false), false); dd = [Math.min(e[0], 0), e[1]]; }
    return { k, V: g.V, S, R, L, dom, ddom: dd };
  }
  function lerpArr(a, b, t) {
    const out = new Float64Array(b.length);
    for (let i = 0; i < b.length; i++) { const x = a ? a[i] : NaN; out[i] = Number.isFinite(x) && Number.isFinite(b[i]) ? x + (b[i] - x) * t : b[i]; }
    return out;
  }
  const ease = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);
  function showSweep(k, animate) {
    const to = targetFrame(k), from = iv.cur;
    const dur = (!animate || reduceMotion || !from) ? 0 : (state.playing ? Math.min(380, 700 / state.speed) : 320);
    if (iv.tween) cancelAnimationFrame(iv.tween.raf);
    if (dur === 0) { iv.cur = to; drawIV(to, to.k); return; }
    const lg = state.yscale === 'log', t0 = performance.now();
    const step = now => {
      const t = Math.min(1, (now - t0) / dur), e = ease(t);
      const dom = lg ? [Math.exp(Math.log(from.dom[0]) + (Math.log(to.dom[0]) - Math.log(from.dom[0])) * e), Math.exp(Math.log(from.dom[1]) + (Math.log(to.dom[1]) - Math.log(from.dom[1])) * e)]
        : [from.dom[0] + (to.dom[0] - from.dom[0]) * e, from.dom[1] + (to.dom[1] - from.dom[1]) * e];
      const frame = { k: to.k, V: lerpArr(from.V, to.V, e), S: lerpArr(from.S, to.S, e), R: t < 1 ? lerpArr(from.R, to.R, e) : to.R, L: lerpArr(from.L, to.L, e), dom,
        ddom: [from.ddom[0] + (to.ddom[0] - from.ddom[0]) * e, from.ddom[1] + (to.ddom[1] - from.ddom[1]) * e] };
      drawIV(frame, to.k, t < 1);
      if (t < 1) iv.tween.raf = requestAnimationFrame(step); else { iv.cur = to; iv.tween = null; }
    };
    iv.tween = { raf: requestAnimationFrame(step) };
  }

  function drawIV(f, k, inFlight) {
    const w = iv.w, log = state.yscale === 'log', showD = state.layers.deriv, r = results[k], g = grid(k), teLocal = state.lower === 'telocal';
    const x = iv.x.domain([Math.min(r.Vmin, -10) - 0.2, Math.max(r.Vmax, 10) + 0.2].map(v => Number.isFinite(v) ? v : 0));
    const unit = log ? { s: 1e9, u: 'nA' } : unitFor(Math.max(Math.abs(f.dom[0]), Math.abs(f.dom[1])));
    const y = (log ? d3.scaleLog() : d3.scaleLinear()).domain(f.dom).range(iv.yMainRange);
    const dunit = teLocal ? { s: 1, u: 'eV' } : unitFor(Math.max(Math.abs(f.ddom[0]), Math.abs(f.ddom[1])));
    const yd = d3.scaleLinear().domain(f.ddom).range(iv.yLowerRange);
    const yv = v => log ? (Math.abs(v) > 0 ? y(Math.abs(v)) : NaN) : y(v);

    /* grid + axes */
    const xt = x.ticks(w < 520 ? 5 : 9), ydt = yd.ticks(3);
    let yt = y.ticks(6);
    if (log) { const decades = Math.log10(f.dom[1] / f.dom[0]); const pow10 = d => Math.abs(Math.log10(d) - Math.round(Math.log10(d))) < 1e-9; yt = y.ticks().filter(d => decades > 2.5 ? pow10(d) : (pow10(d) || /^[25]/.test(d.toExponential(0)))); }
    const bottomY = showD ? iv.yLowerRange[0] : iv.yMainRange[0];
    gGrid.selectAll('line.gx').data(xt).join('line').attr('class', 'grid-line gx').attr('x1', d => x(d)).attr('x2', d => x(d)).attr('y1', iv.m.t).attr('y2', bottomY);
    gGrid.selectAll('line.gy').data(yt).join('line').attr('class', 'grid-line gy').attr('x1', iv.m.l).attr('x2', w - iv.m.r).attr('y1', d => y(d)).attr('y2', d => y(d));
    gGrid.selectAll('line.gyd').data(showD ? ydt : []).join('line').attr('class', 'grid-line gyd').attr('x1', iv.m.l).attr('x2', w - iv.m.r).attr('y1', d => yd(d)).attr('y2', d => yd(d));
    gGrid.selectAll('line.zero').data(!log && f.dom[0] < 0 && f.dom[1] > 0 ? [0] : []).join('line').attr('class', 'zero-line zero').attr('x1', iv.m.l).attr('x2', w - iv.m.r).attr('y1', y(0)).attr('y2', y(0));
    gGrid.selectAll('line.zerod').data(showD && f.ddom[0] < 0 ? [0] : []).join('line').attr('class', 'zero-line zerod').attr('x1', iv.m.l).attr('x2', w - iv.m.r).attr('y1', yd(0)).attr('y2', yd(0));
    gAxes.selectAll('line.ax').data([0]).join('line').attr('class', 'axis-line ax').attr('x1', iv.m.l).attr('x2', w - iv.m.r).attr('y1', bottomY).attr('y2', bottomY);
    gAxes.selectAll('text.tx').data(xt).join('text').attr('class', 'tick-text tx').attr('x', d => x(d)).attr('y', bottomY + 14).attr('text-anchor', 'middle').text(d => (d < 0 ? '−' : d > 0 ? '+' : '') + Math.abs(d) + ' V');
    const decimalsFor = (ticks, scale) => { if (ticks.length < 2) return 1; const step = Math.abs(ticks[1] - ticks[0]) * scale; return Math.max(0, Math.min(4, Math.ceil(-Math.log10(step) - 1e-9))); };
    const yDec = decimalsFor(yt, unit.s), ydDec = decimalsFor(ydt, dunit.s);
    const fmtSI = d => { const u = unitFor(d), v = d * u.s; return (v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(0) : v.toPrecision(1)) + ' ' + u.u; };
    const fmtY = d => log ? fmtSI(d) : (d * unit.s).toFixed(yDec);
    gAxes.selectAll('text.ty').data(yt).join('text').attr('class', 'tick-text ty').attr('x', iv.m.l - 8).attr('y', d => y(d) + 3.5).attr('text-anchor', 'end').text(d => (d < 0 ? '−' : '') + fmtY(Math.abs(d)));
    gAxes.selectAll('text.tyd').data(showD ? ydt : []).join('text').attr('class', 'tick-text tyd').attr('x', iv.m.l - 8).attr('y', d => yd(d) + 3.5).attr('text-anchor', 'end').text(d => { const v = d * dunit.s; return (d < 0 ? '−' : '') + Math.abs(v).toFixed(ydDec); });
    const qLabel = state.quantity === 'Ie' ? 'Electron current Ie' : 'Probe current';
    gAxes.selectAll('text.yl').data([0]).join('text').attr('class', 'axis-title yl').attr('transform', `translate(14,${(iv.yMainRange[0] + iv.yMainRange[1]) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(log ? `|${qLabel}|` : `${qLabel} (${unit.u})`);
    gAxes.selectAll('text.ydl').data(showD ? [0] : []).join('text').attr('class', 'axis-title ydl').attr('transform', `translate(14,${(iv.yLowerRange[0] + iv.yLowerRange[1]) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(teLocal ? 'local Te (eV)' : `dI/dV (${dunit.u}/V)`);
    gAxes.selectAll('text.xl').data([0]).join('text').attr('class', 'axis-title xl').attr('x', (iv.m.l + w - iv.m.r) / 2).attr('y', bottomY + 28).attr('text-anchor', 'middle').text('Probe bias V₂ (V)');

    /* lines */
    const line = d3.line().defined(d => Number.isFinite(d[0]) && Number.isFinite(d[1])).x(d => x(d[0])).y(d => d[1]);
    const pts = (V, Y, fy) => { const out = new Array(V.length); for (let i = 0; i < V.length; i++) out[i] = [V[i], Number.isFinite(Y[i]) ? fy(Y[i]) : NaN]; return out; };
    pathMain.attr('d', state.layers.smooth ? line(pts(f.V, f.S, yv)) : null);
    if (state.overlay === 'pair') { const pg = grid(partnerOf(k)); pathPair.attr('d', line(pts(pg.V, yValues(pg, 'smooth'), yv))); } else pathPair.attr('d', null);
    pathLower.attr('class', teLocal ? 'line-telocal' : 'line-deriv').attr('d', showD ? line(pts(f.V, f.L, v => yd(v))) : null);

    /* fits, shading, markers, readout (target sweep, not tweened) */
    const analysed = r.Ie && !r.flags.includes('no_plasma_signal');
    const showF = state.layers.fits && analysed;
    gFits.style('display', showF ? null : 'none'); gShade.style('display', showF ? null : 'none'); gMark.style('display', state.layers.fits ? null : 'none'); gRead.style('display', state.layers.fits ? null : 'none');
    const readout = [];
    if (showF) {
      const vIon = r.Vmin + state.params.fIon * (r.Vmax - r.Vmin);
      shadeIon.attr('x', x(r.Vmin)).attr('y', iv.m.t).attr('width', Math.max(0, x(vIon) - x(r.Vmin))).attr('height', iv.mainH);
      const tw = r.Te_window, hasTw = Number.isFinite(tw[0]) && Number.isFinite(tw[1]);
      shadeTe.attr('x', hasTw ? x(tw[0]) : 0).attr('y', iv.m.t).attr('width', hasTw ? Math.max(0, x(tw[1]) - x(tw[0])) : 0).attr('height', iv.mainH);
      shadeIonD.attr('x', x(r.Vmin)).attr('y', iv.yLowerRange[1]).attr('width', showD ? Math.max(0, x(vIon) - x(r.Vmin)) : 0).attr('height', iv.derivH);
      shadeTeD.attr('x', hasTw ? x(tw[0]) : 0).attr('y', iv.yLowerRange[1]).attr('width', showD && hasTw ? Math.max(0, x(tw[1]) - x(tw[0])) : 0).attr('height', iv.derivH);
      if (state.quantity === 'I' && Number.isFinite(r.a_i)) {
        const vs = d3.range(r.Vmin, r.Vmax + 0.01, (r.Vmax - r.Vmin) / 60);
        pathIon.attr('d', line(vs.map(v => [v, yv(r.a_i + r.b_i * v)])));
      } else pathIon.attr('d', null);
      if (r.Te !== null && hasTw) {
        const lo = tw[0] - 0.8, hi = Math.min(tw[1] + 1.2, r.Vp !== null ? r.Vp + 0.4 : r.Vmax);
        const vs = d3.range(lo, hi + 0.001, (hi - lo) / 80);
        pathTe.attr('d', line(vs.map(v => { const ie = Math.exp(r.te_alpha + r.te_beta * v); return [v, yv(state.quantity === 'I' ? ie + r.a_i + r.b_i * v : ie)]; })));
      } else pathTe.attr('d', null);
      if (r.Te !== null) readout.push({ t: `Te = ${r.Te.toFixed(2)} eV  R² ${Number.isFinite(r.Te_r2) ? r.Te_r2.toFixed(3) : '—'}  window ${fmtV(tw[0], 1)} → ${fmtV(tw[1], 1)}`, strong: true });
      if (r.Te_lo !== null && r.Te_hi !== null) readout.push({ t: `lower / upper half of window: ${r.Te_lo.toFixed(2)} / ${r.Te_hi.toFixed(2)} eV` });
      if (Number.isFinite(r.a_i)) readout.push({ t: `ion fit  I = ${fmtAstr(r.a_i)} ${r.b_i < 0 ? '−' : '+'} ${fmtAstr(Math.abs(r.b_i))}/V · V` });
      if (r.Vp === null) readout.push({ t: `Vp beyond +10 V` + (r.Vp_expected !== null ? `  (sheath estimate ${fmtV(r.Vp_expected, 1)})` : '') });
      if (r.hysteresis !== null) readout.push({ t: `up/down hysteresis ${r.hysteresis >= 0 ? '+' : '−'}${Math.abs(r.hysteresis).toFixed(2)} V` });
    } else if (state.layers.fits) {
      if (r.flags.includes('no_plasma_signal')) readout.push({ t: 'no plasma signal above noise', strong: true });
      else if (r.flags.includes('partial_sweep')) readout.push({ t: 'partial sweep: not analysed', strong: true });
    }
    // markers: Vf, Vp
    const marks = [];
    if (analysed && r.Vf !== null) marks.push({ id: 'vf', v: r.Vf, label: 'Vf', sub: fmtV(r.Vf) + (r.dVf !== null ? ` ± ${r.dVf.toFixed(2)}` : '') });
    if (analysed && r.Vp !== null) marks.push({ id: 'vp', v: r.Vp, label: 'Vp', sub: fmtV(r.Vp) });
    const mg = gMark.selectAll('g.mk').data(marks, d => d.id).join(enter => {
      const gg = enter.append('g').attr('class', 'mk');
      gg.append('line').attr('class', 'marker-line'); gg.append('text').attr('class', 'marker-text'); gg.append('text').attr('class', 'marker-sub'); return gg;
    });
    mg.attr('transform', d => `translate(${x(d.v)},0)`);
    mg.select('line').attr('y1', iv.m.t).attr('y2', bottomY);
    const flip = d => x(d.v) > w - iv.m.r - 110;
    mg.select('.marker-text').attr('x', d => flip(d) ? -4 : 4).attr('text-anchor', d => flip(d) ? 'end' : 'start').attr('y', iv.yMainRange[0] - 26).text(d => d.label);
    mg.select('.marker-sub').attr('x', d => flip(d) ? -4 : 4).attr('text-anchor', d => flip(d) ? 'end' : 'start').attr('y', iv.yMainRange[0] - 14).text(d => d.sub);
    // readout block (top-left of the main panel)
    const rg = gRead.selectAll('g.rd').data([0]).join(enter => { const gg = enter.append('g').attr('class', 'rd'); gg.append('rect').attr('class', 'readout-bg').attr('rx', 4); return gg; });
    const lines = rg.selectAll('text.readout').data(readout).join('text').attr('class', d => 'readout' + (d.strong ? ' strong' : '')).attr('x', iv.m.l + 10).attr('y', (d, i) => iv.m.t + 16 + i * 14).text(d => d.t);
    const maxLen = readout.reduce((m, d) => Math.max(m, d.t.length), 0);
    rg.select('rect').attr('x', iv.m.l + 4).attr('y', iv.m.t + 4).attr('width', readout.length ? Math.min(w - iv.m.l - iv.m.r - 8, maxLen * 6.4 + 12) : 0).attr('height', readout.length ? readout.length * 14 + 8 : 0);

    /* canvas: ghosts + raw samples + excluded */
    const dpr = window.devicePixelRatio || 1;
    ctx.save(); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, iv.h);
    if (state.overlay === 'all') {
      const key = [w, iv.h, state.yscale, state.quantity, f.dom.join(','), state.channel, JSON.stringify(state.baseline), state.params.window, state.params.fIon].join('|');
      if (key !== ghostKey || !ghostCanvas) {
        ghostCanvas = document.createElement('canvas'); ghostCanvas.width = el.ivCanvas.width; ghostCanvas.height = el.ivCanvas.height;
        const g2 = ghostCanvas.getContext('2d'); g2.setTransform(dpr, 0, 0, dpr, 0, 0);
        g2.strokeStyle = tok('--ghost'); g2.globalAlpha = 0.45; g2.lineWidth = 1;
        g2.beginPath(); g2.rect(iv.m.l, iv.m.t, w - iv.m.l - iv.m.r, iv.mainH); g2.clip();
        for (let j = 0; j < NSEG; j++) {
          const gj = grid(j), Y = yValues(gj, 'smooth'); let pen = false; g2.beginPath();
          for (let i = 0; i < NG; i++) { const yy = yv(Y[i]); if (!Number.isFinite(gj.V[i]) || !Number.isFinite(yy)) { pen = false; continue; } const xx = x(gj.V[i]); if (pen) g2.lineTo(xx, yy); else { g2.moveTo(xx, yy); pen = true; } }
          g2.stroke();
        }
        ghostKey = key;
      }
      ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.drawImage(ghostCanvas, 0, 0); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    ctx.beginPath(); ctx.rect(iv.m.l, iv.m.t, w - iv.m.l - iv.m.r, iv.mainH); ctx.clip();
    if (state.layers.points) {
      ctx.fillStyle = hexA(tok('--s1'), inFlight ? 0.25 : 0.42);
      const R = f.R;
      for (let i = 0; i < NG; i++) { const yy = yv(R[i]); if (!Number.isFinite(f.V[i]) || !Number.isFinite(yy)) continue; ctx.beginPath(); ctx.arc(x(f.V[i]), yy, 2, 0, 6.2832); ctx.fill(); }
      if (g.ex.length && state.quantity === 'I') {
        ctx.strokeStyle = tok('--ghost'); ctx.lineWidth = 1;
        for (const d of g.ex) { const yy = yv(d[1]); if (!Number.isFinite(yy)) continue; ctx.beginPath(); ctx.arc(x(d[0]), Math.max(iv.m.t, Math.min(iv.yMainRange[0], yy)), 2.5, 0, 6.2832); ctx.stroke(); }
      }
    }
    ctx.restore();
    el.lgExcl.hidden = !(g.ex.length && state.layers.points && state.quantity === 'I');

    iv.scales = { x, y, yd, yv, unit, dunit, showD };
    if (hoverG >= 0) hoverAt(hoverG);
  }

  /* hover / tooltip */
  let hoverG = -1;
  function hoverAt(g) {
    const s = iv.scales; if (!s || !iv.cur) return;
    const gr = grid(state.sweep);
    if (g < 0 || !Number.isFinite(gr.V[g])) { crossX.style('opacity', 0); hoverDot.style('opacity', 0); hoverDotD.style('opacity', 0); el.ivTip.classList.remove('on'); return; }
    const xx = s.x(gr.V[g]);
    crossX.attr('x1', xx).attr('x2', xx).attr('y1', iv.m.t).attr('y2', s.showD ? iv.yLowerRange[0] : iv.yMainRange[0]).style('opacity', 1);
    const ys = yValues(gr, 'smooth')[g], yy = s.yv(ys);
    hoverDot.attr('cx', xx).attr('cy', yy).style('opacity', Number.isFinite(yy) && state.layers.smooth ? 1 : 0);
    const lv = lowerValues(gr)[g], dy = s.yd(lv); hoverDotD.attr('cx', xx).attr('cy', dy).style('fill', state.lower === 'telocal' ? 'var(--s4)' : 'var(--s2)').style('opacity', s.showD && Number.isFinite(dy) ? 1 : 0);
    const rows = [
      ['I measured', fmtAstr(gr.R[g]), 'var(--s1)'],
      ['I smoothed', fmtAstr(gr.S[g]), 'var(--s1)'],
    ];
    if (state.overlay === 'pair') { const q = partnerOf(state.sweep); rows.push([`I sweep ${q + 1} (${results[q].dir})`, fmtAstr(grid(q).S[g]), 'var(--s2)']); }
    if (Number.isFinite(gr.ES[g])) rows.push(['Ie = I − Iion', fmtAstr(gr.ES[g]), 'var(--s4)']);
    rows.push(['dI/dV', Number.isFinite(gr.D[g]) ? fmtAstr(gr.D[g]) + '/V' : '—', 'var(--s2)']);
    if (Number.isFinite(gr.TL[g])) rows.push(['local Te', gr.TL[g].toFixed(2) + ' eV', 'var(--s4)']);
    const tip = el.ivTip; tip.innerHTML = '';
    const h = document.createElement('div'); h.className = 't-head'; h.textContent = `V = ${fmtV(gr.V[g], 3)} · sweep ${state.sweep + 1}`; tip.appendChild(h);
    for (const [lab, val, col] of rows) {
      const row = document.createElement('div'); row.className = 't-row';
      const l = document.createElement('span'); const i = document.createElement('i'); i.style.setProperty('--k', col); l.appendChild(i); l.appendChild(document.createTextNode(lab));
      const b = document.createElement('b'); b.textContent = val; row.appendChild(l); row.appendChild(b); tip.appendChild(row);
    }
    tip.classList.add('on');
    const tw = tip.offsetWidth || 160; const left = xx + 14 + tw > iv.w ? xx - 14 - tw : xx + 14;
    tip.style.left = left + 'px'; tip.style.top = Math.max(0, Math.min(iv.h - 140, (Number.isFinite(yy) ? yy : iv.m.t + 40) - 30)) + 'px';
  }
  function nearestG(V, v) {
    let best = -1, bd = Infinity;
    for (let g = 0; g < NG; g++) { if (!Number.isFinite(V[g])) continue; const d = Math.abs(V[g] - v); if (d < bd) { bd = d; best = g; } }
    return best;
  }
  hit.on('pointermove', ev => { const [mx] = d3.pointer(ev, el.ivSvg); hoverG = nearestG(grid(state.sweep).V, iv.scales.x.invert(mx)); hoverAt(hoverG); })
    .on('pointerleave', () => { hoverG = -1; hoverAt(-1); });

  /* ======================= strip (whole log) ======================= */
  const sctx = el.strip.getContext('2d');
  let stripCache = null;
  function drawStrip() {
    const wrap = el.strip.parentElement, w = wrap.clientWidth, h = wrap.clientHeight, dpr = window.devicePixelRatio || 1;
    if (!w || !h) return;
    el.strip.width = Math.round(w * dpr); el.strip.height = Math.round(h * dpr); sctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const ts = DATA.cols.ts, I = DATA.cols[state.channel], n = ts.length, T = ts[n - 1];
    if (!stripCache || stripCache.w !== w || stripCache.ch !== state.channel) {
      const vals = []; for (let i = 0; i < n; i += 7) if (Number.isFinite(I[i])) vals.push(I[i]);
      vals.sort((a, b) => a - b);
      const lo = Math.min(0, vals[Math.floor(vals.length * 0.002)]), hi = vals[Math.floor(vals.length * 0.998)] * 1.05;
      const mins = new Float64Array(w).fill(Infinity), maxs = new Float64Array(w).fill(-Infinity);
      for (let i = 0; i < n; i++) { if (!Number.isFinite(I[i])) continue; const px = Math.min(w - 1, Math.floor(ts[i] / T * w)); if (I[i] < mins[px]) mins[px] = I[i]; if (I[i] > maxs[px]) maxs[px] = I[i]; }
      stripCache = { w, ch: state.channel, lo, hi, mins, maxs };
    }
    const c = stripCache, ys = v => h - 4 - (Math.min(c.hi, Math.max(c.lo, v)) - c.lo) / (c.hi - c.lo) * (h - 8);
    sctx.clearRect(0, 0, w, h);
    const k = state.sweep, b = DATA.bounds, tx = t => t / T * w;
    // plasma-off reference span
    if (state.baseline.mode === 'ref' && refLine) { sctx.fillStyle = hexA(tok('--s3'), 0.10); sctx.fillRect(tx(ts[b[refLine.from]]), 0, Math.max(2, tx(ts[b[refLine.to + 1] - 1]) - tx(ts[b[refLine.from]])), h); }
    if (state.overlay === 'pair') { const p = partnerOf(k); sctx.fillStyle = hexA(tok('--s2'), 0.14); sctx.fillRect(tx(ts[b[p]]), 0, Math.max(2, tx(ts[b[p + 1] - 1]) - tx(ts[b[p]])), h); }
    sctx.fillStyle = tok('--accent-soft'); sctx.fillRect(tx(ts[b[k]]), 0, Math.max(2, tx(ts[b[k + 1] - 1]) - tx(ts[b[k]])), h);
    sctx.strokeStyle = tok('--grid'); sctx.lineWidth = 1; sctx.beginPath();
    for (let j = 1; j < NSEG; j++) { const xx = Math.round(tx(ts[b[j]])) + 0.5; sctx.moveTo(xx, h - 6); sctx.lineTo(xx, h); }
    sctx.stroke();
    sctx.strokeStyle = tok('--s1'); sctx.lineWidth = 1; sctx.beginPath();
    for (let px = 0; px < w; px++) { if (c.mins[px] === Infinity) continue; sctx.moveTo(px + 0.5, ys(c.maxs[px])); sctx.lineTo(px + 0.5, ys(c.mins[px]) + 0.75); }
    sctx.stroke();
    sctx.strokeStyle = tok('--accent'); sctx.lineWidth = 1.5; sctx.beginPath();
    const x0 = Math.floor(tx(ts[b[k]])), x1 = Math.ceil(tx(ts[b[k + 1] - 1]));
    for (let px = x0; px <= Math.min(w - 1, x1); px++) { if (c.mins[px] === Infinity) continue; sctx.moveTo(px + 0.5, ys(c.maxs[px])); sctx.lineTo(px + 0.5, ys(c.mins[px]) + 0.75); }
    sctx.stroke();
    el.stripLeft.textContent = `whole log · ${state.channel} vs time · ${NSEG} half-sweeps` + (state.baseline.mode === 'ref' && refLine ? ` · green = plasma-off reference` : '');
    el.stripRight.textContent = `${clockAt(0)} → ${clockAt(T)} UTC`;
  }
  function sweepAtTime(t) { const ts = DATA.cols.ts, b = DATA.bounds; for (let k = 0; k < NSEG; k++) if (t < ts[b[k + 1] - 1]) return k; return LAST; }
  el.strip.addEventListener('click', ev => { const rect = el.strip.getBoundingClientRect(); setSweep(sweepAtTime((ev.clientX - rect.left) / rect.width * DATA.cols.ts[DATA.cols.ts.length - 1]), true); });

  /* ======================= tiles & flags ======================= */
  function tile(l, v, u, n, dim) { return { l, v, u, n, dim }; }
  function renderTiles() {
    const r = results[state.sweep], gas = GAS_LABEL[String(+state.params.mass)] || state.params.mass + ' u';
    el.anSub.textContent = `sweep ${state.sweep + 1} · ${r.dir === 'up' ? 'up (−10 → +10 V)' : 'down (+10 → −10 V)'} · ${fmtT(r.tStart)}`;
    const gated = !r.Ie || r.flags.includes('no_plasma_signal');
    const tiles = [];
    tiles.push(tile('Electron temperature (slope)', r.Te !== null ? r.Te.toFixed(2) : '—', 'eV',
      r.Te !== null ? `R² ${r.Te_r2.toFixed(3)} · ${r.Te_npts} pts` + (r.Te_lo !== null && r.Te_hi !== null ? ` · halves ${r.Te_lo.toFixed(1)} / ${r.Te_hi.toFixed(1)} eV` : '') : (gated ? 'not analysed' : 'no fit'), r.Te === null));
    tiles.push(tile('Floating potential', r.Vf !== null ? fmtV(r.Vf).replace(' V', '') : (gated ? '—' : (r.flags.includes('no_zero_crossing') ? 'n/a' : '—')), 'V',
      r.Vf !== null ? `± ${r.dVf !== null ? r.dVf.toFixed(2) + ' V' : '—'} · zero crossing` : (gated ? 'not analysed' : 'no zero crossing in ±10 V'), r.Vf === null));
    tiles.push(tile('Plasma potential', r.Vp !== null ? fmtV(r.Vp).replace(' V', '') : (gated ? '—' : '> +10'), 'V',
      r.Vp !== null ? 'max of dI/dV' : (r.flags.includes('top_of_sweep_unreliable') ? 'top of sweep lost' : (r.Vp_expected !== null ? `no knee · Vf + 5.2 Te = ${fmtV(r.Vp_expected, 1)}` : (gated ? 'not analysed' : 'no knee inside the sweep'))), r.Vp === null));
    const ies = fmtA(r.Ies); tiles.push(tile(r.Vp !== null ? 'Electron saturation current' : 'Peak electron current', ies ? ies.v : '—', ies ? ies.u : '', r.Vp !== null ? 'Ie at Vp' : (r.flags.includes('top_of_sweep_unreliable') ? 'top of sweep lost' : (gated ? 'not analysed' : 'max Ie in sweep (near +10 V) · lower bound for Ies')), r.Ies === null));
    const iis = fmtA(r.Iis); tiles.push(tile('Ion current at −10 V', iis ? iis.v : '—', iis ? iis.u : '',
      r.Iis === null ? (gated ? 'not analysed' : '—') : (r.flags.includes('iis_below_noise') ? `below 3σ noise (σ = ${fmtAstr(r.sigma_I)})` : (r.Iis >= 0 ? 'positive: not an ion current' : (r.Iis_vf !== null ? `lowest 1 V mean · fit at Vf: ${fmtAstr(r.Iis_vf)}` : 'mean over lowest 1 V'))), !(r.Iis < 0)));
    tiles.push(tile('Up/down hysteresis', r.hysteresis !== null ? (r.hysteresis >= 0 ? '+' : '−') + Math.abs(r.hysteresis).toFixed(2) : '—', 'V', r.hysteresis !== null ? `pair ${Math.min(state.sweep, partnerOf(state.sweep)) + 1}+${Math.max(state.sweep, partnerOf(state.sweep)) + 1}: down-sweep ${r.hysteresis >= 0 ? 'lags' : 'leads'} by this much` : (r.flags.includes('hysteresis_unresolved') ? 'pair could not be aligned' : 'needs a clean up + down pair'), r.hysteresis === null));
    tiles.push(Object.assign(tile('Electron density n\u2091' + (r.Vp === null ? ' (lower bound)' : ''), fmtSci(r.ne), r.ne ? 'm⁻³' : '', r.ne ? `λD = ${fmtLen(r.lambdaD)} · A = ${(state.params.area * 1e6).toPrecision(3)} mm²` : 'needs Te and Ies > 0', !r.ne), { wide: true }));
    tiles.push(Object.assign(tile('Ion density n\u1d62 (Bohm, ' + gas + ')', fmtSci(r.ni), r.ni ? 'm⁻³' : '', r.ni ? (r.ne ? `n\u1d62/n\u2091 = ${(r.ni / r.ne).toPrecision(2)} (n\u2091 is a lower bound)` : '') : 'needs Te and an ion current above noise', !r.ni), { wide: true }));
    el.tiles.innerHTML = '';
    for (const t of tiles) {
      const d = document.createElement('div'); d.className = 'tile' + (t.dim ? ' dim' : '') + (t.wide ? ' wide' : '');
      const l = document.createElement('div'); l.className = 'l'; l.textContent = t.l;
      const v = document.createElement('div'); v.className = 'v'; v.textContent = t.v; if (t.u) { const s = document.createElement('small'); s.textContent = t.u; v.appendChild(s); }
      const n = document.createElement('div'); n.className = 'n'; n.textContent = t.n; n.title = t.n;
      d.appendChild(l); d.appendChild(v); d.appendChild(n); el.tiles.appendChild(d);
    }
    el.flags.innerHTML = '';
    const info = document.createElement('span'); info.className = 'flag ok';
    info.textContent = `${r.n} samples · ${r.duration.toFixed(1)} s` + (r.nGlitch ? ` · ${r.nGlitch} excluded` : '') + (r.nDropped - r.nGlitch > 0 ? ` · ${r.nDropped - r.nGlitch} NaN` : '') + (Number.isFinite(r.sigma_I) ? ` · noise σ ${fmtAstr(r.sigma_I)}` : '');
    el.flags.appendChild(info);
    for (const f of r.flags) { const t = FLAG_TEXT[f] || [f, 'warn']; const s = document.createElement('span'); s.className = 'flag ' + t[1]; s.textContent = t[0]; el.flags.appendChild(s); }
  }

  /* ======================= trends (small multiples) ======================= */
  const METRICS = [
    { id: 'Te', label: 'Electron temperature', get: r => r.Te, fmt: v => v.toFixed(2) + ' eV', tick: d => d },
    { id: 'Vf', label: 'Floating potential', get: r => r.Vf, fmt: v => fmtV(v), tick: d => (d < 0 ? '−' : '') + Math.abs(d) },
    { id: 'Ihi', label: 'Current at +10 V', get: r => (r.Ie && !r.flags.includes('no_plasma_signal') && !r.flags.includes('top_of_sweep_unreliable')) ? r.I_hi : null, fmt: v => fmtAstr(v), tick: d => fmtAstr(d, 0) },
    { id: 'Iis', label: 'Ion current at −10 V', get: r => r.Iis, fmt: v => fmtAstr(v), tick: d => fmtAstr(d, 0) },
    { id: 'ne', label: 'Electron density (lower bound)', get: r => r.ne, fmt: v => fmtSci(v) + ' m⁻³', tick: d => d === 0 ? '0' : d3.format('.1~e')(d).replace('e+', 'e') },
    { id: 'hys', label: 'Up/down hysteresis', get: r => r.hysteresis, fmt: v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(2) + ' V', tick: d => (d < 0 ? '−' : '') + Math.abs(d) },
  ];
  const minis = {};
  function buildMinis() {
    el.multiples.innerHTML = '';
    for (const m of METRICS) {
      const card = document.createElement('div'); card.className = 'card mini';
      const t = document.createElement('div'); t.className = 't';
      const b = document.createElement('b'); b.textContent = m.label; const cur = document.createElement('span'); cur.className = 'cur';
      t.appendChild(b); t.appendChild(cur); card.appendChild(t);
      const svgm = d3.select(card).append('svg').attr('role', 'img').attr('aria-label', m.label + ' per sweep');
      const tip = document.createElement('div'); tip.className = 'tip'; card.appendChild(tip); card.style.position = 'relative';
      el.multiples.appendChild(card);
      minis[m.id] = { m, card, cur, svg: svgm, tip, g: svgm.append('g') };
    }
  }
  function drawMinis() {
    const ts = DATA.cols.ts, T = ts[ts.length - 1];
    for (const id in minis) {
      const mm = minis[id], m = mm.m, w = Math.max(160, mm.card.clientWidth - 24), h = 110, ml = 44, mr = 10, mt = 8, mb = 18;
      mm.svg.attr('viewBox', `0 0 ${w} ${h}`).attr('width', w).attr('height', h);
      const pts = results.map(r => ({ k: r.k, t: r.tStart, v: m.get(r), dir: r.dir })).filter(d => d.v !== null && Number.isFinite(d.v));
      const x = d3.scaleLinear().domain([0, T / 60]).range([ml, w - mr]);
      const ext = d3.extent(pts, d => d.v); const lo = Math.min(0, ext[0] || 0), hi = Math.max(0, ext[1] || 0);
      const y = d3.scaleLinear().domain([lo, hi === lo ? lo + 1 : hi]).nice(4).range([h - mb, mt]);
      const g = mm.g; g.selectAll('*').remove();
      const yt = y.ticks(3);
      g.selectAll('line.gy').data(yt).join('line').attr('class', 'grid-line').attr('x1', ml).attr('x2', w - mr).attr('y1', d => y(d)).attr('y2', d => y(d));
      g.selectAll('text.ty').data(yt).join('text').attr('class', 'tick-text').attr('x', ml - 6).attr('y', d => y(d) + 3.5).attr('text-anchor', 'end').text(d => m.tick(d));
      g.append('line').attr('class', 'axis-line').attr('x1', ml).attr('x2', w - mr).attr('y1', h - mb).attr('y2', h - mb);
      const xt = x.ticks(w < 220 ? 2 : w < 340 ? 3 : 5).filter(d => d > 0);
      g.selectAll('text.tx').data(xt).join('text').attr('class', 'tick-text').attr('x', d => x(d)).attr('y', h - mb + 12).attr('text-anchor', 'middle').text(d => d + ' min');
      const line = d3.line().x(d => x(d.t / 60)).y(d => y(d.v));
      const runs = []; let run = [];
      for (const d of pts) { if (run.length && d.k !== run[run.length - 1].k + 1) { runs.push(run); run = []; } run.push(d); }
      if (run.length) runs.push(run);
      g.selectAll('path.line').data(runs).join('path').attr('class', 'line').attr('d', d => line(d));
      g.selectAll('circle.pt').data(pts, d => d.k).join('circle').attr('class', d => 'pt' + (d.dir === 'down' ? ' dn' : '') + (d.k === state.sweep ? ' cur' : '')).attr('cx', d => x(d.t / 60)).attr('cy', d => y(d.v)).attr('r', d => d.k === state.sweep ? 5 : 3);
      g.append('line').attr('class', 'cursor').attr('x1', x(results[state.sweep].tStart / 60)).attr('x2', x(results[state.sweep].tStart / 60)).attr('y1', mt).attr('y2', h - mb).style('opacity', 0.6);
      const curV = m.get(results[state.sweep]); mm.cur.textContent = curV !== null && Number.isFinite(curV) ? m.fmt(curV) : '—';
      const hitR = g.append('rect').attr('class', 'hit').attr('x', ml).attr('y', 0).attr('width', w - ml - mr).attr('height', h);
      const nearest = mx => { const t = x.invert(mx) * 60; let best = 0, bd = Infinity; for (const r of results) { const d = Math.abs(r.tStart + r.duration / 2 - t); if (d < bd) { bd = d; best = r.k; } } return best; };
      hitR.on('pointermove', ev => {
        const [mx] = d3.pointer(ev); const k = nearest(mx); const r = results[k]; const v = m.get(r);
        mm.tip.innerHTML = ''; const hd = document.createElement('div'); hd.className = 't-head'; hd.textContent = `sweep ${k + 1} · ${r.dir} · ${fmtT(r.tStart)}`; mm.tip.appendChild(hd);
        const row = document.createElement('div'); row.className = 't-row'; const s = document.createElement('span'); s.textContent = m.label; const b = document.createElement('b'); b.textContent = v !== null && Number.isFinite(v) ? m.fmt(v) : '—'; row.appendChild(s); row.appendChild(b); mm.tip.appendChild(row);
        mm.tip.classList.add('on'); const tw = mm.tip.offsetWidth || 150; mm.tip.style.left = Math.max(0, Math.min(mx + 12, w - tw)) + 'px'; mm.tip.style.top = '18px';
      }).on('pointerleave', () => mm.tip.classList.remove('on')).on('click', ev => { const [mx] = d3.pointer(ev); setSweep(nearest(mx), true); });
    }
  }

  /* ======================= results table ======================= */
  const COLS = [
    ['#', r => r.k + 1], ['dir', r => r.dir], ['t start', r => fmtT(r.tStart)], ['Vf', r => r.Vf !== null ? fmtV(r.Vf) : '—'], ['Vp', r => r.Vp !== null ? fmtV(r.Vp) : (r.Te !== null ? '> +10 V' : '—')],
    ['Te (eV)', r => r.Te !== null ? r.Te.toFixed(2) : '—'], ['R²', r => Number.isFinite(r.Te_r2) ? r.Te_r2.toFixed(3) : '—'], ['Te lo/hi', r => r.Te_lo !== null && r.Te_hi !== null ? `${r.Te_lo.toFixed(2)} / ${r.Te_hi.toFixed(2)}` : '—'],
    ['Ie(+10 V)', r => fmtAstr(r.Ies)], ['Ii(−10 V)', r => fmtAstr(r.Iis)], ['ne (m⁻³)', r => fmtSci(r.ne)], ['ni (m⁻³)', r => fmtSci(r.ni)], ['λD', r => fmtLen(r.lambdaD)],
    ['hyst.', r => r.hysteresis !== null ? fmtV(r.hysteresis) : '—'], ['σ noise', r => fmtAstr(r.sigma_I)], ['flags', r => r.flags.join(' ')],
  ];
  function renderTable() {
    if (!state.tableOpen) return;
    const thead = el.tbl.tHead, tbody = el.tbl.tBodies[0];
    thead.innerHTML = ''; const trh = document.createElement('tr'); for (const c of COLS) { const th = document.createElement('th'); th.textContent = c[0]; trh.appendChild(th); } thead.appendChild(trh);
    tbody.innerHTML = '';
    for (const r of results) {
      const tr = document.createElement('tr'); if (r.k === state.sweep) tr.className = 'cur';
      for (const c of COLS) { const td = document.createElement('td'); td.textContent = String(c[1](r)); tr.appendChild(td); }
      tr.addEventListener('click', () => setSweep(r.k, true));
      tbody.appendChild(tr);
    }
  }
  function tableCSV() {
    const head = ['sweep', 'dir', 't_start_s', 'duration_s', 'n', 'n_excluded', 'sigma_I_A', 'Vf_V', 'dVf_V', 'Vp_V', 'Vp_sheath_estimate_V', 'Te_eV', 'Te_r2', 'Te_window_lo_V', 'Te_window_hi_V', 'Te_lower_half_eV', 'Te_upper_half_eV', 'Ies_A', 'Iis_A', 'Iis_fit_at_Vf_A', 'ne_m3', 'ni_m3', 'lambdaD_m', 'I_lo_A', 'I_hi_A', 'dIdV_max_A_per_V', 'hysteresis_V', 'flags'];
    const num = v => (v === null || v === undefined || !Number.isFinite(v)) ? '' : String(v);
    const lines = ['# Langmuir Sweep Bench export · channel ' + state.channel + ' · baseline ' + JSON.stringify(state.baseline) + ' · params ' + JSON.stringify(state.params), head.join(',')];
    for (const r of results) lines.push([r.k + 1, r.dir, num(r.tStart), num(r.duration), r.n, r.nGlitch, num(r.sigma_I), num(r.Vf), num(r.dVf), num(r.Vp), num(r.Vp_expected), num(r.Te), num(r.Te_r2), num(r.Te_window[0]), num(r.Te_window[1]), num(r.Te_lo), num(r.Te_hi), num(r.Ies), num(r.Iis), num(r.Iis_vf), num(r.ne), num(r.ni), num(r.lambdaD), num(r.I_lo), num(r.I_hi), num(r.dIdV_max), num(r.hysteresis), '"' + r.flags.join(' ') + '"'].join(','));
    return lines.join('\n');
  }
  el.tblToggle.addEventListener('click', () => { state.tableOpen = !state.tableOpen; syncTable(); persist(); });
  function syncTable() {
    el.tblWrap.hidden = !state.tableOpen; el.tblCopy.hidden = !state.tableOpen;
    el.tblToggle.setAttribute('aria-expanded', String(state.tableOpen)); el.tblToggle.textContent = state.tableOpen ? 'Hide table' : 'Show table';
    renderTable();
  }
  el.tblCopy.addEventListener('click', () => {
    const csv = tableCSV(), done = () => { el.tblCopy.textContent = 'Copied'; setTimeout(() => { el.tblCopy.textContent = 'Copy as CSV'; }, 1500); };
    const fallback = () => { const ta = document.createElement('textarea'); ta.value = csv; ta.style.position = 'fixed'; ta.style.opacity = '0'; document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); done(); } catch (e) { el.tblCopy.textContent = 'Select the table and copy'; } document.body.removeChild(ta); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(csv).then(done, fallback); else fallback();
  });

  /* ======================= sweep selection & playback ======================= */
  function setSweep(k, animate) {
    k = Math.max(0, Math.min(LAST, k | 0));
    state.sweep = k; el.sweep.value = String(k); el.sweep.style.setProperty('--pct', (k / LAST * 100) + '%');
    const r = results[k];
    el.sweepLabel.textContent = `Sweep ${k + 1} of ${NSEG} · ${r.dir === 'up' ? '−10 → +10 V' : '+10 → −10 V'}`;
    el.sweepTime.textContent = `t = ${fmtT(r.tStart)} · ${clockAt(r.tStart)} UTC`;
    showSweep(k, animate !== false);
    renderTiles(); drawStrip(); drawMinis(); renderTable(); persist();
  }
  el.sweep.addEventListener('input', () => setSweep(+el.sweep.value, true));
  let playTimer = null, lastTick = 0;
  function tick(now) {
    if (!state.playing) return;
    if (now - lastTick >= 900 / state.speed) {
      lastTick = now;
      let next = state.sweep + 1;
      if (next > LAST) { if (state.loop) next = 0; else { setPlaying(false); return; } }
      setSweep(next, true);
    }
    playTimer = requestAnimationFrame(tick);
  }
  function setPlaying(on) {
    state.playing = on; el.play.setAttribute('aria-pressed', String(on)); el.play.setAttribute('aria-label', on ? 'Pause' : 'Play through sweeps');
    if (on) { lastTick = performance.now() - 1e9; playTimer = requestAnimationFrame(tick); } else if (playTimer) cancelAnimationFrame(playTimer);
  }
  el.play.addEventListener('click', () => setPlaying(!state.playing));
  segControl('speed', () => state.speed, v => { state.speed = +v; persist(); });
  el.loop.checked = state.loop; el.loop.addEventListener('change', () => { state.loop = el.loop.checked; persist(); });
  document.addEventListener('keydown', ev => {
    if (ev.target && /^(INPUT|SELECT|TEXTAREA|BUTTON)$/.test(ev.target.tagName) && ev.target.type !== 'range') return;
    if (ev.key === 'ArrowRight') { ev.preventDefault(); setSweep(state.sweep + 1, true); }
    else if (ev.key === 'ArrowLeft') { ev.preventDefault(); setSweep(state.sweep - 1, true); }
    else if (ev.key === ' ' && ev.target.tagName !== 'INPUT') { ev.preventDefault(); setPlaying(!state.playing); }
  });

  /* ======================= option controls ======================= */
  function refreshAll(recompute) {
    if (recompute) { recomputeAll(); ghostKey = ''; syncBaselineInputs(); }
    iv.cur = null; layout(); setSweep(state.sweep, false);
  }
  segControl('channel', () => state.channel, v => { state.channel = v; stripCache = null; persist(); refreshAll(true); });
  segControl('overlay', () => state.overlay, v => { state.overlay = v; el.lgPair.hidden = v !== 'pair'; el.lgGhost.hidden = v !== 'all'; persist(); refreshAll(false); });
  el.lgPair.hidden = state.overlay !== 'pair'; el.lgGhost.hidden = state.overlay !== 'all';
  segControl('yscale', () => state.yscale, v => { state.yscale = v; persist(); refreshAll(false); });
  segControl('quantity', () => state.quantity, v => { state.quantity = v; persist(); refreshAll(false); });
  segControl('lower', () => state.lower, v => { state.lower = v; persist(); refreshAll(false); });
  for (const [id, key] of [['ly-points', 'points'], ['ly-smooth', 'smooth'], ['ly-fits', 'fits'], ['ly-deriv', 'deriv']]) {
    const cb = $('#' + id); cb.checked = state.layers[key];
    cb.addEventListener('change', () => { state.layers[key] = cb.checked; persist(); refreshAll(false); });
  }
  function syncBaselineInputs() {
    el.baseRefWrap.hidden = state.baseline.mode !== 'ref'; el.baseConstWrap.hidden = state.baseline.mode !== 'const';
    el.baseFrom.value = String(state.baseline.refFrom + 1); el.baseTo.value = String(state.baseline.refTo + 1);
    el.baseFrom.max = el.baseTo.max = String(NSEG);
    el.baseConst.value = String(Math.round(state.baseline.value * 1e9));
    el.baseFit.textContent = refLine ? `→ ${fmtAstr(refLine.a)} ${refLine.b < 0 ? '−' : '+'} ${fmtAstr(Math.abs(refLine.b))}/V · V` : (state.baseline.mode === 'ref' ? '→ too few samples' : '');
  }
  segControl('baseline', () => state.baseline.mode, v => { state.baseline.mode = v; persist(); refreshAll(true); });
  const clampIdx = v => Math.max(0, Math.min(LAST, (v | 0) - 1));
  el.baseFrom.addEventListener('change', () => { state.baseline.refFrom = clampIdx(+el.baseFrom.value); if (state.baseline.refTo < state.baseline.refFrom) state.baseline.refTo = state.baseline.refFrom; persist(); refreshAll(true); });
  el.baseTo.addEventListener('change', () => { state.baseline.refTo = clampIdx(+el.baseTo.value); if (state.baseline.refFrom > state.baseline.refTo) state.baseline.refFrom = state.baseline.refTo; persist(); refreshAll(true); });
  el.baseConst.addEventListener('change', () => { state.baseline.value = (+el.baseConst.value || 0) * 1e-9; persist(); refreshAll(true); });

  /* fit parameter sliders */
  function bindParam(input, valueEl, key, fmt) {
    input.value = String(state.params[key]);
    const show = () => { valueEl.textContent = fmt(state.params[key]); };
    show();
    let raf = 0;
    input.addEventListener('input', () => {
      state.params[key] = +input.value; show();
      if (raf) cancelAnimationFrame(raf);
      raf = requestAnimationFrame(() => { raf = 0; persist(); refreshAll(true); });
    });
  }
  bindParam(el.pWindow, el.vWindow, 'window', v => `${v} pts · ${(v * 20 / NG).toFixed(2)} V`);
  bindParam(el.pFion, el.vFion, 'fIon', v => `${Math.round(v * 100)} % · to ${fmtV(-10 + v * 20, 1)}`);
  bindParam(el.pTelo, el.vTelo, 'teLo', v => `${Math.round(v * 100)} %`);
  bindParam(el.pTehi, el.vTehi, 'teHi', v => `${Math.round(v * 100)} %`);
  bindParam(el.pVlo, el.vVlo, 'teVlo', v => fmtV(v, 1));
  bindParam(el.pVhi, el.vVhi, 'teVhi', v => fmtV(v, 1));
  segControl('tewin', () => state.params.teMode, v => { state.params.teMode = v; el.tewinAuto.hidden = v !== 'auto'; el.tewinManual.hidden = v !== 'manual'; persist(); refreshAll(true); });
  el.tewinAuto.hidden = state.params.teMode !== 'auto'; el.tewinManual.hidden = state.params.teMode !== 'manual';
  el.pArea.value = String(+(state.params.area * 1e6).toPrecision(4));
  el.pArea.addEventListener('change', () => { const v = +el.pArea.value; if (v > 0) { state.params.area = v * 1e-6; persist(); refreshAll(true); } });
  (function initGas() {
    const known = Array.from(el.pGas.options).some(o => +o.value === +state.params.mass);
    el.pGas.value = known ? Array.from(el.pGas.options).find(o => +o.value === +state.params.mass).value : 'custom'; el.pMassWrap.hidden = known; el.pMass.value = String(state.params.mass);
  })();
  el.pGas.addEventListener('change', () => { if (el.pGas.value === 'custom') { el.pMassWrap.hidden = false; } else { el.pMassWrap.hidden = true; state.params.mass = +el.pGas.value; persist(); refreshAll(true); } });
  el.pMass.addEventListener('change', () => { const v = +el.pMass.value; if (v > 0) { state.params.mass = v; persist(); refreshAll(true); } });

  /* ======================= boot ======================= */
  recomputeAll();
  syncBaselineInputs();
  buildMinis();
  layout();
  syncTable();
  setSweep(state.sweep, false);
  let resizeRaf = 0;
  const onResize = () => { if (resizeRaf) return; resizeRaf = requestAnimationFrame(() => { resizeRaf = 0; stripCache = null; ghostKey = ''; iv.cur = null; layout(); setSweep(state.sweep, false); }); };
  if (window.ResizeObserver) new ResizeObserver(onResize).observe(el.ivChart); else window.addEventListener('resize', onResize);
  const onTheme = () => { ghostKey = ''; iv.cur = null; layout(); setSweep(state.sweep, false); };
  if (window.matchMedia) matchMedia('(prefers-color-scheme: dark)').addEventListener('change', onTheme);
  new MutationObserver(onTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
})();
