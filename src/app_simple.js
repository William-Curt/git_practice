/* Langmuir Probe Basics — the simple learner view. Depends on d3 (global) and LP (analysis.js). */
(function () {
  'use strict';
  const $ = s => document.querySelector(s);
  const reduceMotion = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
  const NG = 1457, DAC_STEP = 45;
  const FULL_URL = window.FULL_BENCH_URL || '';

  /* ---------- data ---------- */
  function b64ToBytes(b64) { const bin = atob(b64), u8 = new Uint8Array(bin.length); for (let i = 0; i < u8.length; i++) u8[i] = bin.charCodeAt(i); return u8; }
  function decode(raw) {
    const enc = raw.meta.encodings || {}, cols = {};
    for (const name of raw.meta.columns) {
      const bytes = b64ToBytes(raw.cols[name]), e = enc[name] || {};
      const dtype = (e.dtype || 'float32').toLowerCase();
      let arr = dtype.indexOf('16') >= 0 ? new Uint16Array(bytes.buffer, 0, bytes.length >> 1) : new Float32Array(bytes.buffer, 0, bytes.length >> 2);
      if (e.scale !== undefined || e.offset !== undefined) { const f = new Float32Array(arr.length); for (let i = 0; i < f.length; i++) f[i] = arr[i] * (e.scale === undefined ? 1 : e.scale) + (e.offset || 0); arr = f; }
      cols[name] = arr;
    }
    return { meta: raw.meta, bounds: raw.bounds, cols, nSeg: raw.bounds.length - 1 };
  }
  const DATA = decode(window.SWEEP_DATA), NSEG = DATA.nSeg, LAST = NSEG - 1;
  const P = Object.assign({}, LP.DEFAULTS, { baseline: { mode: 'ref', refFrom: Math.max(0, NSEG - 7), refTo: Math.max(0, NSEG - 2), value: 1.48e-7 } });
  const { results } = LP.analyzeAll(DATA, P);
  const T0 = new Date(DATA.meta.t0_iso + (/Z|[+-]\d\d:\d\d$/.test(DATA.meta.t0_iso) ? '' : 'Z')).getTime();
  const clockAt = s => new Date(T0 + s * 1000).toISOString().slice(11, 19);
  const fmtT = s => (s < 60 ? s.toFixed(0) + ' s' : Math.floor(s / 60) + ' min ' + (s % 60).toFixed(0).padStart(2, '0') + ' s');
  const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
  const sup = s => String(s).split('').map(c => SUP[c] || c).join('');
  const fmtSci = n => { if (!(n > 0)) return '—'; const ex = Math.floor(Math.log10(n)); return (n / Math.pow(10, ex)).toFixed(1) + '×10' + sup(ex); };
  const fmtA = a => { if (!Number.isFinite(a)) return '—'; const s = Math.abs(a) >= 1e-6 ? [1e6, 'µA'] : [1e9, 'nA']; const v = a * s[0]; return (v < 0 ? '−' : '') + Math.abs(v).toFixed(Math.abs(v) >= 100 ? 0 : Math.abs(v) >= 10 ? 1 : 2) + ' ' + s[1]; };
  const fmtV = (v, d) => !Number.isFinite(v) ? '—' : (v < 0 ? '−' : '') + Math.abs(v).toFixed(d === undefined ? 1 : d) + ' V';

  $('#when').textContent = new Date(T0).toISOString().slice(0, 16).replace('T', ' ') + ' UTC';
  $('#nsweeps').textContent = String(NSEG);
  if (FULL_URL) $('#full-link').href = FULL_URL;

  /* grid arrays on the common DAC step, ascending V */
  const gridCache = new Array(NSEG);
  function grid(k) {
    if (gridCache[k]) return gridCache[k];
    const r = results[k], V = new Float64Array(NG).fill(NaN), R = new Float64Array(NG).fill(NaN), S = new Float64Array(NG).fill(NaN);
    for (let i = 0; i < r.n; i++) { const g = NG - 1 - Math.round(DATA.cols.DAC[r.rows[i]] / DAC_STEP); if (g < 0 || g >= NG) continue; V[g] = r.V[i]; R[g] = r.I[i]; if (r.Is) S[g] = r.Is[i]; }
    return (gridCache[k] = { V, R, S });
  }
  const analysed = r => !!r.Ie && !r.flags.includes('no_plasma_signal');

  /* ---------- state ---------- */
  const state = { sweep: 1, playing: false, yscale: 'lin' };
  try { const s = JSON.parse(localStorage.getItem('lpb-state') || 'null'); if (s) { state.sweep = Math.min(LAST, Math.max(0, s.sweep | 0)); if (s.yscale) state.yscale = s.yscale; } } catch (e) { /* ignore */ }
  const persist = () => { try { localStorage.setItem('lpb-state', JSON.stringify({ sweep: state.sweep, yscale: state.yscale })); } catch (e) { /* ignore */ } };

  /* ---------- chart ---------- */
  const svg = d3.select('#svg'), chartEl = $('#chart'), tip = $('#tip');
  const m = { l: 66, r: 20, t: 46, b: 40 };
  const gRegions = svg.append('g'), gGrid = svg.append('g'), gAxes = svg.append('g'), gData = svg.append('g'), gFit = svg.append('g'), gMark = svg.append('g'), gHover = svg.append('g');
  const pathMain = gData.append('path').attr('class', 'line-main'), pathTe = gFit.append('path').attr('class', 'line-te');
  const cross = gHover.append('line').attr('class', 'crosshair'), hdot = gHover.append('circle').attr('class', 'hover-dot').attr('r', 5);
  const hit = svg.append('rect').attr('class', 'hit');
  let W = 800, H = 420, cur = null, tween = null, scales = null;

  function layout() {
    W = Math.max(320, chartEl.clientWidth || 800); H = W < 560 ? 340 : 420; m.l = (W < 560 ? 54 : 66) + (state.yscale === 'log' ? 14 : 0);
    svg.attr('viewBox', `0 0 ${W} ${H}`).attr('width', W).attr('height', H);
    hit.attr('x', m.l).attr('y', m.t).attr('width', W - m.l - m.r).attr('height', H - m.t - m.b);
  }
  function frameFor(k) {
    const g = grid(k), log = state.yscale === 'log';
    let lo = Infinity, hi = -Infinity;
    for (const arr of [g.S, g.R]) for (let i = 0; i < NG; i++) { const v = log ? Math.abs(arr[i]) : arr[i]; if (!Number.isFinite(v) || (log && v <= 0)) continue; if (v < lo) lo = v; if (v > hi) hi = v; }
    if (!(lo < hi)) { lo = log ? 1e-9 : 0; hi = log ? 1e-6 : 1e-7; }
    const dom = log ? [Math.max(lo, 1e-10) / 1.5, hi * 1.6] : [Math.min(lo, 0) - 0.06 * (hi - Math.min(lo, 0)), hi + 0.10 * (hi - Math.min(lo, 0))];
    return { k, V: g.V, R: g.R, S: g.S, dom };
  }
  const lerpArr = (a, b, t) => { const o = new Float64Array(b.length); for (let i = 0; i < b.length; i++) { const x = a ? a[i] : NaN; o[i] = Number.isFinite(x) && Number.isFinite(b[i]) ? x + (b[i] - x) * t : b[i]; } return o; };
  const ease = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);
  function show(k, animate) {
    const to = frameFor(k), from = cur, log = state.yscale === 'log';
    const dur = (!animate || reduceMotion || !from) ? 0 : (state.playing ? 380 : 340);
    if (tween) cancelAnimationFrame(tween);
    if (!dur) { cur = to; draw(to, false); return; }
    const t0 = performance.now();
    const step = now => {
      const t = Math.min(1, (now - t0) / dur), e = ease(t);
      const dom = log ? [Math.exp(Math.log(from.dom[0]) + (Math.log(to.dom[0]) - Math.log(from.dom[0])) * e), Math.exp(Math.log(from.dom[1]) + (Math.log(to.dom[1]) - Math.log(from.dom[1])) * e)] : [from.dom[0] + (to.dom[0] - from.dom[0]) * e, from.dom[1] + (to.dom[1] - from.dom[1]) * e];
      draw({ k: to.k, V: lerpArr(from.V, to.V, e), R: t < 1 ? lerpArr(from.R, to.R, e) : to.R, S: lerpArr(from.S, to.S, e), dom }, t < 1);
      if (t < 1) tween = requestAnimationFrame(step); else { cur = to; tween = null; }
    };
    tween = requestAnimationFrame(step);
  }
  const unitFor = v => (v >= 1e-6 ? { s: 1e6, u: 'µA' } : { s: 1e9, u: 'nA' });
  const fmtSI = d => { const u = unitFor(d), v = d * u.s; return (v >= 1 ? v.toFixed(0) : v.toPrecision(1)) + ' ' + u.u; };
  function draw(f, inFlight) {
    const r = results[f.k], log = state.yscale === 'log';
    const x = d3.scaleLinear().domain([-10.3, 10.3]).range([m.l, W - m.r]);
    const y = (log ? d3.scaleLog() : d3.scaleLinear()).domain(f.dom).range([H - m.b, m.t]);
    const yv = v => log ? (Math.abs(v) > 0 ? y(Math.abs(v)) : NaN) : y(v);
    const unit = unitFor(Math.max(Math.abs(f.dom[0]), Math.abs(f.dom[1])));
    const xt = x.ticks(W < 560 ? 5 : 11);
    let yt = y.ticks(6);
    if (log) { const pow10 = d => Math.abs(Math.log10(d) - Math.round(Math.log10(d))) < 1e-9; yt = y.ticks().filter(pow10); }
    const step = yt.length > 1 ? Math.abs(yt[1] - yt[0]) * unit.s : 1, dec = Math.max(0, Math.min(3, Math.ceil(-Math.log10(step) - 1e-9)));
    gGrid.selectAll('line.gx').data(xt).join('line').attr('class', 'grid-line gx').attr('x1', d => x(d)).attr('x2', d => x(d)).attr('y1', m.t).attr('y2', H - m.b);
    gGrid.selectAll('line.gy').data(yt).join('line').attr('class', 'grid-line gy').attr('x1', m.l).attr('x2', W - m.r).attr('y1', d => y(d)).attr('y2', d => y(d));
    gGrid.selectAll('line.zero').data(!log && f.dom[0] < 0 ? [0] : []).join('line').attr('class', 'zero-line zero').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(0)).attr('y2', y(0));
    gAxes.selectAll('line.ax').data([0]).join('line').attr('class', 'axis-line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', H - m.b).attr('y2', H - m.b);
    gAxes.selectAll('text.tx').data(xt).join('text').attr('class', 'tick-text tx').attr('x', d => x(d)).attr('y', H - m.b + 15).attr('text-anchor', 'middle').text(d => (d < 0 ? '−' : d > 0 ? '+' : '') + Math.abs(d) + ' V');
    gAxes.selectAll('text.ty').data(yt).join('text').attr('class', 'tick-text ty').attr('x', m.l - 8).attr('y', d => y(d) + 3.5).attr('text-anchor', 'end').text(d => log ? fmtSI(d) : (d < 0 ? '−' : '') + (Math.abs(d) * unit.s).toFixed(dec));
    gAxes.selectAll('text.yl').data([0]).join('text').attr('class', 'axis-title').attr('transform', `translate(14,${(H - m.b + m.t) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text(log ? '|current| (log)' : `current (${unit.u})`);
    gAxes.selectAll('text.xl').data([0]).join('text').attr('class', 'axis-title').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').text('probe bias (V)');

    /* regions: ion side / transition / electron side */
    const ok = analysed(r);
    const vf = ok && r.Vf !== null ? r.Vf : null, tw = ok ? r.Te_window : [NaN, NaN];
    const regs = [];
    if (ok && vf !== null) {
      regs.push({ id: 'ion', a: -10.3, b: vf, k: 'var(--ion)', l: 'ion side', s: 'only ions arrive: small negative current' });
      const teHi = Number.isFinite(tw[1]) ? tw[1] : Math.min(vf + 5, 10.3);
      regs.push({ id: 'trans', a: vf, b: teHi, k: 'var(--trans)', l: 'transition', s: 'electron current grows exponentially' });
      regs.push({ id: 'elec', a: teHi, b: 10.3, k: 'var(--elec)', l: 'electron side', s: r.Vp === null ? 'still rising: knee is beyond +10 V' : 'knee (plasma potential) at ' + fmtV(r.Vp) });
    }
    const rg = gRegions.selectAll('g.reg').data(regs, d => d.id).join(enter => { const g = enter.append('g').attr('class', 'reg'); g.append('rect').attr('class', 'region'); g.append('text').attr('class', 'region-label'); g.append('text').attr('class', 'region-sub'); return g; });
    rg.select('rect').attr('x', d => x(d.a)).attr('y', m.t).attr('width', d => Math.max(0, x(d.b) - x(d.a))).attr('height', H - m.t - m.b).attr('fill', d => d.k);
    const narrow = W < 560;
    rg.select('.region-label').attr('x', d => (x(d.a) + x(d.b)) / 2).attr('y', m.t - (narrow ? 22 : 24)).attr('text-anchor', 'middle').text(d => d.l);
    rg.select('.region-sub').attr('x', d => (x(d.a) + x(d.b)) / 2).attr('y', m.t - (narrow ? 9 : 10)).attr('text-anchor', 'middle').text(d => narrow ? '' : d.s);

    /* data */
    const line = d3.line().defined(d => Number.isFinite(d[0]) && Number.isFinite(d[1])).x(d => x(d[0])).y(d => d[1]);
    const pts = (V, Y) => { const o = new Array(NG); for (let i = 0; i < NG; i++) o[i] = [V[i], Number.isFinite(Y[i]) ? yv(Y[i]) : NaN]; return o; };
    pathMain.attr('d', line(pts(f.V, f.S)));
    const sample = []; for (let i = 0; i < NG; i += 2) { const yy = Number.isFinite(f.R[i]) ? yv(f.R[i]) : NaN; if (Number.isFinite(f.V[i]) && Number.isFinite(yy)) sample.push([x(f.V[i]), yy]); }
    gData.selectAll('circle.pt').data(sample).join('circle').attr('class', 'pt').attr('r', 1.8).attr('cx', d => d[0]).attr('cy', d => d[1]).style('opacity', inFlight ? 0.18 : null);

    /* Te fit + markers */
    if (ok && r.Te !== null && Number.isFinite(tw[0])) {
      const lo = tw[0] - 0.6, hi = Math.min(tw[1] + 1.0, 10);
      const vs = d3.range(lo, hi + 0.001, (hi - lo) / 80);
      pathTe.attr('d', line(vs.map(v => [v, yv(Math.exp(r.te_alpha + r.te_beta * v) + r.a_i + r.b_i * v)])));
    } else pathTe.attr('d', null);
    const marks = [];
    if (ok && vf !== null) marks.push({ id: 'vf', v: vf, l: 'Vf', s: fmtV(vf, 2), y: yv(0) });
    if (ok && r.Vp !== null) marks.push({ id: 'vp', v: r.Vp, l: 'Vp', s: fmtV(r.Vp, 2), y: m.t + 14 });
    const mk = gMark.selectAll('g.mk').data(marks, d => d.id).join(enter => { const g = enter.append('g').attr('class', 'mk'); g.append('line').attr('class', 'marker-line'); g.append('circle').attr('r', 5).attr('fill', 'var(--surface)').attr('stroke', 'var(--ink)').attr('stroke-width', 2); g.append('text').attr('class', 'marker-text'); g.append('text').attr('class', 'marker-sub'); return g; });
    mk.attr('transform', d => `translate(${x(d.v)},0)`);
    mk.select('line').attr('y1', m.t).attr('y2', H - m.b);
    mk.select('circle').attr('cy', d => Number.isFinite(d.y) ? d.y : m.t + 14).style('display', d => (d.id === 'vf' && !log) ? null : 'none');
    mk.select('.marker-text').attr('x', d => x(d.v) > W - 120 ? -8 : 8).attr('text-anchor', d => x(d.v) > W - 120 ? 'end' : 'start').attr('y', H - m.b - 26).text(d => d.l);
    mk.select('.marker-sub').attr('x', d => x(d.v) > W - 120 ? -8 : 8).attr('text-anchor', d => x(d.v) > W - 120 ? 'end' : 'start').attr('y', H - m.b - 12).text(d => d.s);
    gMark.selectAll('text.note').data(ok ? [] : [r.flags.includes('partial_sweep') ? 'partial sweep: the log starts or ends here, so it is not analysed' : r.flags.includes('no_plasma_signal') ? 'no plasma signal: this sweep was recorded with the plasma off' : 'not analysed']).join('text').attr('class', 'note').attr('x', m.l + 12).attr('y', m.t + 22).text(d => d);
    scales = { x, y, yv };
    if (hoverG >= 0) hoverAt(hoverG);
  }
  let hoverG = -1;
  function hoverAt(g) {
    if (!scales || !cur) return;
    const gr = grid(state.sweep);
    if (g < 0 || !Number.isFinite(gr.V[g])) { cross.style('opacity', 0); hdot.style('opacity', 0); tip.classList.remove('on'); return; }
    const xx = scales.x(gr.V[g]), yy = scales.yv(gr.S[g]);
    cross.attr('x1', xx).attr('x2', xx).attr('y1', m.t).attr('y2', H - m.b).style('opacity', 1);
    hdot.attr('cx', xx).attr('cy', yy).style('opacity', Number.isFinite(yy) ? 1 : 0);
    tip.textContent = ''; const a = document.createElement('div'); a.textContent = 'bias ' + fmtV(gr.V[g], 2); const b = document.createElement('div'); const bb = document.createElement('b'); bb.textContent = fmtA(gr.S[g]); b.appendChild(document.createTextNode('current ')); b.appendChild(bb); tip.appendChild(a); tip.appendChild(b);
    tip.classList.add('on'); const tw = tip.offsetWidth || 120; tip.style.left = (xx + 14 + tw > W ? xx - 14 - tw : xx + 14) + 'px'; tip.style.top = Math.max(0, (Number.isFinite(yy) ? yy : m.t + 30) - 24) + 'px';
  }
  hit.on('pointermove', ev => { const [mx] = d3.pointer(ev, svg.node()); const v = scales.x.invert(mx), V = grid(state.sweep).V; let best = -1, bd = Infinity; for (let g = 0; g < NG; g++) { if (!Number.isFinite(V[g])) continue; const d = Math.abs(V[g] - v); if (d < bd) { bd = d; best = g; } } hoverG = best; hoverAt(best); }).on('pointerleave', () => { hoverG = -1; hoverAt(-1); });

  /* ---------- the three numbers ---------- */
  function renderNumbers() {
    const r = results[state.sweep], ok = analysed(r), box = $('#numbers'); box.innerHTML = '';
    const card = (label, value, unit, meaning, dim, badge) => {
      const d = document.createElement('div'); d.className = 'card num';
      const l = document.createElement('div'); l.className = 'l'; l.textContent = label;
      const v = document.createElement('div'); v.className = 'v' + (dim ? ' dim' : ''); v.textContent = value; if (unit) { const s = document.createElement('small'); s.textContent = unit; v.appendChild(s); }
      const mm = document.createElement('div'); mm.className = 'm'; mm.innerHTML = meaning;
      d.appendChild(l); d.appendChild(v); d.appendChild(mm);
      if (badge) { const b = document.createElement('span'); b.className = 'badge' + (badge[1] ? ' crit' : ''); b.textContent = badge[0]; d.appendChild(b); }
      box.appendChild(d);
    };
    const why = r.flags.includes('partial_sweep') ? 'This sweep is cut off at the start or end of the log.' : r.flags.includes('no_plasma_signal') ? 'No plasma signal: the plasma was off for this sweep.' : 'This sweep could not be analysed.';
    if (!ok) { card('Floating potential', 'not analysed', '', why, true); card('Electron temperature', 'not analysed', '', why, true); card('Electron density', 'not analysed', '', why, true); return; }
    card('Floating potential  Vf', r.Vf !== null ? fmtV(r.Vf, 2).replace(' V', '') : '—', r.Vf !== null ? 'V' : '',
      r.Vf !== null ? `The bias at which the probe collects <b>zero net current</b>: electrons and ions arrive at the same rate. An insulated object in this plasma would sit at this voltage${r.dVf !== null ? ` (uncertainty about ±${r.dVf.toFixed(2)} V from noise)` : ''}.` : 'The current never crossed zero inside the sweep.', r.Vf === null);
    const kelvin = r.Te !== null ? Math.round(r.Te * 11604.5 / 1000) : null;
    card('Electron temperature  Te', r.Te !== null ? r.Te.toFixed(2) : '—', r.Te !== null ? 'eV' : '',
      r.Te !== null ? `About <b>${kelvin.toLocaleString()},000 K</b>. From the slope of the exponential part of the curve (the dashed fit, from ${fmtV(r.Te_window[0])} to ${fmtV(r.Te_window[1])}). Fit quality R² = ${r.Te_r2.toFixed(2)}${r.flags.includes('te_nonexponential') ? '; the slope is not constant, so this is an average over that window' : ''}.` : 'Too few points for a fit in this sweep.', r.Te === null,
      r.flags.includes('top_of_sweep_unreliable') ? ['recording glitch near +10 V', true] : null);
    const lower = r.Vp === null;
    card('Electron density  ne' + (lower ? '  (at least)' : ''), r.ne ? fmtSci(r.ne) : '—', r.ne ? 'm⁻³' : '',
      r.ne ? `That is about <b>${fmtSci(r.ne / 1e6)} per cm³</b>, assuming a 10 mm² probe tip. ${lower ? 'A <b>lower bound</b>: the curve is still rising at +10 V, so the true electron saturation current (and density) is higher.' : 'From the electron saturation current at the knee.'}` : (r.flags.includes('top_of_sweep_unreliable') ? 'The top of this sweep was lost to a recording glitch, so no density estimate.' : 'Needs both a temperature and a current at +10 V.'), !r.ne,
      r.flags.includes('top_of_sweep_unreliable') ? ['recording glitch near +10 V', true] : null);
  }

  /* ---------- run overview (mini chart) ---------- */
  const msvg = d3.select('#mini-svg'), miniEl = $('#mini'), mtip = $('#mini-tip');
  const gm = msvg.append('g');
  function drawMini() {
    const w = Math.max(280, miniEl.clientWidth || 800), h = 140, ml = 56, mr = 12, mt = 10, mb = 22;
    msvg.attr('viewBox', `0 0 ${w} ${h}`).attr('width', w).attr('height', h);
    const T = DATA.cols.ts[DATA.cols.ts.length - 1] / 60;
    const x = d3.scaleLinear().domain([0, T]).range([ml, w - mr]);
    const pts = results.map(r => ({ k: r.k, t: (r.tStart + r.duration / 2) / 60, v: r.ne })).filter(d => d.v);
    const y = d3.scaleLinear().domain([0, (d3.max(pts, d => d.v) || 1) * 1.1]).nice(4).range([h - mb, mt]);
    gm.selectAll('*').remove();
    gm.selectAll('rect.band').data(results.filter(r => !analysed(r) || r.flags.includes('top_of_sweep_unreliable'))).join('rect').attr('class', 'band').attr('x', d => x(d.tStart / 60)).attr('y', mt).attr('width', d => Math.max(1, x((d.tStart + d.duration) / 60) - x(d.tStart / 60))).attr('height', h - mt - mb);
    const yt = y.ticks(3);
    gm.selectAll('line.gy').data(yt).join('line').attr('class', 'grid-line').attr('x1', ml).attr('x2', w - mr).attr('y1', d => y(d)).attr('y2', d => y(d));
    gm.selectAll('text.ty').data(yt).join('text').attr('class', 'tick-text').attr('x', ml - 6).attr('y', d => y(d) + 3.5).attr('text-anchor', 'end').text(d => d === 0 ? '0' : d3.format('.1~e')(d).replace('e+', 'e'));
    gm.append('text').attr('class', 'axis-title').attr('transform', `translate(11,${(h - mb + mt) / 2}) rotate(-90)`).attr('text-anchor', 'middle').text('nₑ (m⁻³)');
    gm.append('line').attr('class', 'axis-line').attr('x1', ml).attr('x2', w - mr).attr('y1', h - mb).attr('y2', h - mb);
    const xt = x.ticks(w < 500 ? 3 : 5).filter(d => d > 0);
    gm.selectAll('text.tx').data(xt).join('text').attr('class', 'tick-text').attr('x', d => x(d)).attr('y', h - mb + 14).attr('text-anchor', 'middle').text(d => d + ' min');
    const runs = []; let run = []; for (const d of pts) { if (run.length && d.k !== run[run.length - 1].k + 1) { runs.push(run); run = []; } run.push(d); } if (run.length) runs.push(run);
    const line = d3.line().x(d => x(d.t)).y(d => y(d.v)), area = d3.area().x(d => x(d.t)).y0(h - mb).y1(d => y(d.v));
    gm.selectAll('path.area').data(runs).join('path').attr('class', 'area').attr('d', d => area(d));
    gm.selectAll('path.line').data(runs).join('path').attr('class', 'line').attr('d', d => line(d));
    gm.selectAll('circle.pt').data(pts, d => d.k).join('circle').attr('class', d => 'pt' + (d.k === state.sweep ? ' cur' : '')).attr('cx', d => x(d.t)).attr('cy', d => y(d.v)).attr('r', d => d.k === state.sweep ? 6 : 3);
    const rc = results[state.sweep]; gm.append('line').attr('class', 'cursor').attr('x1', x((rc.tStart + rc.duration / 2) / 60)).attr('x2', x((rc.tStart + rc.duration / 2) / 60)).attr('y1', mt).attr('y2', h - mb).style('opacity', 0.6);
    $('#ov-cur').textContent = rc.ne ? `sweep ${state.sweep + 1}: ${fmtSci(rc.ne)} m⁻³` : `sweep ${state.sweep + 1}: not analysed`;
    const nearest = mx => { const t = x.invert(mx); let best = 0, bd = Infinity; for (const r of results) { const d = Math.abs((r.tStart + r.duration / 2) / 60 - t); if (d < bd) { bd = d; best = r.k; } } return best; };
    gm.append('rect').attr('class', 'hit').attr('x', ml).attr('y', 0).attr('width', w - ml - mr).attr('height', h)
      .on('pointermove', ev => { const [mx] = d3.pointer(ev); const k = nearest(mx), r = results[k]; mtip.textContent = `sweep ${k + 1} · ${fmtT(r.tStart)} · ${r.ne ? fmtSci(r.ne) + ' m⁻³' : 'not analysed'}`; mtip.classList.add('on'); const tw = mtip.offsetWidth || 150; mtip.style.left = Math.max(0, Math.min(mx + 12, w - tw)) + 'px'; mtip.style.top = '4px'; })
      .on('pointerleave', () => mtip.classList.remove('on')).on('click', ev => { const [mx] = d3.pointer(ev); setSweep(nearest(mx), true); });
  }

  /* ---------- selection & playback ---------- */
  const elSweep = $('#sweep'), elPlay = $('#play');
  elSweep.max = String(LAST);
  function setSweep(k, animate) {
    k = Math.max(0, Math.min(LAST, k | 0)); state.sweep = k; elSweep.value = String(k); elSweep.style.setProperty('--pct', (k / LAST * 100) + '%');
    const r = results[k];
    $('#sweep-label').textContent = `Sweep ${k + 1} of ${NSEG} · ${r.dir === 'up' ? '−10 → +10 V' : '+10 → −10 V'}`;
    $('#sweep-time').textContent = `${fmtT(r.tStart)} into the run · ${clockAt(r.tStart)} UTC`;
    show(k, animate !== false); renderNumbers(); drawMini(); persist();
  }
  elSweep.addEventListener('input', () => setSweep(+elSweep.value, true));
  let raf = null, last = 0;
  function tick(now) { if (!state.playing) return; if (now - last >= 1000) { last = now; setSweep(state.sweep >= LAST ? 0 : state.sweep + 1, true); } raf = requestAnimationFrame(tick); }
  function setPlaying(on) { state.playing = on; elPlay.setAttribute('aria-pressed', String(on)); elPlay.setAttribute('aria-label', on ? 'Pause' : 'Play through the sweeps'); if (on) { last = performance.now() - 1e9; raf = requestAnimationFrame(tick); } else if (raf) cancelAnimationFrame(raf); }
  elPlay.addEventListener('click', () => setPlaying(!state.playing));
  const ys = $('#yscale');
  const syncY = () => ys.querySelectorAll('button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.v === state.yscale)));
  ys.addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; state.yscale = b.dataset.v; syncY(); persist(); cur = null; layout(); setSweep(state.sweep, false); });
  syncY();
  document.addEventListener('keydown', ev => {
    if (ev.target && /^(INPUT|SELECT|TEXTAREA|BUTTON)$/.test(ev.target.tagName) && ev.target.type !== 'range') return;
    if (ev.key === 'ArrowRight') { ev.preventDefault(); setSweep(state.sweep + 1, true); }
    else if (ev.key === 'ArrowLeft') { ev.preventDefault(); setSweep(state.sweep - 1, true); }
    else if (ev.key === ' ' && ev.target.tagName !== 'INPUT') { ev.preventDefault(); setPlaying(!state.playing); }
  });

  /* ---------- boot ---------- */
  layout(); setSweep(state.sweep, false);
  let rr = 0; const onResize = () => { if (rr) return; rr = requestAnimationFrame(() => { rr = 0; cur = null; layout(); setSweep(state.sweep, false); }); };
  if (window.ResizeObserver) new ResizeObserver(onResize).observe(chartEl); else window.addEventListener('resize', onResize);
})();
