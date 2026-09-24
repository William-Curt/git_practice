#!/usr/bin/env python3
"""
Langmuir sweep analysis -- reference implementation of ANALYSIS_SPEC_v2.md
including the "v2.1 amendments" section (A: over-range sweeps get
'top_of_sweep_unreliable' and no Vp/Ies/ne/lambdaD/hysteresis; B: hysteresis
h from a coarse 0.1 V / fine 0.02 V search of Is_up(V) vs Is_down(V + h);
C (v2.2): Vp is accepted only if additionally Is[hi] >= Is[ip] - 3*sigma_I;
D (v2.3): partner(k) = k+1 for even k, k-1 for odd k, partner > nSeg-1 -> nSeg-2;
E (v2.3): when Vp is not accepted iRef = first maximum of Ie_s over [lo, hi];
F (v2.3): gate derivative threshold 5*sigma_b; G (v2.3): the gate removes the
step-5 flags 'no_zero_crossing' / 'vf_multiple_crossings' when it fires;
H (v2.4): Vp acceptance is a four-condition flattening test (see A17), with
condition (b) in its baseline-invariant v2.4b form (amendment K);
I (v2.4): per-segment gate_reason 'span' | 'slope' | null;
J (v2.4): the pair step copies a partner's value only for symmetric pairs).

Only numpy + pandas are used (no scipy).  Every numerical step is plain
arithmetic (window means, ordinary least squares, argmax, clamped linear
interpolation) so the algorithm can be ported line-for-line to JavaScript.
To make a bit-exact port possible, EVERY sum in this file is accumulated
sequentially left-to-right in index order (no pairwise / Kahan summation):
the window primitives are vectorised across windows but accumulate column by
column, which is arithmetically identical to a scalar `for` loop.

CLI:
  python langmuir_analysis.py <csv> [--channel I2|I1]
        [--baseline none|const:<A>|ref:<from>,<to>]      (default ref:74,79)
        [--f-ion 0.2] [--window 31] [--te-lo 0.02] [--te-hi 0.30] [--te-manual Vlo,Vhi]
        [--area 1e-5] [--mass 39.948] [--out results.json]
        [--dump-segments k1,k2,...] [--csv-out table.csv] [--quantize-float32]

--quantize-float32 rounds V2, I1, I2 and the seconds-since-first-timestamp
through float32 right after loading (t_start/duration then come from the
float32 seconds), so the reference sees exactly the values the applet's
Float32Array columns hold.

Resolved spec ambiguities (all choices are the most literal reading; each is
also listed in the final report):

 A1  Step 0a uses I2 for the screening even when --channel I1.  A NaN in I2
     never satisfies any of the three comparisons, so such a row is not `bad`
     (I2 has no NaN in the reference data anyway).
 A2  n_dropped = (r1-r0) - n counts BOTH bad rows and non-finite channel rows;
     n_glitch counts the bad rows only.
 A3  'baseline_unavailable' (ref pool < 50 points) is put on EVERY segment,
     including those later flagged 'insufficient'.  Ref indices outside
     [0, nSeg-1] are ignored when building the pool.  With mode 'none' or
     'const', or after the fall-back, top-level a_r / b_r are null.
 A4  For an 'insufficient' segment the algorithm stops at step 0b, so Vmin,
     Vmax, sigma_*, I_lo, I_hi are null as well (only k, dir, t_start,
     duration, n, n_dropped, n_glitch, flags are filled).
 A5  dIdV_max is defined in step 6 (dIdV[ip]); for a 'partial_sweep' segment
     steps 3-11 are skipped, so dIdV_max is null.  Is / dIdV / sigma_I /
     sigma_b / I_lo / I_hi are still reported.  In a dump the arrays that are
     only produced by skipped steps (Ie, Ie_s, Iion, Te_local, te_mask,
     ion_mask) are null, not empty arrays.
 A6  Step 5 is evaluated before the step-6 gate.  When 'no_plasma_signal'
     fires, Vf and dVf are reset to null and (amendment G) the step-5 flags
     'no_zero_crossing' / 'vf_multiple_crossings' are removed.  The Vp
     acceptance test is not run in that case, so 'vp_beyond_range' /
     'not_saturated' are NOT added by the gate.
 A7  iVf = index of the point with V closest to Vf: first index on a tie
     (argmin of |V - Vf|).
 A8  In the multiple-crossings branch the local OLS uses every point with
     |V - Vc_m| <= 0.5; if that fit has no x-variance (b is NaN) it is
     treated as "slope not > 0" -> Vf = Vc_m.
 A9  Te_npts = |teMask| is always reported (also when < 8); Te_window is
     [min V, max V] over teMask when teMask is non-empty, else null.  In
     manual mode Te_window is ALSO [min V, max V] over teMask (the spec gives
     one definition for both modes), not the requested [Vlo, Vhi].
 A10 te_alpha / te_beta are the OLS coefficients and are reported whenever
     the fit was performed (|teMask| >= 8), including when beta <= 0.
     'te_poor_fit' is evaluated from Te_r2 independently of the sign of beta.
     The split diagnostic (Te_lo / Te_hi) is only run when the main fit was
     performed (|teMask| >= 8).
 A11 Auto Te walk: j starts at iRef.  Phase 1 skips while Ie_s[j] > teHi*Ie_ref;
     phase 2 collects while (Ie_s[j] >= teLo*Ie_ref and Ie_s[j] > 3*sigma_I
     and (Vf null or V[j] > Vf)); the first failure ends the walk (j reaching
     -1 also ends it).  teMask is stored ascending by V.
 A12 Te_local is part of step 11 and is therefore null for a segment that is
     partial / insufficient / no_plasma_signal / ion_fit_insufficient.
 A13 Hysteresis (v2.1 amendment B, pairing per v2.3 amendment D):
     partner(k) = k+1 for even k, k-1 for odd k; if partner > nSeg-1 then
     partner = nSeg-2.  The value is computed per unordered pair
     {k, partner(k)} and stored on k; partner() is symmetric for every pair
     (2m, 2m+1), which is the "same value on both members" of the spec.  The
     only asymmetric case is an even last segment (nSeg odd): its partner
     nSeg-2 is paired with nSeg-3 by its own rule, so the last segment's
     value never overwrites its partner's.  If the partner is k itself or
     both members have the same direction, hysteresis = null.  A pair is evaluated only when neither
     member has 'overrange_event', 'partial_sweep', 'insufficient' or
     'no_plasma_signal' (an 'ion_fit_insufficient' segment still counts
     because Is exists).  rms(h) is taken over the up-sweep points with
     -5 <= V <= 9 (inclusive) of (Is_up(V) - interp(V_down, Is_down, V+h));
     if there are none -> null.  Grids: coarse h = (i-20)/10, i = 0..40;
     fine h = (c + 2*(j-5)) / 100, j = 0..10, where c = 100*h0 as an exact
     integer -- integer hundredths divided by 100 IS the "rounded to 2
     decimals" value, bit-exactly, in Python and JS alike.  "First minimum"
     = the smallest h among equal minima (strict '<' scan).  |h| >= 1.98 on
     the fine pass -> null and 'hysteresis_unresolved' appended to the flags
     of BOTH members (literal reading; only matters for the asymmetric (0,1)
     pair, which cannot occur here because segment 0 is partial).
 A13b Amendment A: 'top_of_sweep_unreliable' is added right after the
     signal gate for any segment with n_glitch > 0; the Vp test is not run;
     Ies/ne/lambdaD/hysteresis are null.  A glitch segment that fails the
     gate is handled by the gate alone ('no_plasma_signal', no
     'top_of_sweep_unreliable').  Amendment E (later) supersedes A's
     "iRef = hi": whenever Vp is not accepted, including glitch segments,
     iRef = lo + first argmax of Ie_s over [lo, hi].
 A14 interp(xp, fp, x) has numpy.interp semantics: j = largest index with
     xp[j] <= x; x == xp[j] -> fp[j]; otherwise fp[j] + (x - xp[j]) *
     ((fp[j+1] - fp[j]) / (xp[j+1] - xp[j])); x < xp[0] -> fp[0];
     x >= xp[-1] -> fp[-1].
 A15 Direction: 'up' iff V2 of the segment's last row > V2 of its first row.
 A16 Dump extras: `excluded_rows` = absolute 0-based data-row indices of the
     bad (overrange) rows, `excluded_VI` = their [V2, I2] pairs (always I2,
     as requested, regardless of --channel).  NaN-channel rows that are not
     bad are not listed (they have no I value to plot).
 A17 Amendment H (flattening test), non-glitch segments only: Vp = V[ip] iff
     (a) V[ip] <= Vmax - 1.0; (b) Is[hi] >= Is[ip] - 3*sigma_I and
     Is[hi] - Is[ip] <= 0.4*(Is[ip] - Is[lo]) (amendment K, v2.4b; the
     original H form was Is[hi] <= 1.4*Is[ip]); (c) the tail set T = {dIdV[j] : ip < j <= hi, V[j] >= V[ip] + 0.3} is
     non-empty and median(T) <= 0.7*dIdV[ip], median = middle element of the
     sorted values, or (s[m/2-1] + s[m/2]) / 2 for even m; (d) Vf is null or
     V[ip] > Vf + 1.0, where Vf is the step-5 value.  The four booleans, the
     tail median and the tail count are written to dump.idx.vp_test.
 A18 Amendment I: gate_reason is a per-segment output key ('span' if the
     span test failed, else 'slope' if the derivative test failed, else
     null); it is null for segments that never reach step 6.
 A19 Amendment J: a segment copies its partner's cached value (and shares
     the 'hysteresis_unresolved' flag) only when partner(partner(k)) == k;
     an asymmetric partner (only possible for an even last segment) is
     evaluated on its own and rejected if the pair is not eligible; its
     unresolved flag, if any, is put on that segment only.
"""
import argparse
import json
import math

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view

# ---------------------------------------------------------------------------
# Physical constants (spec step 10)
# ---------------------------------------------------------------------------
E_CHARGE = 1.602176634e-19      # C
M_E = 9.1093837015e-31          # kg
AMU = 1.66053906660e-27         # kg
EPS0 = 8.8541878128e-12         # F/m

# Screening thresholds (spec step 0a)
OVERRANGE_NEG = -5e-8           # A
OVERRANGE_JUMP = 1e-6           # A
# Signal gate (step 6): dIdV_max must reach GATE_DERIV_FACTOR * sigma_b (v2.3 amendment F)
GATE_DERIV_FACTOR = 5.0

SCALAR_KEYS = ["k", "dir", "t_start", "duration", "n", "n_dropped", "n_glitch",
               "Vmin", "Vmax", "sigma_I", "sigma_b", "a_i", "b_i",
               "Vf", "dVf", "Vp", "Te", "Te_r2", "Te_window", "Te_npts", "Te_lo", "Te_hi",
               "te_alpha", "te_beta", "Ies", "Iis", "Iis_vf", "ne", "ni", "lambdaD",
               "Vp_expected", "I_lo", "I_hi", "dIdV_max", "gate_reason", "hysteresis", "flags"]

DUMP_KEYS = ["V", "I", "Is", "dIdV", "Ie", "Ie_s", "Iion", "Te_local", "te_mask", "ion_mask",
             "excluded_rows", "excluded_VI", "idx", "crossings"]


# ---------------------------------------------------------------------------
# Plain-arithmetic primitives (port these 1:1 to JS)
# ---------------------------------------------------------------------------
def seq_sum(values):
    """Sequential left-to-right sum: ((v0 + v1) + v2) + ..."""
    s = 0.0
    for v in values:
        s += v
    return s


def seq_sum_rows(M):
    """Sequential left-to-right sum along axis 1 of a 2-D array.
    Row r yields ((M[r,0] + M[r,1]) + M[r,2]) + ... -- bit-identical to seq_sum."""
    acc = M[:, 0].copy()
    for j in range(1, M.shape[1]):
        acc = acc + M[:, j]
    return acc


def median_plain(values):
    """Median of a plain list: middle of the sorted values, or the mean of the
    two middle values for an even count (spec amendment H)."""
    s = sorted(values)
    m = len(s)
    if m % 2 == 1:
        return s[m // 2]
    return (s[m // 2 - 1] + s[m // 2]) / 2.0


def window_bounds(i, n, h):
    """Clipped index window [lo, hi] (inclusive) around i: max(0,i-h) .. min(n-1,i+h)."""
    lo = i - h
    if lo < 0:
        lo = 0
    hi = i + h
    if hi > n - 1:
        hi = n - 1
    return lo, hi


def ols_fit(x, y):
    """Ordinary least squares y = a + b*x on plain sequences.  Returns (a, b, r2).
    b is NaN if x has zero variance; r2 is NaN if y has zero variance.
    All sums sequential in index order."""
    m = len(x)
    xbar = 0.0
    ybar = 0.0
    for j in range(m):
        xbar += x[j]
        ybar += y[j]
    xbar /= m
    ybar /= m
    sxx = 0.0
    sxy = 0.0
    syy = 0.0
    for j in range(m):
        dx = x[j] - xbar
        dy = y[j] - ybar
        sxx += dx * dx
        sxy += dx * dy
        syy += dy * dy
    if sxx == 0.0:
        return float("nan"), float("nan"), float("nan")
    b = sxy / sxx
    a = ybar - b * xbar
    if syy == 0.0:
        r2 = float("nan")
    else:
        ss_res = 0.0
        for j in range(m):
            r = y[j] - (a + b * x[j])
            ss_res += r * r
        r2 = 1.0 - ss_res / syy
    return a, b, r2


def smooth_mean(I, h):
    """Step 1: Is[i] = mean(I[max(0,i-h) .. min(n-1,i+h)]) (sequential sums)."""
    n = len(I)
    w = 2 * h + 1
    Is = np.empty(n)
    if n >= w:
        W = sliding_window_view(I, w)              # row r = I[r : r+w], centre i = r + h
        Is[h:n - h] = seq_sum_rows(W) / w
        ends = list(range(0, h)) + list(range(n - h, n))
    else:
        ends = list(range(n))
    Il = I.tolist()
    for i in ends:
        lo, hi = window_bounds(i, n, h)
        Is[i] = seq_sum(Il[lo:hi + 1]) / (hi - lo + 1)
    return Is


def window_slope(V, I, h):
    """Step 2: dIdV[i] = OLS slope of I vs V over the same clipped window
    (same arithmetic and summation order as ols_fit)."""
    n = len(I)
    w = 2 * h + 1
    d = np.empty(n)
    if n >= w:
        Vw = sliding_window_view(V, w)
        Iw = sliding_window_view(I, w)
        xbar = seq_sum_rows(Vw) / w
        ybar = seq_sum_rows(Iw) / w
        dx = Vw - xbar[:, None]
        dy = Iw - ybar[:, None]
        sxx = seq_sum_rows(dx * dx)
        sxy = seq_sum_rows(dx * dy)
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = np.where(sxx == 0.0, np.nan, sxy / sxx)
        d[h:n - h] = slope
        ends = list(range(0, h)) + list(range(n - h, n))
    else:
        ends = list(range(n))
    Vl = V.tolist()
    Il = I.tolist()
    for i in ends:
        lo, hi = window_bounds(i, n, h)
        _, b, _ = ols_fit(Vl[lo:hi + 1], Il[lo:hi + 1])
        d[i] = b
    return d


def interp_clamped(x, xp, fp):
    """Linear interpolation of (xp, fp) at x, clamped to fp[0]/fp[-1] outside
    [xp[0], xp[-1]].  xp must be ascending (ties allowed).  See docstring A14."""
    return np.interp(x, xp, fp)


# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------
def segment_boundaries(dac):
    """d[i] = sign(DAC[i+1]-DAC[i]); turning points T = {i+1 : d[i+1] != d[i]}; B = [0]+T+[N]."""
    N = len(dac)
    d = np.sign(np.diff(dac))
    T = [i + 1 for i in range(len(d) - 1) if d[i + 1] != d[i]]
    return [0] + T + [N]


def _q32(a):
    """Round a float64 array through float32 and back (the applet stores these
    columns as Float32Array; --quantize-float32 makes the reference see the
    same values)."""
    return a.astype(np.float32).astype(np.float64)


def load_csv(path, quantize_float32=False):
    df = pd.read_csv(path)
    ts = pd.to_datetime(df["timestamp"])
    t_rel = (ts - ts.iloc[0]).dt.total_seconds().to_numpy(dtype=float)
    data = {
        "DAC": df["DAC"].to_numpy(dtype=float),
        "V2": df["V2"].to_numpy(dtype=float),
        "I1": df["I1"].to_numpy(dtype=float),
        "I2": df["I2"].to_numpy(dtype=float),
        "t_rel": t_rel,
        "N": len(df),
    }
    if quantize_float32:
        # V2/I1/I2 and the seconds-since-first-timestamp are rounded to float32
        # right after loading; DAC stays as is (integers, exact in Uint16).
        for key in ("V2", "I1", "I2", "t_rel"):
            data[key] = _q32(data[key])
    return data


# ---------------------------------------------------------------------------
# Step 0a/0b helpers
# ---------------------------------------------------------------------------
def screen_overrange(I2seg):
    """Step 0a: row i is bad if I2[i] < -5e-8, or |I2[i]-I2[i-1]| > 1e-6 (i > r0),
    or |I2[i+1]-I2[i]| > 1e-6 (i < r1-1).  Time order.  NaN never compares true."""
    bad = I2seg < OVERRANGE_NEG
    jump = np.abs(I2seg[1:] - I2seg[:-1]) > OVERRANGE_JUMP
    bad[1:] |= jump
    bad[:-1] |= jump
    return bad


def select_rows(data, r0, r1, channel):
    """Step 0a + 0b (without sorting): returns (keep mask, bad mask) over rows r0..r1-1."""
    bad = screen_overrange(data["I2"][r0:r1])
    finite = np.isfinite(data[channel][r0:r1])
    keep = (~bad) & finite
    return keep, bad


def ref_baseline(data, B, ref_from, ref_to, channel):
    """Step 0c mode 'ref': pool all kept rows of segments ref_from..ref_to (time order),
    OLS I = a_r + b_r*V.  Returns (a_r, b_r, n_pool); a_r/b_r are None if n_pool < 50."""
    nseg = len(B) - 1
    Vp, Ip = [], []
    for k in range(ref_from, ref_to + 1):
        if k < 0 or k >= nseg:
            continue
        r0, r1 = B[k], B[k + 1]
        keep, _ = select_rows(data, r0, r1, channel)
        Vp.extend(data["V2"][r0:r1][keep].tolist())
        Ip.extend(data[channel][r0:r1][keep].tolist())
    if len(Vp) < 50:
        return None, None, len(Vp)
    a_r, b_r, _ = ols_fit(Vp, Ip)
    return a_r, b_r, len(Vp)


# ---------------------------------------------------------------------------
# Per-segment analysis (steps 0-10 + Te_local of step 11)
# ---------------------------------------------------------------------------
def analyze_segment(data, B, k, params, baseline, dump=False):
    r0, r1 = B[k], B[k + 1]
    channel = params["channel"]
    w = params["window"]
    h = (w - 1) >> 1
    f_ion = params["f_ion"]
    te_lo_f, te_hi_f = params["te_lo"], params["te_hi"]
    te_manual = params["te_manual"]
    A = params["area"]
    M = params["mass"]

    flags = []
    V2seg = data["V2"][r0:r1]
    I2seg = data["I2"][r0:r1]
    Iseg = data[channel][r0:r1]
    t_rel = data["t_rel"][r0:r1]
    direction = "up" if V2seg[-1] > V2seg[0] else "down"          # A15

    out = {key: None for key in SCALAR_KEYS}
    out.update({"k": k, "dir": direction,
                "t_start": float(t_rel[0]), "duration": float(t_rel[-1] - t_rel[0])})
    state = {"full": False, "overrange": False, "dir": direction, "V": None, "Is": None}   # for hysteresis
    dmp = {key: None for key in DUMP_KEYS} if dump else None

    def finish():
        out["flags"] = flags
        if baseline["unavailable"]:
            out["flags"] = ["baseline_unavailable"] + flags                 # A3
        if dump:
            out["dump"] = dmp
        return out, state

    # ---- Step 0a: overrange screening (time order, I2) ----------------------
    keep, bad = select_rows(data, r0, r1, channel)
    n_glitch = int(bad.sum())
    out["n_glitch"] = n_glitch
    if n_glitch > 0:
        flags.append("overrange_event")
        state["overrange"] = True
    if dump:
        bad_idx = np.nonzero(bad)[0]
        dmp["excluded_rows"] = [int(r0 + i) for i in bad_idx]
        dmp["excluded_VI"] = [[float(V2seg[i]), float(I2seg[i])] for i in bad_idx]   # A16

    # ---- Step 0b: select, drop NaN, sort -------------------------------------
    V = V2seg[keep]
    I = Iseg[keep]
    order = np.argsort(V, kind="stable")
    V = V[order].astype(float)
    I = I[order].astype(float)
    n = len(V)
    out["n"] = n
    out["n_dropped"] = int((r1 - r0) - n)                                    # A2
    if n < 50:
        flags.append("insufficient")
        return finish()                                                      # A4

    # ---- Step 0c: baseline -----------------------------------------------------
    if baseline["mode"] == "const":
        I = I - baseline["value"]
    elif baseline["mode"] == "ref":
        I = I - (baseline["a_r"] + baseline["b_r"] * V)
    # mode 'none': nothing

    Vmin = float(V[0])
    Vmax = float(V[-1])
    out["Vmin"], out["Vmax"] = Vmin, Vmax

    # ---- Steps 1-2: smoothing, derivative, noise ------------------------------
    Is = smooth_mean(I, h)
    dIdV = window_slope(V, I, h)
    lo, hi = h, n - 1 - h
    resid = I - Is
    sigma_I = math.sqrt(seq_sum((resid * resid).tolist()) / n)
    dV = (Vmax - Vmin) / (n - 1)
    sigma_b = sigma_I / math.sqrt(dV * dV * (w ** 3 - w) / 12.0)
    out["sigma_I"], out["sigma_b"] = sigma_I, sigma_b
    out["I_lo"] = float(Is[lo])
    out["I_hi"] = float(Is[hi])
    state["V"], state["Is"] = V, Is
    if dump:
        dmp["V"], dmp["I"], dmp["Is"], dmp["dIdV"] = V.tolist(), I.tolist(), Is.tolist(), dIdV.tolist()
        dmp["idx"] = {"lo": lo, "hi": hi, "ip": None, "iRef": None, "iVf": None}

    # ---- Step 2b: coverage -------------------------------------------------------
    if Vmin > -8.0 or Vmax < 8.0:
        flags.append("partial_sweep")
        return finish()                                                      # A5
    state["full"] = True

    # ---- Step 3: ion-saturation fit ----------------------------------------------
    ion_mask = V <= Vmin + f_ion * (Vmax - Vmin)
    if dump:
        dmp["ion_mask"] = ion_mask.tolist()
    if int(ion_mask.sum()) < 5:
        flags.append("ion_fit_insufficient")
        return finish()
    a_i, b_i, _ = ols_fit(V[ion_mask].tolist(), I[ion_mask].tolist())
    out["a_i"], out["b_i"] = a_i, b_i
    Iion = a_i + b_i * V

    # ---- Step 4: electron current ---------------------------------------------------
    Ie = I - Iion
    Ie_s = Is - Iion
    if dump:
        dmp["Ie"], dmp["Ie_s"], dmp["Iion"] = Ie.tolist(), Ie_s.tolist(), Iion.tolist()

    # ---- Step 5: floating potential ---------------------------------------------------
    crossings = []
    for i in range(n - 1):
        if Is[i] < 0 and Is[i + 1] >= 0:
            crossings.append(float(V[i] + (0.0 - Is[i]) * (V[i + 1] - V[i]) / (Is[i + 1] - Is[i])))
    if dump:
        dmp["crossings"] = list(crossings)
    Vf = None
    dVf = None
    if len(crossings) == 0:
        flags.append("no_zero_crossing")
    elif len(crossings) == 1:
        Vf = crossings[0]
    else:
        flags.append("vf_multiple_crossings")
        Vc_m = crossings[(len(crossings) - 1) // 2]
        near = np.abs(V - Vc_m) <= 0.5
        a_l, b_l, _ = ols_fit(V[near].tolist(), Is[near].tolist())
        if b_l > 0:                                                          # A8 (NaN -> False)
            Vf = float(-a_l / b_l)
        else:
            Vf = Vc_m
    iVf = None
    if Vf is not None:
        iVf = int(np.argmin(np.abs(V - Vf)))                                 # A7
        if dIdV[iVf] > 0:
            dVf = float(sigma_I / dIdV[iVf])
    out["Vf"], out["dVf"] = Vf, dVf

    # ---- Step 6: signal gate and plasma potential -------------------------------------
    ip = lo + int(np.argmax(dIdV[lo:hi + 1]))
    dIdV_max = float(dIdV[ip])
    out["dIdV_max"] = dIdV_max
    span = float(Is[hi] - Is[lo])
    if dump:
        dmp["idx"]["ip"] = ip
        dmp["idx"]["iVf"] = iVf
    span_fail = span < 10.0 * sigma_I
    slope_fail = dIdV_max < GATE_DERIV_FACTOR * sigma_b                      # v2.3 amendment F
    out["gate_reason"] = "span" if span_fail else ("slope" if slope_fail else None)   # v2.4 amendment I
    if span_fail or slope_fail:
        for f in ("no_zero_crossing", "vf_multiple_crossings"):              # v2.3 amendment G
            if f in flags:
                flags.remove(f)
        flags.append("no_plasma_signal")
        out["Vf"], out["dVf"] = None, None                                   # A6
        if dump:
            dmp["idx"]["iVf"] = None
        state["full"] = False
        return finish()

    Vp = None
    if n_glitch > 0:                                                         # v2.1 amendment A
        flags.append("top_of_sweep_unreliable")
    else:
        # v2.4 amendment H: flattening test, all four conditions must hold (A17)
        cond_a = bool(V[ip] <= Vmax - 1.0)
        # v2.4b amendment K: baseline-invariant growth limit (relative to the span below the knee)
        cond_b = bool(Is[hi] >= Is[ip] - 3.0 * sigma_I
                      and Is[hi] - Is[ip] <= 0.4 * (Is[ip] - Is[lo]))
        tail = [float(dIdV[j]) for j in range(ip + 1, hi + 1) if V[j] >= V[ip] + 0.3]
        tail_median = median_plain(tail) if tail else None
        cond_c = bool(tail_median is not None and tail_median <= 0.7 * dIdV[ip])
        cond_d = bool(Vf is None or V[ip] > Vf + 1.0)
        if cond_a and cond_b and cond_c and cond_d:
            Vp = float(V[ip])
        else:
            flags.append("vp_beyond_range")
            flags.append("not_saturated")
        if dump:
            dmp["idx"]["vp_test"] = {"a": cond_a, "b": cond_b, "c": cond_c, "d": cond_d,
                                     "tail_median": tail_median, "n_tail": len(tail)}
    out["Vp"] = Vp

    # ---- Step 7: electron temperature ----------------------------------------------------
    # v2.3 amendment E: iRef = ip if Vp accepted, else first maximum of Ie_s over [lo, hi]
    iRef = ip if Vp is not None else lo + int(np.argmax(Ie_s[lo:hi + 1]))
    Ie_ref = float(Ie_s[iRef])
    three_sigma = 3.0 * sigma_I
    if te_manual is None:                                                    # A11
        j = iRef
        while j >= 0 and Ie_s[j] > te_hi_f * Ie_ref:
            j -= 1
        te_idx = []
        while j >= 0 and (Ie_s[j] >= te_lo_f * Ie_ref and Ie_s[j] > three_sigma
                          and (Vf is None or V[j] > Vf)):
            te_idx.append(j)
            j -= 1
        te_idx.reverse()
    else:
        vlo_m, vhi_m = te_manual
        te_idx = [j for j in range(n) if vlo_m <= V[j] <= vhi_m and Ie_s[j] > three_sigma]
    te_mask = np.zeros(n, dtype=bool)
    te_mask[te_idx] = True
    m_te = len(te_idx)
    Te = None
    Te_r2 = None
    te_alpha = te_beta = None
    Te_lo = Te_hi = None
    out["Te_npts"] = m_te                                                    # A9
    out["Te_window"] = [float(V[te_idx[0]]), float(V[te_idx[-1]])] if m_te > 0 else None
    if m_te < 8:
        flags.append("te_insufficient")
    else:
        xs = [float(V[j]) for j in te_idx]
        ys = [math.log(float(Ie_s[j])) for j in te_idx]
        te_alpha, te_beta, Te_r2 = ols_fit(xs, ys)
        if te_beta > 0:
            Te = 1.0 / te_beta
        else:
            flags.append("te_negative_slope")
        if Te_r2 is not None and math.isfinite(Te_r2) and Te_r2 < 0.9:
            flags.append("te_poor_fit")                                      # A10
        # split diagnostic
        half = m_te // 2
        first, second = te_idx[:half], te_idx[half:]
        if len(first) >= 4 and len(second) >= 4:
            _, b_lo, _ = ols_fit([float(V[j]) for j in first], [math.log(float(Ie_s[j])) for j in first])
            _, b_hi, _ = ols_fit([float(V[j]) for j in second], [math.log(float(Ie_s[j])) for j in second])
            Te_lo = 1.0 / b_lo if b_lo > 0 else None
            Te_hi = 1.0 / b_hi if b_hi > 0 else None
            if Te_lo is not None and Te_hi is not None and Te_hi / Te_lo > 1.3:
                flags.append("te_nonexponential")
    out["Te"], out["Te_r2"], out["te_alpha"], out["te_beta"] = Te, Te_r2, te_alpha, te_beta
    out["Te_lo"], out["Te_hi"] = Te_lo, Te_hi
    if dump:
        dmp["te_mask"] = te_mask.tolist()
        dmp["idx"]["iRef"] = iRef

    # ---- Step 8: electron saturation current ----------------------------------------------
    Ies = None if n_glitch > 0 else float(Ie_s[iRef])                       # v2.1 amendment A
    out["Ies"] = Ies

    # ---- Step 9: ion saturation current ----------------------------------------------------
    vmin_mask = V <= Vmin + 1.0
    Iis_vmin = seq_sum(Is[vmin_mask].tolist()) / int(vmin_mask.sum())
    Iis_vf = float(a_i + b_i * Vf) if Vf is not None else None
    iis_below_noise = abs(Iis_vmin) < three_sigma
    if iis_below_noise:
        flags.append("iis_below_noise")
    if Iis_vmin >= 0:
        flags.append("ion_current_positive")
    out["Iis"], out["Iis_vf"] = Iis_vmin, Iis_vf

    # ---- Step 10: densities ----------------------------------------------------------------
    ne = ni = lambdaD = Vp_expected = None
    if Te is not None:
        if Ies is not None and Ies > 0:
            ne = Ies / (E_CHARGE * A * math.sqrt(E_CHARGE * Te / (2.0 * math.pi * M_E)))
            lambdaD = math.sqrt(EPS0 * Te / (ne * E_CHARGE))
        if Iis_vmin < 0 and not iis_below_noise:
            ni = abs(Iis_vmin) / (0.61 * E_CHARGE * A * math.sqrt(E_CHARGE * Te / (M * AMU)))
        if Vf is not None:
            Vp_expected = Vf + Te * math.log(math.sqrt(M * AMU / (2.0 * math.pi * M_E)) / 0.61)
    out["ne"], out["ni"], out["lambdaD"], out["Vp_expected"] = ne, ni, lambdaD, Vp_expected

    # ---- Step 11a: local slope temperature (array, for plotting) ----------------------------
    if dump:                                                                 # A12
        Te_local = []
        for i in range(n):
            den = dIdV[i] - b_i
            if Ie_s[i] > three_sigma and den > 0:
                Te_local.append(float(Ie_s[i] / den))
            else:
                Te_local.append(None)
        dmp["Te_local"] = Te_local

    return finish()


# ---------------------------------------------------------------------------
# Step 11b: hysteresis (v2.1 amendment B)
# ---------------------------------------------------------------------------
HYST_UNRESOLVED = 1.98


def _rms_curve(x, y, VD, IsD, H):
    """rms(h) = sqrt(mean_x (y - interp(VD, IsD, x + h))^2) for every h in H (sequential sums)."""
    cnt = len(x)
    D = np.empty((len(H), cnt))
    for i, h in enumerate(H):
        D[i] = y - interp_clamped(x + h, VD, IsD)
    return np.sqrt(seq_sum_rows(D * D) / cnt)


def _first_min(rms):
    best = 0
    for i in range(1, len(rms)):
        if rms[i] < rms[best]:
            best = i
    return best


def hysteresis_shift(VU, IsU, VD, IsD):
    """Shift h minimising rms over up-sweep points with -5 <= V <= 9 of
    (Is_up(V) - interp(V_down, Is_down, V + h)); positive h = the down-sweep
    reaches each current at a higher bias.  Coarse pass -2.0..2.0 step 0.1
    (41 values), fine pass h0-0.10..h0+0.10 step 0.02 (11 values, exact
    hundredths); first minimum wins in both passes.
    Returns (h, unresolved): h is None if no up points qualify; unresolved is
    True (and h None) when the fine-pass minimum has |h| >= 1.98."""
    m = (VU >= -5.0) & (VU <= 9.0)
    cnt = int(m.sum())
    if cnt == 0:
        return None, False
    x = VU[m]
    y = IsU[m]
    coarse = [(i - 20) / 10.0 for i in range(41)]
    i0 = _first_min(_rms_curve(x, y, VD, IsD, coarse))
    c0 = (i0 - 20) * 10                                   # h0 in exact hundredths
    fine = [(c0 + 2 * (j - 5)) / 100.0 for j in range(11)]
    j0 = _first_min(_rms_curve(x, y, VD, IsD, fine))
    h = fine[j0]
    if abs(h) >= HYST_UNRESOLVED:
        return None, True
    return h, False


def partner_of(k, nseg):
    """v2.3 amendment D: up-sweep k (even) pairs with the following down-sweep k+1;
    odd k pairs with k-1; a partner beyond the last segment becomes nSeg-2."""
    p = k + 1 if k % 2 == 0 else k - 1
    if p > nseg - 1:
        p = nseg - 2
    return p


def _eligible(st):
    return st["full"] and not st["overrange"]


def _pairable(states, k, p):
    return p != k and _eligible(states[k]) and _eligible(states[p]) and states[k]["dir"] != states[p]["dir"]


def compute_hysteresis(segments, states):
    nseg = len(segments)
    cache = {}
    for k in range(nseg):
        p = partner_of(k, nseg)
        symmetric = partner_of(p, nseg) == k                                  # v2.4 amendment J
        key = (min(k, p), max(k, p))
        if symmetric and key in cache:
            val, unresolved = cache[key]                                      # copy branch
        elif _pairable(states, k, p):
            u, d = (k, p) if states[k]["dir"] == "up" else (p, k)
            val, unresolved = hysteresis_shift(states[u]["V"], states[u]["Is"],
                                               states[d]["V"], states[d]["Is"])
            if symmetric:
                cache[key] = (val, unresolved)
        else:
            val, unresolved = None, False
        if unresolved:
            for member in ((k, p) if symmetric else (k,)):
                if "hysteresis_unresolved" not in segments[member]["flags"]:
                    segments[member]["flags"].append("hysteresis_unresolved")
        segments[k]["hysteresis"] = val


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run(csv_path, params, dump_segments=()):
    data = load_csv(csv_path, quantize_float32=params.get("quantize_float32", False))
    B = segment_boundaries(data["DAC"])
    nseg = len(B) - 1
    baseline = {"mode": params["baseline"], "value": params["baseline_value"],
                "a_r": None, "b_r": None, "unavailable": False, "n_pool": None}
    if params["baseline"] == "ref":
        a_r, b_r, n_pool = ref_baseline(data, B, params["baseline_ref"][0], params["baseline_ref"][1],
                                        params["channel"])
        baseline["n_pool"] = n_pool
        if a_r is None:
            baseline["mode"] = "none"
            baseline["unavailable"] = True
        else:
            baseline["a_r"], baseline["b_r"] = a_r, b_r
    segments, states = [], []
    for k in range(nseg):
        seg, st = analyze_segment(data, B, k, params, baseline, dump=(k in dump_segments))
        segments.append(seg)
        states.append(st)
    compute_hysteresis(segments, states)
    return {
        "params": dict(params, n_segments=nseg, boundaries=B,
                       baseline_effective=baseline["mode"], ref_pool_n=baseline["n_pool"]),
        "a_r": baseline["a_r"], "b_r": baseline["b_r"],
        "segments": segments,
    }


def _json_clean(obj):
    """Replace NaN/inf with None recursively so the JSON is strict."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _json_clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_clean(v) for v in obj]
    if isinstance(obj, np.floating):
        return _json_clean(float(obj))
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Langmuir sweep analysis (spec v2 reference implementation)")
    p.add_argument("csv")
    p.add_argument("--channel", choices=["I2", "I1"], default="I2")
    p.add_argument("--baseline", default="ref:74,79",
                   help="none | const:<A> | ref:<from>,<to> (0-based inclusive segment range; "
                        "'ref' alone = ref:74,79; ref:<k> = ref:k,k)")
    p.add_argument("--f-ion", type=float, default=0.20)
    p.add_argument("--window", type=int, default=31, help="odd smoothing window (points)")
    p.add_argument("--te-lo", type=float, default=0.02)
    p.add_argument("--te-hi", type=float, default=0.30)
    p.add_argument("--te-manual", default=None, help="Vlo,Vhi manual Te window (volts)")
    p.add_argument("--area", type=float, default=1.0e-5, help="probe area [m^2]")
    p.add_argument("--mass", type=float, default=39.948, help="ion mass [amu]")
    p.add_argument("--out", default="results.json")
    p.add_argument("--dump-segments", default="", help="comma-separated segment indices to dump full arrays for")
    p.add_argument("--csv-out", default=None, help="write one row per segment with the scalar outputs")
    p.add_argument("--quantize-float32", action="store_true",
                   help="round V2/I1/I2 and the relative timestamps through float32 right after "
                        "loading (matches the applet's Float32Array columns)")
    a = p.parse_args(argv)

    if a.window < 3 or a.window % 2 == 0:
        p.error("--window must be an odd integer >= 3")
    params = {
        "csv": a.csv, "channel": a.channel, "f_ion": a.f_ion, "window": a.window,
        "te_lo": a.te_lo, "te_hi": a.te_hi, "area": a.area, "mass": a.mass,
        "baseline": "none", "baseline_value": None, "baseline_ref": None, "te_manual": None,
        "quantize_float32": bool(a.quantize_float32),
    }
    b = a.baseline.strip()
    if b == "none":
        pass
    elif b.startswith("const:"):
        params["baseline"] = "const"
        params["baseline_value"] = float(b[len("const:"):])
    elif b == "ref" or b.startswith("ref:"):
        params["baseline"] = "ref"
        if b == "ref":
            params["baseline_ref"] = [74, 79]
        else:
            parts = [int(s) for s in b[len("ref:"):].split(",") if s.strip() != ""]
            if len(parts) == 1:
                parts = [parts[0], parts[0]]
            if len(parts) != 2 or parts[0] > parts[1] or parts[0] < 0:
                p.error("--baseline ref:<from>,<to> needs 0 <= from <= to")
            params["baseline_ref"] = parts
    else:
        p.error("--baseline must be none | const:<A> | ref:<from>,<to>")
    if a.te_manual:
        parts = [float(s) for s in a.te_manual.split(",")]
        if len(parts) != 2 or parts[0] >= parts[1]:
            p.error("--te-manual must be Vlo,Vhi with Vlo < Vhi")
        params["te_manual"] = parts
    dump = set()
    if a.dump_segments:
        dump = {int(s) for s in a.dump_segments.split(",") if s.strip() != ""}
    return a, params, dump


def write_csv(result, path):
    rows = []
    for s in result["segments"]:
        row = {}
        for key in SCALAR_KEYS:
            v = s.get(key)
            if key == "Te_window":
                row["Te_Vlo"] = None if v is None else v[0]
                row["Te_Vhi"] = None if v is None else v[1]
            elif key == "flags":
                row["flags"] = ";".join(v or [])
            else:
                row[key] = v
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def main(argv=None):
    a, params, dump = parse_args(argv)
    result = run(a.csv, params, dump)
    result = _json_clean(result)
    with open(a.out, "w") as f:
        json.dump(result, f, indent=1)
    if a.csv_out:
        write_csv(result, a.csv_out)
    nseg = result["params"]["n_segments"]
    msg = f"analysed {nseg} segments (baseline {result['params']['baseline_effective']}"
    if result["a_r"] is not None:
        msg += f", a_r={result['a_r']:.6g} b_r={result['b_r']:.6g}"
    msg += f") -> {a.out}" + (f", {a.csv_out}" if a.csv_out else "")
    print(msg)


if __name__ == "__main__":
    main()
