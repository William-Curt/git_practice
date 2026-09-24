# Langmuir sweep analysis — algorithm spec v2 (supersedes v1; unchanged parts are restated)

Constants: e=1.602176634e-19, me=9.1093837015e-31, amu=1.66053906660e-27, eps0=8.8541878128e-12.
Defaults: channel=I2; baseline={mode:'ref', refFrom:74, refTo:79} (0-based segment indices, inclusive) | {mode:'const', value:1.48e-7} | {mode:'none'};
f_ion=0.20; w=31; teMode='auto' with teLo=0.02, teHi=0.30 (fractions of Ie_ref) or 'manual' [Vlo,Vhi]; area A=1.0e-5 m²; mass M=39.948 amu.

## Segmentation (unchanged)
d[i]=sign(DAC[i+1]−DAC[i]); turning points where d changes; bounds=[0]+turns+[N]. Direction 'up' if V2 rises through the segment.

## Step 0a — overrange screening (TIME ORDER, before anything else; uses I2 regardless of channel)
Within the segment rows r0..r1−1 (time order), row i is `bad` if I2[i] < −5e-8, or (i>r0 and |I2[i]−I2[i−1]| > 1e-6), or (i<r1−1 and |I2[i+1]−I2[i]| > 1e-6).
n_glitch = count(bad). If n_glitch>0 flag 'overrange_event'. Bad rows are excluded from ALL further steps (they are kept only for plotting).
## Step 0b — select, drop NaN, sort
Keep rows that are not bad and have finite channel value. n = count; n_dropped = (r1−r0) − n. Stable sort ascending by V2 (ties keep row order). If n<50 → flag 'insufficient', stop.
## Step 0c — baseline
mode 'const': I := I − value.  mode 'none': nothing.
mode 'ref': pool every row of segments refFrom..refTo that is not bad and has a finite channel value; OLS fit I_ref = a_r + b_r·V over the pool; I := I − (a_r + b_r·V). Report a_r, b_r once (top level). If the pool has <50 points, fall back to mode 'none' and flag 'baseline_unavailable' on every segment.

## Step 1–2 — smoothing, derivative, noise (unchanged formulas)
h=(w−1)>>1; Is[i]=mean(I[max(0,i−h)..min(n−1,i+h)]); dIdV[i]=OLS slope of I vs V over the same clipped window. lo=h, hi=n−1−h.
sigma_I = sqrt(mean((I−Is)²)) over all n points.  dV=(Vmax−Vmin)/(n−1).  sigma_b = sigma_I / sqrt(dV²·(w³−w)/12).
Report sigma_I, sigma_b, Vmin, Vmax, I_lo=Is[lo], I_hi=Is[hi].
## Step 2b — coverage
If Vmin > −8 or Vmax < +8 → flag 'partial_sweep'; skip steps 3–11 (all physics outputs null) but still return Is, dIdV, sigma.

## Step 3 — ion-saturation fit (unchanged)
ionMask: V ≤ Vmin + f_ion·(Vmax−Vmin); need ≥5 points else flag 'ion_fit_insufficient' and stop. OLS on raw I → a_i, b_i. Iion(V)=a_i+b_i·V.
## Step 4 — electron current (unchanged)
Ie = I − Iion(V); Ie_s = Is − Iion(V).
## Step 5 — floating potential
crossings = all i in 0..n−2 with Is[i] < 0 and Is[i+1] ≥ 0, each with Vc_i = V[i] + (0−Is[i])·(V[i+1]−V[i])/(Is[i+1]−Is[i]).
0 crossings → Vf=null, flag 'no_zero_crossing'.  1 crossing → Vf=Vc.
>1 crossings → flag 'vf_multiple_crossings'; m = crossings[floor((count−1)/2)] (the lower-middle one); OLS fit Is vs V over points with |V − Vc_m| ≤ 0.5; if slope b>0 then Vf = −a/b else Vf = Vc_m.
If Vf found: iVf = index of the point with V closest to Vf; dVf = sigma_I / dIdV[iVf] if dIdV[iVf] > 0 else null.
## Step 6 — signal gate and plasma potential
ip = argmax dIdV over [lo,hi] (first maximum). dIdV_max = dIdV[ip]. span = Is[hi] − Is[lo].
If span < 10·sigma_I or dIdV_max < 6·sigma_b → flag 'no_plasma_signal'; set Vf=null (and dVf), Vp, Te, Ies, Iis, ne, ni, lambdaD, Vp_expected, hysteresis = null; skip steps 7–11.
Vp accepted iff V[ip] ≤ Vmax − 1.0 AND there exists j in (ip, hi] with dIdV[j] ≤ 0.7·dIdV[ip]. Then Vp = V[ip]. Otherwise Vp=null, flags 'vp_beyond_range' and 'not_saturated'.
## Step 7 — electron temperature
iRef = ip if Vp found else hi. Ie_ref = Ie_s[iRef].
auto: walk j from iRef downward to 0: first skip while Ie_s[j] > teHi·Ie_ref; then collect j while (Ie_s[j] ≥ teLo·Ie_ref and Ie_s[j] > 3·sigma_I and (Vf null or V[j] > Vf)); stop at the first j that fails. teMask = collected indices (contiguous).
manual: teMask = { j : Vlo ≤ V[j] ≤ Vhi and Ie_s[j] > 3·sigma_I }.
If |teMask| < 8 → Te=null, flag 'te_insufficient'. Else OLS ln(Ie_s) = alpha + beta·V over teMask; Te = 1/beta if beta>0 else null + flag 'te_negative_slope'; Te_r2 = R²; if Te_r2 < 0.9 flag 'te_poor_fit'. Te_window = [min V, max V] over teMask; Te_npts.
Split diagnostic: sort teMask ascending by V; first half (floor(m/2) points) and second half (the rest); if each half has ≥4 points fit OLS separately → Te_lo=1/beta_lo, Te_hi=1/beta_hi (null if slope ≤0). If both exist and Te_hi/Te_lo > 1.3 → flag 'te_nonexponential'.
## Step 8 — electron saturation current
Ies = Ie_s[iRef]. (When Vp is null it is a lower bound: 'not_saturated' already set.)
## Step 9 — ion saturation current
Iis_vmin = mean of Is over points with V ≤ Vmin + 1.0 (primary; report as Iis). Iis_vf = Iion(Vf) if Vf found else null.
If |Iis_vmin| < 3·sigma_I → flag 'iis_below_noise'. If Iis_vmin ≥ 0 → flag 'ion_current_positive'.
## Step 10 — densities
ne = Ies/(e·A·sqrt(e·Te/(2π·me))) if Te and Ies>0; lambdaD = sqrt(eps0·Te/(ne·e)).
ni = |Iis_vmin|/(0.61·e·A·sqrt(e·Te/(M·amu))) only if Te and Iis_vmin<0 and not 'iis_below_noise'.
Vp_expected = Vf + Te·ln( sqrt(M·amu/(2π·me)) / 0.61 ) if Vf and Te.
## Step 11 — local slope temperature and hysteresis
Te_local[i] = Ie_s[i]/(dIdV[i] − b_i) where Ie_s[i] > 3·sigma_I and dIdV[i] − b_i > 0, else null (array, for plotting).
Hysteresis (per pair; partner(k) = k+1 if k odd else k−1, clamped to [0, nSeg−1], partner(0)=1): only if both members are full (not partial_sweep, not insufficient, not no_plasma_signal): let U be the up member, D the down member (by dir). For s in −2.0..+2.0 step 0.02: rms(s) = sqrt(mean over U points with −5 ≤ V ≤ 9 of (Is_U(V) − interp(V_D + s, Is_D, V))²) (interp linear, clamped, V_D+s ascending). hysteresis = s minimizing rms (first minimum). Same value stored on both members. Positive s means the down-sweep must be shifted to higher V to align with the up-sweep.

Outputs per segment: {k, dir, t_start, duration, n, n_dropped, n_glitch, Vmin, Vmax, sigma_I, sigma_b, a_i, b_i, Vf, dVf, Vp, Te, Te_r2, Te_window, Te_npts, Te_lo, Te_hi, te_alpha, te_beta, Ies, Iis (=Iis_vmin), Iis_vf, ne, ni, lambdaD, Vp_expected, I_lo, I_hi, dIdV_max, hysteresis, flags}. Top level: {a_r, b_r} for the ref baseline.

## v2.1 amendments (supersede the corresponding text above)
A. Over-range sweeps. If a segment has 'overrange_event' (n_glitch>0), the top of the sweep is untrustworthy: after step 6's signal gate, set Vp=null, Ies=null, ne=null, lambdaD=null, hysteresis=null and add flag 'top_of_sweep_unreliable' INSTEAD of 'vp_beyond_range'/'not_saturated' (do not evaluate the Vp test at all). Te, Te window, Vf, Iis, ni, Vp_expected are still computed as specified (iRef = hi for the Te window since Vp is null).
B. Hysteresis definition and search. h is the shift that minimises rms(h) = sqrt(mean over up-sweep points with −5 ≤ V ≤ 9 of (Is_up(V) − interp(V_down, Is_down, V + h))²), i.e. the down-sweep curve is sampled at V + h. Positive h means the down-sweep reaches each current at a higher bias than the up-sweep (a lag). Search: coarse pass h = −2.0, −1.9, …, +2.0 (step 0.1, 41 values, first minimum wins); fine pass h = h0−0.10, h0−0.08, …, h0+0.10 (step 0.02, 11 values, first minimum wins; values are rounded to 2 decimals to avoid float drift). If the fine-pass minimum satisfies |h| ≥ 1.98 → hysteresis=null and flag 'hysteresis_unresolved' on both members. A pair is only evaluated if neither member has 'overrange_event', 'partial_sweep', 'insufficient' or 'no_plasma_signal', and the members have opposite directions.
C. (v2.2) Vp acceptance gains a third condition: Is[hi] ≥ Is[ip] − 3·sigma_I (the current must not drop after the knee). All three must hold: V[ip] ≤ Vmax − 1.0, a later dIdV[j] ≤ 0.7·dIdV[ip] for some j in (ip, hi], and Is[hi] ≥ Is[ip] − 3·sigma_I. Otherwise 'vp_beyond_range' + 'not_saturated' as before.
D. (v2.3) Pairing: partner(k) = k+1 for even k, k−1 for odd k (an up-sweep with the down-sweep that follows it); if partner > nSeg−1 use nSeg−2. Everything else about the pair step is unchanged.
E. (v2.3) Step 7 reference index: iRef = ip if Vp was accepted; otherwise iRef = the index of the first maximum of Ie_s over [lo, hi]. Ie_ref = Ie_s[iRef]; Step 8 Ies = Ie_s[iRef] as before.
F. (v2.3) Signal gate derivative threshold is 5·sigma_b (was 6·sigma_b). The span criterion (10·sigma_I) is unchanged.
G. (v2.3) When the gate fires 'no_plasma_signal' it also removes any 'no_zero_crossing' and 'vf_multiple_crossings' flags set in step 5 (Vf/dVf are null as before).
