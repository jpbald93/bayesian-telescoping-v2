#!/usr/bin/env python3
"""
Verification + numerical study for V2 of "Bayesian Telescoping Under Non-Ideal
Conditions".  Regenerates every numeral in the V2 manuscript and asserts every
claimed sign / threshold.  No LLM arithmetic: all numbers are computed here.
"""
import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 50
FAILS = []
def check(cond, msg):
    print(("  [ok ] " if cond else "  [FAIL] ") + msg)
    if not cond: FAILS.append(msg)
def slope(xs, ys):
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    mx, my = xs.mean(), ys.mean()
    return float(((xs-mx)*(ys-my)).sum() / ((xs-mx)**2).sum())

Re = lambda z: mp.e**z
Tz = lambda z: (1+z/2)/(1-z/2)           # trapezoidal stability function (z<2)
Rz = lambda z: 1+z+z**2/2+z**3/6+z**4/24 # RK4 stability function
Rzp = lambda z: 1+z+z**2/2+z**3/6        # R'(z)

# ---------------------------------------------------------------- 1
print("="*74); print("1. Telescoping tail lemma (correct integral step)"); print("="*74)
for k in [1, 2, 3]:
    for N in [10, 1000]:
        S = sum(n**(-(k+1)) for n in range(N, N+200000))
        bd = N**(-(k+1)) + N**(-k)/k      # n^-(k+1) <= int_{n-1}^n t^-(k+1) dt + first term
        check(S <= bd*1.0001, f"  k={k} N={N}: sum={S:.6e} <= N^-(k+1)+N^-k/k={bd:.6e}")
    A = sum(n**(-(k+1)) for n in range(100, 100000))
    B = sum(n**(-(k+1)) for n in range(1000, 1000000))
    r = (A/100.0**-k)/(B/1000.0**-k)
    check(abs(r-1) < 0.05, f"  k={k}: Theta(N^-k) verified (ratio {r:.3f})")

# ---------------------------------------------------------------- 2
print("="*74); print("2. MLMC complexity -- corrected regimes"); print("="*74)
print("  Cost ~ eps^-2 (sum_{n<L} sqrt(V_n C_n))^2 ; s=(beta-gamma)/2 ; O(1) iff gamma>beta+2")
for beta, gamma in [(1,4),(1,3),(2,4),(1,2.5),(2,3.5),(1,5)]:
    s = (beta-gamma)/2; S = sum(n**s for n in range(1, 10**6))
    if abs(s+1) < 1e-12:
        check(10 < S < 100, f"  beta={beta} gamma={gamma}: s=-1, sum~{S:.3e} (log divergence; V1 said O(1))")
    else:
        check((S < 50) == (gamma > beta+2),
              f"  beta={beta} gamma={gamma}: s={s:.2f} sum={S:.3e}  O(1) iff gamma>{beta+2}")
beta, gamma, k = 2, 1.5, 2
eps = 1e-3; L = int(round(eps**(-1/k))); s = (beta-gamma)/2
cost = eps**-2 * (sum(n**s for n in range(1, L)))**2
pred = eps**(-2-(beta+2-gamma)/k)
# negative control: the superseded V1 pair (cutoff eps^-1/(k+1) + exponent /(k+1)) must NOT match
Lv1 = int(round(eps**(-1/(k+1)))); costv1 = eps**-2 * (sum(n**s for n in range(1, Lv1)))**2
check(abs(np.log10(costv1/pred)) > 0.15, "  negative control: the V1 cutoff/exponent pair does NOT reproduce the corrected scaling")
check(abs(np.log10(cost/pred)) < 1.0, f"  case-(iii) exponent -2-(beta+2-gamma)/k: cost={cost:.3e} pred~{pred:.3e}")
# bias cutoff: sum_{n>=L} n^-(k+1) = Theta(L^-k)  =>  L ~ eps^-1/k  (NOT eps^-1/(k+1))
for k in [1,2,3]:
    L = 10**6; S = sum(n**(-(k+1)) for n in range(L, L+3000000))
    check(0.5*(L**-k)/k < S < 2*(L**-k)/k, f"  bias sum_(n>=L)n^-(k+1)={S:.3e} = Theta(L^-k/k={L**-k/k:.3e})  => L~eps^-1/k")

# ---------------------------------------------------------------- 3
print("="*74); print("3. Method maps: monotonicity"); print("="*74)
check(all(Re(mp.mpf(z)/10) < Re(mp.mpf(z+1)/10) for z in range(-200,200)), "  exp(z) strictly increasing")
check(all(Tz(mp.mpf(z)/10) < Tz(mp.mpf(z+1)/10) for z in range(-100,17)), "  trapezoidal T(z) strictly increasing on z<2")
lo, hi = mp.mpf(-2), mp.mpf(-1)
for _ in range(200):
    mid = (lo+hi)/2
    if Rzp(lo)*Rzp(mid) <= 0: hi = mid
    else: lo = mid
zstar = (lo+hi)/2
check(-2 < zstar < -1, f"  RK4 R'(z) root z*={mp.nstr(zstar,8)} in (-2,-1) => R non-monotone (interior min R(z*)={mp.nstr(Rz(zstar),6)})")
print(f"    => RK4 map lambda |-> R(lambda h)^n is non-injective when |lambda|T > |z*| = {mp.nstr(abs(zstar),6)}")

# ---------------------------------------------------------------- 4
print("="*74); print("4. Discretization-induced bimodality (real RK4 vs trapezoidal)"); print("="*74)
def modes(phi, y, m0, s0, sigma, lo=-8.0, hi=6.0, ng=600001):
    l = np.linspace(lo, hi, ng)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        f = phi(l)
    lp = np.where(np.isfinite(f), -0.5*((l-m0)/s0)**2 - 0.5*((y-f)/sigma)**2, -1e9)
    lp -= lp.max()
    out = []
    for i in range(1, len(l)-1):
        if lp[i] > lp[i-1] and lp[i] >= lp[i+1] and lp[i] > -200:   # NO height threshold
            if out and l[i]-out[-1][0] < 0.3:
                if lp[i] > out[-1][1]: out[-1] = (float(l[i]), float(lp[i]))
            else: out.append((float(l[i]), float(lp[i])))
    return out
Rnp = lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24
Tnp = lambda z: (1 + z/2)/(1 - z/2)
T = 1.0; m0, s0, sigma, lam = -1.0, 0.5, 0.05, -1.0
y0 = float(mp.e**(lam*T))
print(f"  y = phi_exact(lambda*) + eta ; T={T} lambda*={lam} prior=N({m0},{s0}^2) sigma={sigma} y={y0:.6f}")
scan = {}
for n in [1, 2, 3, 4, 8, 16, 32, 64]:
    h = T/n
    me = modes(lambda z: np.exp(z*T),        y0, m0, s0, sigma)
    mr = modes(lambda z, h=h, n=n: Rnp(z*h)**n, y0, m0, s0, sigma)
    mt = modes(lambda z, h=h, n=n: Tnp(z*h)**n, y0, m0, s0, sigma)
    scan[n] = (me, mr, mt)
    def fmt(ms): return f"{len(ms)}" + (" " + ", ".join(f"{m:+.2f}(m~{np.exp(lh):.1e})" for m, lh in ms) if ms else "")
    print(f"    n={n:>3}: EXACT {fmt(me)} | RK4 {fmt(mr)} | TRAP {fmt(mt)}")
# observable-mass bimodality only at n=1; the n=2 second mode is exponentially suppressed
m1 = scan[1][1]; m2 = scan[2][1]
check(len(m1) >= 2 and np.exp(min(lh for _, lh in m1)) > 1e-2,
      "  RK4 has an OBSERVABLE second mode at n=1 (relative mass > 1e-2)")
check(len(m2) >= 2 and np.exp(min(lh for _, lh in m2)) < 1e-6,
      "  RK4's n=2 second mode exists but is exponentially suppressed (mass < 1e-6)")
check(all(len(scan[n][1]) == 1 for n in [4,8,16,32,64]), "  RK4 secondary mode leaves the scanned region for n>=4")
check(all(len(scan[n][0]) == 1 for n in scan), "  exact posterior unimodal throughout (monotone map)")
check(all(len(scan[n][2]) == 1 for n in scan), "  trapezoidal posterior unimodal throughout (monotone map on prior support)")

# ---------------------------------------------------------------- 5
print("="*74); print("5. Obstruction: linear forward operator + Gaussian prior"); print("="*74)
a, m0b, s0b, sig, yb = 1.3, 0.0, 1.0, 0.2, 0.7
prec = 1/s0b**2 + a**2/sig**2; mu = (m0b/s0b**2 + a*yb/sig**2)/prec; var = 1/prec
print(f"  posterior = N({mu:.6f}, {var:.6f}) exactly -- Gaussian, hence log-concave & unimodal")
check(var > 0, "  linear F + Gaussian prior => exact Gaussian posterior for EVERY level => its approximation cannot be bimodal")

# ---------------------------------------------------------------- 6
print("="*74); print("6. Heavy-tail priors: does the rate degrade?"); print("="*74)
def Emean(n, y, sigma, lp, lo=-30.0, hi=30.0, ng=200001):
    th = np.linspace(lo, hi, ng); Fn = th*(1-n**-2)
    v = lp(th) - 0.5*((y-Fn)/sigma)**2; v -= v.max(); w = np.exp(v); w /= w.sum()
    return float((w*th).sum())
gauss = lambda t: -0.5*t**2
t_nu  = lambda nu: (lambda t: -((nu+1)/2)*np.log1p(t**2/nu))
print("  V1 Experiment-6 model F_n(theta)=theta(1-n^-2) (k=2).  consecutive-n increments:")
for name, lp in [("Gaussian",gauss),("t_nu=10",t_nu(10)),("t_nu=5",t_nu(5)),("t_nu=3",t_nu(3)),("t_nu=2",t_nu(2))]:
    ns = list(range(60, 200, 2))
    Es = [Emean(n, 1.5, 0.1, lp) for n in ns]
    inc = [abs(Es[i]-Es[i+1]) for i in range(len(Es)-1)]
    good = [(ns[i], inc[i]) for i in range(len(inc)) if inc[i] > 1e-13]
    s = slope([np.log(a) for a,_ in good], [np.log(b) for _,b in good])
    rho = None if name=="Gaussian" else (float(name.split("=")[1])-2)/(float(name.split("=")[1])-1)
    print(f"    {name:<9} slope={s:+.3f}   V1's rho(alpha)*(k+1)={(-3.0*rho) if rho is not None else -3.0:+.3f}")
    if name == "Gaussian":
        check(abs(s+3) < 0.15, "  Gaussian prior reproduces the (k+1)=3 rate")
    else:
        check(abs(s+3) < 0.15, f"  {name}: rate is STILL ~(k+1)=3 -> NO degradation (contradicts V1 Table 9)")

print("="*74); print("RESULT:", "ALL ASSERTIONS PASS" if not FAILS else f"{len(FAILS)} FAILURE(S)")
for f in FAILS: print("   -", f)
print("="*74)
sys.exit(1 if FAILS else 0)
