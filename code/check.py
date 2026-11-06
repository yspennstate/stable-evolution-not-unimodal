"""Checks for the two counterexamples (classical and free, index 1/2).

Part 1: exact rational checks of the constants in the proofs.
Part 2: numerical density of mu * L_t (classical) at the three test points.
Part 3: numerical density of mu_0 boxplus T_t (free) across the gap, by
        Denjoy-Wolff iteration of the subordination map.
"""
from fractions import Fraction as Fr
import math, cmath

out = []
P = out.append

# ---------------- Part 1: exact constants ----------------
# classical: sup|H| bound, H(1/6) > 4, tail ratio, margin
P("Part 1 (exact)")
P("  (1/2)(7/2)^4 + (3/2)(5/2)^3 = %s < 99: %s" % (Fr(1, 2) * Fr(7, 2) ** 4 + Fr(3, 2) * Fr(5, 2) ** 3,
                                                  Fr(1, 2) * Fr(7, 2) ** 4 + Fr(3, 2) * Fr(5, 2) ** 3 < 99))
H16 = 54 * math.sqrt(6) * math.exp(-3)
P("  H(1/6) = 54 sqrt6 e^-3 = %.6f > 4: %s ; lower bound 108/27 = %s" % (H16, H16 > 4, Fr(108, 27)))
q = Fr(1, 2 ** 7)
P("  tail ratio at n=3: q/(1-q) = %s ; margin 4 - 99/127 = %s" % (q / (1 - q), 4 - Fr(99, 127)))
# check (1-3u)/(2u^{7/2}) e^{-1/(2u)} formula against numerical derivative of K_1
K = lambda r: r ** -1.5 * math.exp(-1 / (2 * r)) if r > 0 else 0.0
Hf = lambda u: (1 - 3 * u) / (2 * u ** 3.5) * math.exp(-1 / (2 * u))
for u in (0.1, 1 / 6, 0.5, 2.0):
    h = 1e-6
    P("  K_1'(%.4f): numeric %.8f formula %.8f" % (u, (K(u + h) - K(u - h)) / (2 * h), Hf(u)))
P("  sup over a grid of |H| = %.4f (< 99)" % max(abs(Hf(k / 1000)) for k in range(1, 200000)))
# free: W bound, w0, I2 constant
W_bound = Fr(1, 8) + Fr(1, 32 ** 2) / (1 - Fr(1, 32))
P("  W < 1/8 + 1/992 = %s < 1/4: %s" % (W_bound, W_bound < Fr(1, 4)))
P("  I2 bound 4 + 16/4 + 4/4 = %s (uses 16W<4, 4W<1)" % (4 + 4 + 1))
P("  sqrt(A) > 9t  =>  9t/(2 sqrt A) < 1/2")

# ---------------- Part 2: classical ----------------
P("")
P("Part 2 (classical, mu = sum Z^-1 2^{-n^2} delta_{-2^n}, density of mu * L_t up to the factor t/sqrt(2 pi))")
N = 12
Z = sum(2.0 ** (-n * n) for n in range(1, N + 1))
atoms = [(-(2.0 ** n), 2.0 ** (-n * n) / Z) for n in range(1, N + 1)]


def Kp(r, s):  # derivative of r^{-3/2} e^{-s/(2r)}
    return (s - 3 * r) / (2 * r ** 3.5) * math.exp(-s / (2 * r)) if r > 0 else 0.0


def dens(x, s):
    return sum(w * (x - a) ** -1.5 * math.exp(-s / (2 * (x - a))) for a, w in atoms if x > a)


def ddens(x, s):
    return sum(w * Kp(x - a, s) for a, w in atoms)


for t in (0.5, 1.0, 3.0, 10.0):
    s = t * t
    n = 3
    while not 2 ** (n - 1) > s:
        n += 1
    a = -(2.0 ** n)
    pts = [a, a + s / 6, a + s]
    d = [ddens(x, s) for x in pts]
    P("  t=%5.1f n=%d: f'(a_n)=%+.3e  f'(a_n+s/6)=%+.3e  f'(a_n+s)=%+.3e  pattern ok: %s"
      % (t, n, d[0], d[1], d[2], d[0] <= 0 and d[1] > 0 and d[2] < 0))
    # count sign changes of f' on a grid (local extrema)
    xs = [-(2.0 ** (n + 2)) + k * (2.0 ** (n + 2)) / 4000 for k in range(1, 4000)]
    sg = [ddens(x, s) for x in xs]
    mins = sum(1 for i in range(len(sg) - 1) if sg[i] < 0 <= sg[i + 1])
    P("          local minima of the density on (-2^{n+2}, 0): %d" % mins)

# ---------------- Part 3: free ----------------
P("")
P("Part 3 (free, mu_0 = w0 delta_0 + sum 2^-n A_n^-n delta_{-A_n}, A_n = 2^(2^n); T_t with phi(z) = -i t sqrt z)")
At = [2.0 ** (2 ** n) for n in range(1, 5)]
wt = [2.0 ** -n * A ** -n for n, A in zip(range(1, 5), At)]
w0 = 1 - sum(wt)
nu = [(0.0, w0)] + [(-A, w) for A, w in zip(At, wt)]
P("  weights: w0=%.6f, w_n=%s" % (w0, ["%.3e" % w for w in wt]))


def G(z):
    return sum(w / (z - x) for x, w in nu)


def omega(z, t, it=4000):
    w = z + 1j
    for _ in range(it):
        w = z + 1j * t * cmath.sqrt(1 / G(w))
    return w


for t in (1.0, 1.7):
    n = next(k for k, A in enumerate(At, start=1) if k >= 2 and math.sqrt(A) > 9 * t)
    A = At[n - 1]
    c = -A / 2
    a = -G(c).real
    I2 = sum(w / (c - x) ** 2 for x, w in nu)
    Hp = 1 - t * I2 / (2 * a ** 1.5)
    P("  t=%.1f n=%d A=%g: -G(c)=%.4e (>1/A=%.4e), I2=%.4e (<9/A^2=%.4e), H'(c)=%.6f"
      % (t, n, A, a, 1 / A, I2, 9 / A ** 2, Hp))
    hc = c + t / math.sqrt(a)  # H(c) = c - i t sqrt(F(c)), sqrt(F(c)) = i/sqrt(a)
    P("    H(c) = %.4f is a point of zero density" % hc)
    # on the real gap around c, H is real-analytic and H' > 0: the zero-density interval h(J)
    for x in (c - A / 8, c, c + A / 8):
        av = -G(x).real
        i2 = sum(w / (x - s_) ** 2 for s_, w in nu)
        P("    x=%9.2f: G(x)=%+.3e  H'(x)=%.6f  H(x)=%.4f" % (x, G(x).real, 1 - t * i2 / (2 * av ** 1.5),
                                                          x + t / math.sqrt(av)))
    lo_g, hi_g = None, None
    for k in range(1, 4000):
        x = -A + k * (A - At[n - 2]) / 4000
        av = -G(x).real
        i2 = sum(w / (x - s_) ** 2 for s_, w in nu)
        ok = av > 0 and 1 - t * i2 / (2 * av ** 1.5) > 0
        if ok and lo_g is None:
            lo_g = x
        if ok:
            hi_g = x
        if not ok and lo_g is not None and x > c:
            break
    P("    grid: G<0 and H'>0 on [%.2f, %.2f] (spacing %.3f); H maps it onto [%.2f, %.2f]"
      % (lo_g, hi_g, (A - At[n - 2]) / 4000, lo_g + t / math.sqrt(-G(lo_g).real),
         hi_g + t / math.sqrt(-G(hi_g).real)))
    # flanking atoms: f(lam) = largest y>0 with Im H(lam+iy) = y - t Re sqrt F(lam+iy) = 0
    for lam in (-A, -At[n - 2]):
        imH = lambda y: y - t * cmath.sqrt(1 / G(lam + 1j * y)).real
        ys = [10.0 ** (k / 20) for k in range(200, -400, -1)]
        hi_y = next(y for y in ys if imH(y) <= 0)
        lo_y, up_y = hi_y, hi_y * 10 ** (1 / 20)
        for _ in range(200):
            m = (lo_y + up_y) / 2
            if imH(m) <= 0:
                lo_y = m
            else:
                up_y = m
        wv = lam + 1j * lo_y
        hval = (wv - 1j * t * cmath.sqrt(1 / G(wv)))
        p = -G(wv).imag / math.pi
        # cross-check: run the subordination iteration at the image point
        w2 = omega(hval.real + 1e-12j, t, it=20000)
        P("    atom %8.1f: f=%.4e, h(atom)=%.6f (Im %.1e), density %.4e; omega(h) = %.4f%+.4ei"
          % (lam, lo_y, hval.real, hval.imag, p, w2.real, w2.imag))

txt = "\n".join(out) + "\n"
print(txt)
open(__file__.replace(".py", ".out"), "w", encoding="utf-8", newline="\n").write(txt)
