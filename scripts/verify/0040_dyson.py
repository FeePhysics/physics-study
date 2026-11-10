"""บทที่ 26 · Dyson series — ตรวจทุกตัวเลขในบท"""
import sympy as sp, numpy as np
from scipy.linalg import expm
t, t1, t2, t3, E, v, lam, phi, phidot = sp.symbols('t t1 t2 t3 E v lambda phi phidot', real=True)
ok = 0
def check(name, cond):
    global ok
    assert cond, name; ok += 1; print("✅", name)

# ① e^{-iEt} แก้ i dψ/dt = E ψ
psi = sp.exp(-sp.I*E*t)
check("Schrödinger", sp.simplify(sp.I*sp.diff(psi, t) - E*psi) == 0)
check("d/dt e^{-3it} = -3i e^{-3it}", sp.simplify(sp.diff(sp.exp(-3*sp.I*t), t)/sp.exp(-3*sp.I*t)) == -3*sp.I)
check("i*(-3i)=3", sp.I*(-3*sp.I) == 3)
# ② กาง e^{-ivt}
ser = sp.series(sp.exp(-sp.I*v*t), v, 0, 4).removeO()
check("order1 = -i v t", sp.expand(ser.coeff(v, 1)) == -sp.I*t)
check("order2 = -1/2 v^2 t^2", sp.expand(ser.coeff(v, 2)) == -t**2/2)
check("(-i)^2 = -1", (-sp.I)**2 == -1)
check("(-i)^3/3! = i/6", sp.simplify((-sp.I)**3/6 - sp.I/6) == 0)
# ③ สามเหลี่ยม = ครึ่งสี่เหลี่ยม (สำหรับ f ใด ๆ ที่สลับกันได้)
f = sp.Function('f')
for g in [sp.cos(t), t**2 + 1, sp.exp(t)]:
    gg = lambda s: g.subs(t, s)
    tri = sp.integrate(sp.integrate(gg(t1)*gg(t2), (t2, 0, t1)), (t1, 0, t))
    sq = sp.integrate(gg(t1), (t1, 0, t))**2
    check(f"triangle=square/2 for {g}", sp.simplify(tri - sq/2) == 0)
check("area tri t=1", sp.integrate(sp.integrate(1, (t2, 0, t1)), (t1, 0, 1)) == sp.Rational(1, 2))
check("vol tetra t=1", sp.integrate(sp.integrate(sp.integrate(1, (t3, 0, t2)), (t2, 0, t1)), (t1, 0, 1)) == sp.Rational(1, 6))
# Dyson กับเวลาแปรผัน v(t) ที่สลับได้: exp(-i∫v) = 1 - i∫v + (-i)^2 ∫∫_{t1>t2} v v + ...
vt = sp.cos(t)
exact = sp.series(sp.exp(-sp.I*sp.Symbol('e')*sp.integrate(vt.subs(t, t1), (t1, 0, t))), sp.Symbol('e'), 0, 3).removeO()
d2 = sp.integrate(sp.integrate(sp.cos(t1)*sp.cos(t2), (t2, 0, t1)), (t1, 0, t))
check("Dyson 2nd order = series", sp.simplify(exact.coeff(sp.Symbol('e'), 2) - (-sp.I)**2*d2) == 0)
# เมทริกซ์ที่สลับไม่ได้: U_I ≈ 1 - i∫V_I dt + O(ε²)
rng = np.random.default_rng(0)
H0 = np.diag([0.0, 1.3]); V = np.array([[0.2, 0.7], [0.7, -0.4]])
eps, T = 1e-3, 2.0
U = expm(-1j*(H0 + eps*V)*T); UI = expm(1j*H0*T) @ U
ts = np.linspace(0, T, 4001)
VI = [expm(1j*H0*s) @ V @ expm(-1j*H0*s) for s in ts]
first = np.trapezoid(np.array(VI), ts, axis=0)
err = np.abs(UI - (np.eye(2) - 1j*eps*first)).max()
check(f"matrix 1st order err {err:.1e} ~ eps^2", err < 10*eps**2)
# ④ H_int = -L_int (ไม่มี phidot) · Legendre
L = sp.Rational(1, 2)*phidot**2 - lam/24*phi**4
pi = sp.diff(L, phidot); H = sp.expand(pi*phidot - L)
check("H_int = +λ/4! φ^4", sp.expand(H - sp.Rational(1, 2)*phidot**2) == lam/24*phi**4)
check("-i λ/4! * 4! = -iλ", sp.simplify(-sp.I*lam/24*24 + sp.I*lam) == 0)
check("2nd order prefactor (-i)^2/2! = -1/2", (-sp.I)**2/2 == sp.Rational(-1, 2))
check("swap 2 vertices cancels 1/2!", sp.Rational(1, 2)*2 == 1)

# ฉบับเขียนใหม่ (รอบสอง): ตัวเลขอุ่นเครื่อง + ประกอบกับบทที่ 25
check("e^0.1 ≈ 1.1052", abs(float(sp.exp(sp.Rational(1,10))) - 1.10517) < 1e-5)
check("0.1^2/2 = 0.005", sp.Rational(1,10)**2/2 == sp.Rational(1,200))
check("e^{-i(5+0.1)t} = e^{-5it} e^{-0.1it}", sp.simplify(sp.exp(-sp.I*(5+sp.Rational(1,10))*t) - sp.exp(-5*sp.I*t)*sp.exp(-sp.I*t/10)) == 0)
check("1st order of e^{-0.1it}: -0.1i t", sp.series(sp.exp(-sp.I*t/10), t, 0, 2).removeO().coeff(t, 1) == -sp.I/10)
check("A: -i * 1/2 = -i/2", -sp.I*sp.Rational(12, 24) == -sp.I/2)
check("X: -i * 24/24 = -i", -sp.I*sp.Rational(24, 24) == -sp.I)
check("B: -i * 3/24 = -i/8", -sp.I*sp.Rational(3, 24) == -sp.I/8)
print(f"{ok}/{ok}")
