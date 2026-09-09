"""ชุดฝึก 'ประกอบนิพจน์ก่อนอินทิเกรต' — ไต่จากง่ายมากไปถึงของจริงที่พลาด"""
import sympy as sp
t, eps = sp.symbols('t epsilon', real=True)

def show(lbl, y, m, g, a, b, eta=None):
    yd = sp.diff(y, t)
    T = sp.Rational(1,2)*m*yd**2
    V = m*g*y
    print(f"\n=== {lbl} ===")
    print(f"  y(t)     = {sp.expand(y)}")
    print(f"  ẏ(t)     = {sp.expand(yd)}")
    print(f"  T = ½m ẏ² = {sp.expand(T)}")
    if g: print(f"  V = mgy   = {sp.expand(V)}")
    print(f"  L = T − V = {sp.expand(T - V)}")
    IT = sp.integrate(T, (t,a,b)); IV = sp.integrate(V, (t,a,b))
    print(f"  ∫T dt = {IT}   ∫V dt = {IV}   ⇒ S = {sp.simplify(IT-IV)}")

show("A · y = t², m = 2, g = 0, ช่วง [0,1]", t**2, 2, 0, 0, 1)
show("B · y = t², m = 2, g = 10, ช่วง [0,1]", t**2, 2, 10, 0, 1)
show("C · y = 4t(2−t), m = 1, g = 10, ช่วง [0,2]", 4*t*(2-t), 1, 10, 0, 2)

print("\n=== D · ใส่การบิด: y = 4t(2−t) + ε·t(2−t) ===")
y0, eta = 4*t*(2-t), t*(2-t)
Y = y0 + eps*eta
yd = sp.expand(sp.diff(Y,t))
print(f"  ẏ_ε      = {yd}")
sq = sp.expand(yd**2)
print(f"  (ẏ_ε)²   = {sq}")
print(f"  พจน์ที่มี ε ตัวเดียว ใน (ẏ_ε)² : {sp.Poly(sq, eps).coeff_monomial(eps)}")
print(f"       = 2·ẏ₀·η̇ = {sp.expand(2*sp.diff(y0,t)*sp.diff(eta,t))}  ← ตัวคูณไขว้")
S = sp.expand(sp.integrate(sp.Rational(1,2)*yd**2 - 10*Y, (t,0,2)))
P = sp.Poly(S, eps)
print(f"\n  S(ε) = {S}")
print(f"  ขั้น 1 · S(0)              = {P.coeff_monomial(1)}")
print(f"  ขั้น 2 · สัมประสิทธิ์ของ ε   = {P.coeff_monomial(eps)}")
print(f"           แยกเป็น ∫ẏ₀η̇dt = {sp.integrate(sp.diff(y0,t)*sp.diff(eta,t),(t,0,2))}"
      f"  ลบ  ∫10η dt = {sp.integrate(10*eta,(t,0,2))}")
print(f"  ขั้น 4 · สัมประสิทธิ์ของ ε²  = {P.coeff_monomial(eps**2)}")
