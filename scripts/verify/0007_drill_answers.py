"""เฉลยทุกข้อของชุดฝึกคำตอบเปิด (บทที่ 4) — คำนวณจริง ไม่ได้พิมพ์มือ"""
import sympy as sp

t, eps = sp.symbols('t epsilon', real=True)
m, k, g = sp.symbols('m k g', positive=True)

def EL(L, q):
    return sp.simplify(sp.diff(sp.diff(L, sp.diff(q,t)), t) - sp.diff(L, q))

print("ข้อ 1 — สปริง L = ½mẏ² − ½ky² ⇒ ÿ = ?")
y = sp.Function('y')(t)
sol = sp.solve(EL(sp.Rational(1,2)*m*sp.diff(y,t)**2 - sp.Rational(1,2)*k*y**2, y), sp.diff(y,t,2))
print("   ÿ =", sol[0], "   ⇒ พิมพ์: -k*y/m หรือ -(k/m)y")

print("\nข้อ 2 — L = ½mẏ² + bẏ  (b คงที่) ⇒ ÿ = ?")
b = sp.symbols('b', positive=True)
sol2 = sp.solve(EL(sp.Rational(1,2)*m*sp.diff(y,t)**2 + b*sp.diff(y,t), y), sp.diff(y,t,2))
print("   ÿ =", sol2[0], "   ⇒ พจน์ bẏ เป็นอนุพันธ์แท้ ไม่มีผลต่อสมการเลย")

print("\nข้อ 3 — พจน์ขอบ [m·ẏ·η] เมื่อ y = 5t(2−t), η = t²(2−t), ช่วง 0→2")
Y = 5*t*(2-t); eta = t**2*(2-t)
bt = m*sp.diff(Y,t)*eta
print("   ที่ t=2:", sp.simplify(bt.subs(t,2)), " · ที่ t=0:", sp.simplify(bt.subs(t,0)))
print("   ⇒ พจน์ขอบ =", sp.simplify(bt.subs(t,2) - bt.subs(t,0)), "  (η(0)=η(2)=0 ⇒ หายทั้งคู่)")

print("\nข้อ 4 — η = t(2−t)² ปลายเป็นศูนย์ไหม")
e4 = t*(2-t)**2
print("   η(0) =", e4.subs(t,0), " · η(2) =", e4.subs(t,2), " ⇒ ใช้ได้")

print("\nข้อ 5 — ∫₀² η·f dt เมื่อ η = t(2−t), f = t−1")
I = sp.integrate(t*(2-t)*(t-1), (t,0,2))
print("   =", I, "  ⇒ ศูนย์! ทั้งที่ f ไม่ใช่ศูนย์ (η สมมาตรรอบ t=1 · f ปฏิสมมาตร)")

print("\nข้อ 6 — ∫₀² η·f dt เมื่อ η = t(2−t)(t−1), f = t−1")
I6 = sp.integrate(t*(2-t)*(t-1)*(t-1), (t,0,2))
print("   =", I6, "=", float(I6), "  ⇒ ไม่เป็นศูนย์ ⇒ จับได้ว่า f ≠ 0")

print("\nข้อ 7 — dS/dε|₀ ของ L = ½mẏ² − mgy รอบเส้นทางจริง y = 5t(2−t), g=10, m=1")
Yv = 5*t*(2-t) + eps*t*(2-t)
Lv = sp.Rational(1,2)*sp.diff(Yv,t)**2 - 10*Yv
Sv = sp.integrate(Lv, (t,0,2))
print("   S(ε) =", sp.expand(Sv))
print("   dS/dε|₀ =", sp.diff(Sv,eps).subs(eps,0), "  ⇒ ศูนย์เพราะเป็นเส้นทางจริง")

print("\nข้อ 8 — เส้นทางที่ไม่ใช่เส้นทางจริง: y = 4t(2−t) (สูงสุด 4 ไม่ใช่ 5)")
Yw = 4*t*(2-t) + eps*t*(2-t)
Sw = sp.integrate(sp.Rational(1,2)*sp.diff(Yw,t)**2 - 10*Yw, (t,0,2))
print("   dS/dε|₀ =", sp.diff(Sw,eps).subs(eps,0), "  ⇒ ไม่เป็นศูนย์ ⇒ ไม่ใช่เส้นทางจริง")

print("\nข้อ 9 — มิติเชิงมวลของ g ในพจน์ g·φ³ (ℒ มีมิติ M⁴, [φ]=M)")
print("   [g]·M³ = M⁴ ⇒ [g] = M¹  ⇒ ตอบ 1")

print("\nข้อ 10 — ℒ = ½(∂φ)² − ½m²φ² − (λ/4!)φ⁴ · EL ให้อะไรเพิ่มจากไคลน์–กอร์ดอน")
tt, xx = sp.symbols('t x', real=True)
lam = sp.symbols('lambda', positive=True)
ph = sp.Function('phi')(tt, xx)
Lg = (sp.Rational(1,2)*(sp.diff(ph,tt)**2 - sp.diff(ph,xx)**2)
      - sp.Rational(1,2)*m**2*ph**2 - lam/24*ph**4)
el = (sp.diff(sp.diff(Lg, sp.diff(ph,tt)), tt) + sp.diff(sp.diff(Lg, sp.diff(ph,xx)), xx)
      - sp.diff(Lg, ph))
print("   ", sp.simplify(el), "= 0")
print("   ⇒ พจน์ใหม่คือ (λ/6)φ³")
