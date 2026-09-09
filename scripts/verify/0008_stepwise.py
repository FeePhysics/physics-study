"""ชุดฝึกไล่ขั้น — เฉลยของ *ทุกขั้นย่อย* ไม่ใช่แค่คำตอบสุดท้าย"""
import sympy as sp

t, eps = sp.symbols('t epsilon', real=True)
m, k = sp.symbols('m k', positive=True)

print("=== โจทย์ A · สปริง L = ½mv² − ½ky²  (v = ẏ, a = ÿ) ===")
y = sp.Function('y')(t)
L = sp.Rational(1,2)*m*sp.diff(y,t)**2 - sp.Rational(1,2)*k*y**2
p = sp.diff(L, sp.diff(y,t))
print("  ขั้น 1 · ∂L/∂v     =", p,            "  ⇒ พิมพ์: m*v")
print("  ขั้น 2 · d/dt ของมัน =", sp.diff(p,t), "  ⇒ พิมพ์: m*a")
print("  ขั้น 3 · ∂L/∂y     =", sp.diff(L, y), "  ⇒ พิมพ์: -k*y")
sol = sp.solve(sp.diff(p,t) - sp.diff(L,y), sp.diff(y,t,2))
print("  ขั้น 4 · a          =", sol[0],       "  ⇒ พิมพ์: -k*y/m")

print("\n=== โจทย์ B · พจน์ขอบ  y = 5t(2−t), m = 1, ช่วง [0,2] ===")
Y = 5*t*(2-t)
for name, e in [("t²(2−t)", t**2*(2-t)), ("t", t)]:
    b = sp.diff(Y,t)*e
    print(f"  η = {name:8}  η(0)={e.subs(t,0)}  η(2)={e.subs(t,2)}"
          f"  ⇒ [ẏη] = {sp.simplify(b.subs(t,2) - b.subs(t,0))}")

print("\n=== โจทย์ C · แปรผันทีละขั้น  m=1, g=10 ===")
for lab, base in [("y = 4t(2−t) (เส้นทางผิด)", 4*t*(2-t)),
                  ("y = 5t(2−t) (เส้นทางจริง)", 5*t*(2-t))]:
    Ye = base + eps*t*(2-t)
    S = sp.expand(sp.integrate(sp.Rational(1,2)*sp.diff(Ye,t)**2 - 10*Ye, (t,0,2)))
    poly = sp.Poly(S, eps)
    A, B, C = poly.coeff_monomial(1), poly.coeff_monomial(eps), poly.coeff_monomial(eps**2)
    print(f"  {lab}")
    print(f"    S(ε) = {S}")
    print(f"    ขั้น 1 · S(0) = {A} = {float(A):.4f}")
    print(f"    ขั้น 2 · สัมประสิทธิ์ของ ε = {B} = {float(B):.4f}")
    print(f"    ขั้น 3 · dS/dε|₀ = {sp.diff(S,eps).subs(eps,0)}")
    print(f"    ขั้น 4 · สัมประสิทธิ์ของ ε² = {C}  (บวก = ค่าต่ำสุด)")

print("\n=== โจทย์ D · ทำไม bẏ ถึงไม่มีผล (ต่อยอดข้อที่ทำถูก) ===")
b = sp.symbols('b', positive=True)
L2 = sp.Rational(1,2)*m*sp.diff(y,t)**2 + b*sp.diff(y,t)
print("  ∂L/∂v =", sp.diff(L2, sp.diff(y,t)), " ⇒ พิมพ์: m*v+b")
print("  d/dt  =", sp.diff(sp.diff(L2, sp.diff(y,t)), t), " ⇒ b เป็นค่าคงที่ อนุพันธ์เป็นศูนย์")
print("  ∂L/∂y =", sp.diff(L2, y))
print("  ⇒ a =", sp.solve(sp.diff(sp.diff(L2,sp.diff(y,t)),t) - sp.diff(L2,y), sp.diff(y,t,2))[0])
