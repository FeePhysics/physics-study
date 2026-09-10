"""ชุดทดสอบรวม — ระบบที่ไม่เคยเจอ ผสมทุกบท · ทุกคำตอบคำนวณจริง"""
import sympy as sp

t = sp.symbols('t', real=True)
m, k, g, l, w, a = sp.symbols('m k g l omega a', positive=True)

def EL(L, q):
    return sp.simplify(sp.diff(sp.diff(L, sp.diff(q,t)), t) - sp.diff(L, q))

print("=== ① สปริงแนวดิ่ง (ระบบใหม่) L = ½mv² − ½ky² − mgy ===")
y = sp.Function('y')(t)
L1 = sp.Rational(1,2)*m*sp.diff(y,t)**2 - sp.Rational(1,2)*k*y**2 - m*g*y
print("  ∂L/∂v      =", sp.diff(L1, sp.diff(y,t)))
print("  ∂L/∂y      =", sp.diff(L1, y), "  ⇒ พิมพ์: -k*y-m*g")
sol = sp.solve(EL(L1,y), sp.diff(y,t,2))
print("  a          =", sol[0], "  ⇒ พิมพ์: -k*y/m-g")
print("  จุดสมดุล (a=0) : y =", sp.solve(sol[0], y)[0], "  ⇒ พิมพ์: -m*g/k")

print("\n=== ② ลูกปัดบนลวดพาราโบลา y = a·x²  (ระบบใหม่ · ข้อท้าทาย) ===")
x = sp.Function('x')(t)
Y = a*x**2
T = sp.Rational(1,2)*m*(sp.diff(x,t)**2 + sp.diff(Y,t)**2)
V = m*g*Y
L2 = sp.expand(T - V)
print("  ẏ = d(ax²)/dt =", sp.diff(Y,t), "  ⇒ พิมพ์: 2*a*x*v")
print("  L =", L2)
print("  ⇒ พจน์จลน์ = ½m(1+4a²x²)v²  ← สัมประสิทธิ์ขึ้นกับ x เอง")
print("  ∂L/∂x =", sp.simplify(sp.diff(L2, x)))

print("\n=== ③ อนุภาคอิสระในพิกัดเชิงขั้ว (transfer จากบทที่ 2) ===")
r, th = sp.Function('r')(t), sp.Function('theta')(t)
L3 = sp.Rational(1,2)*m*(sp.diff(r,t)**2 + r**2*sp.diff(th,t)**2)
print("  ∂L/∂θ  =", sp.diff(L3, th), " ⇒ θ เป็นพิกัดวัฏจักร")
print("  ปริมาณอนุรักษ์ =", sp.diff(L3, sp.diff(th,t)), "  ⇒ พิมพ์: m*r^2*w")

print("\n=== ④ มิติเชิงมวลของค่าคงที่ในพจน์ต่าง ๆ ([ℒ]=M⁴ · [φ]=M) ===")
for term, powr in [("λφ⁶", 6), ("gφ³", 3), ("cφ²(∂φ)²", None)]:
    if powr:
        print(f"  {term:12} : [c]·M^{powr} = M⁴ ⇒ [c] = M^{4-powr}")
print("  cφ²(∂φ)²   : [c]·M²·M⁴ = M⁴ ⇒ [c] = M^-2")

print("\n=== ⑤ สองอนุภาคเชื่อมสปริง (สมมาตรเลื่อน) ===")
x1, x2 = sp.Function('x1')(t), sp.Function('x2')(t)
L5 = sp.Rational(1,2)*m*(sp.diff(x1,t)**2 + sp.diff(x2,t)**2) - sp.Rational(1,2)*k*(x1-x2)**2
print("  เลื่อนทั้งระบบ x→x+c : L ไม่เปลี่ยน (ขึ้นกับ x1−x2 เท่านั้น)")
p = sp.diff(L5, sp.diff(x1,t)) + sp.diff(L5, sp.diff(x2,t))
print("  ปริมาณอนุรักษ์ = ", p, "  ⇒ โมเมนตัมรวม")
print("  ตรวจ d/dt = 0 :", sp.simplify(sp.diff(p,t).subs({
    sp.diff(x1,t,2): sp.solve(EL(L5,x1), sp.diff(x1,t,2))[0],
    sp.diff(x2,t,2): sp.solve(EL(L5,x2), sp.diff(x2,t,2))[0]})) == 0)

print("\n=== ⑥ สปริงที่ค่าคงที่ลดลงตามเวลา  L = ½mv² − ½k·e^(−ωt)·y² ===")
q = sp.Function('q')(t)
L6 = sp.Rational(1,2)*m*sp.diff(q,t)**2 - sp.Rational(1,2)*k*sp.exp(-w*t)*q**2
Q, Qd = sp.symbols('Q Qd')
pt = sp.diff(L6.subs(sp.diff(q,t), Qd).subs(q, Q), t).subs({Q: q, Qd: sp.diff(q,t)})
print("  ∂L/∂t =", sp.simplify(pt))
H = sp.simplify(sp.diff(q,t)*sp.diff(L6, sp.diff(q,t)) - L6)
dH = sp.simplify(sp.diff(H,t).subs(sp.diff(q,t,2), sp.solve(EL(L6,q), sp.diff(q,t,2))[0]))
print("  dH/dt =", dH)
print("  ตรงกับ −∂L/∂t :", sp.simplify(dH + pt) == 0, " ⇒ พลังงานไม่อนุรักษ์ (ลดลง — สปริงอ่อนลงเรื่อย ๆ)")

print("\n=== ⑦ แปรผัน L ใหม่: อนุภาคอิสระ m=2, y₀=3t, η=t(1−t), ช่วง [0,1] ===")
eps = sp.symbols('epsilon', real=True)
y0, eta = 3*t, t*(1-t)
Ye = y0 + eps*eta
S = sp.expand(sp.integrate(sp.Rational(1,2)*2*sp.diff(Ye,t)**2, (t,0,1)))
P = sp.Poly(S, eps)
print("  ẏ₀ =", sp.diff(y0,t), " η̇ =", sp.expand(sp.diff(eta,t)))
print("  S(ε) =", S)
print("  S(0) =", P.coeff_monomial(1), " · สัมประสิทธิ์ ε =", P.coeff_monomial(eps),
      " · ของ ε² =", P.coeff_monomial(eps**2))
print("  ⇒ ∫ẏ₀η̇dt =", sp.integrate(sp.diff(y0,t)*sp.diff(eta,t), (t,0,1)), " (ไม่มีพจน์ V)")
