"""เขียนปริมาณอนุรักษ์ออกมาเอง — รวมกรณีที่ 'ไม่มีตัวแปรไหนหายไป' """
import sympy as sp
t = sp.symbols('t', real=True)
m, k, g, l, M = sp.symbols('m k g l M', positive=True)

def EL(L, q):
    return sp.simplify(sp.diff(sp.diff(L, sp.diff(q,t)), t) - sp.diff(L, q))

print("=== กรณี A · ตัวแปรหายไปจาก L ตรง ๆ (พิกัดวัฏจักร) ===")
for lbl, L, q in [
    ("เชิงขั้ว θ", sp.Rational(1,2)*m*(sp.diff(sp.Function('r')(t),t)**2
        + sp.Function('r')(t)**2*sp.diff(sp.Function('theta')(t),t)**2), sp.Function('theta')(t)),
]:
    print(f"  {lbl}: ∂L/∂q = {sp.diff(L,q)}  ⇒  p = {sp.diff(L, sp.diff(q,t))}")

th, ph = sp.Function('theta')(t), sp.Function('phi')(t)
Lsph = sp.Rational(1,2)*m*l**2*(sp.diff(th,t)**2 + sp.sin(th)**2*sp.diff(ph,t)**2) + m*g*l*sp.cos(th)
print(f"  ลูกตุ้มทรงกลม φ: ∂L/∂φ = {sp.diff(Lsph,ph)}  ⇒  p = {sp.diff(Lsph, sp.diff(ph,t))}")

x, y, z = [sp.Function(c)(t) for c in 'xyz']
Lz = sp.Rational(1,2)*m*(sp.diff(x,t)**2+sp.diff(y,t)**2+sp.diff(z,t)**2) - sp.Function('V')(z)
print(f"  V ขึ้นกับ z เท่านั้น: ∂L/∂x = {sp.diff(Lz,x)}  ⇒  p_x = {sp.diff(Lz, sp.diff(x,t))}"
      f"  ·  p_y = {sp.diff(Lz, sp.diff(y,t))}")

print("\n=== กรณี B · ไม่มีตัวแปรไหนหายไปเลย แต่ยังมีสมมาตร ===")
x1, x2 = sp.Function('x1')(t), sp.Function('x2')(t)
L2 = sp.Rational(1,2)*m*(sp.diff(x1,t)**2+sp.diff(x2,t)**2) - sp.Rational(1,2)*k*(x1-x2)**2
print(f"  ∂L/∂x1 = {sp.simplify(sp.diff(L2,x1))}   ← ไม่เป็นศูนย์!")
print(f"  ∂L/∂x2 = {sp.simplify(sp.diff(L2,x2))}   ← ไม่เป็นศูนย์!")
print("  ⇒ ไม่มีพิกัดวัฏจักรเลย · แต่ ∂L/∂x1 + ∂L/∂x2 =",
      sp.simplify(sp.diff(L2,x1) + sp.diff(L2,x2)), "← ผลรวมเป็นศูนย์")
P = sp.diff(L2, sp.diff(x1,t)) + sp.diff(L2, sp.diff(x2,t))
print(f"  ⇒ ปริมาณอนุรักษ์ = ผลรวมของโมเมนตัม = {P}")
a1 = sp.solve(EL(L2,x1), sp.diff(x1,t,2))[0]; a2 = sp.solve(EL(L2,x2), sp.diff(x2,t,2))[0]
print(f"  ตรวจ d/dt = 0 :", sp.simplify(sp.diff(P,t).subs({sp.diff(x1,t,2):a1, sp.diff(x2,t,2):a2})) == 0)
print(f"  แต่ละตัวไม่คงที่: d(m·v1)/dt = {sp.simplify(m*a1)}  ← ไม่เป็นศูนย์")

print("\n=== กรณี C · มวลต่างกัน (ยังใช้กฎเดิม) ===")
L3 = sp.Rational(1,2)*m*sp.diff(x1,t)**2 + sp.Rational(1,2)*M*sp.diff(x2,t)**2 \
     - sp.Rational(1,2)*k*(x1-x2)**2
P3 = sp.diff(L3, sp.diff(x1,t)) + sp.diff(L3, sp.diff(x2,t))
print(f"  ปริมาณอนุรักษ์ = {P3}  ⇒ พิมพ์: m*v1+M*v2")
b1 = sp.solve(EL(L3,x1), sp.diff(x1,t,2))[0]; b2 = sp.solve(EL(L3,x2), sp.diff(x2,t,2))[0]
print("  ตรวจ d/dt = 0 :", sp.simplify(sp.diff(P3,t).subs({sp.diff(x1,t,2):b1, sp.diff(x2,t,2):b2})) == 0)

print("\n=== กรณี D · แปรผันเมื่อไม่มี V และ m ไม่ใช่ 1 (ข้อ 5 ที่พลาด) ===")
eps = sp.symbols('epsilon', real=True)
for mm, y0, eta, rng in [(2, 3*t, t*(1-t), (0,1)), (5, 2*t, t*(1-t), (0,1)), (3, 4*t, t*(2-t), (0,2))]:
    Ye = y0 + eps*eta
    S = sp.expand(sp.integrate(sp.Rational(1,2)*mm*sp.diff(Ye,t)**2, (t,)+rng))
    P = sp.Poly(S, eps)
    lin = sp.integrate(mm*sp.diff(y0,t)*sp.diff(eta,t), (t,)+rng)
    print(f"  m={mm}, y₀={y0}, ช่วง{rng}: สัมประสิทธิ์ ε = {P.coeff_monomial(eps)}"
          f"  = ∫m·ẏ₀·η̇dt = {lin}  · ∫η̇dt = {sp.integrate(sp.diff(eta,t),(t,)+rng)}")
print("  ⇒ ẏ₀ คงที่ ⇒ ดึงออกนอกอินทิกรัลได้ ⇒ เหลือ m·ẏ₀·∫η̇dt = m·ẏ₀·[η] = 0 เสมอ")
