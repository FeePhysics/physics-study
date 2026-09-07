"""บทที่ 4 — พิสูจน์ EL จากการแปรผัน · ทุกขั้นตรวจด้วย sympy"""
import sympy as sp

t, eps = sp.symbols('t epsilon', real=True)
m, g = sp.symbols('m g', positive=True)
t1, t2 = 0, 2

y  = sp.Function('y')(t)          # เส้นทางที่ยังไม่รู้
eta = sp.Function('eta')(t)       # การบิด — ยังไม่กำหนดรูป

print("=== ① แปรผัน y → y + ε·η แล้วกระจาย S(ε) ตามกำลังของ ε ===")
Y = y + eps*eta
L = sp.Rational(1,2)*m*sp.diff(Y,t)**2 - m*g*Y
S = sp.integrate(L, (t, t1, t2))

dS = sp.diff(L, eps).subs(eps, 0)          # ตัวถูกอินทิเกรตของ dS/dε ที่ ε=0
dS = sp.expand(dS)
print("  dS/dε|₀ =  ∫ (", dS, ") dt")

print("\n=== ② อนุพันธ์ทีละส่วน: ย้ายอนุพันธ์จาก η̇ ไปที่ ẏ ===")
# พจน์ m·ẏ·η̇ = d/dt(m·ẏ·η) − m·ÿ·η
lhs = m*sp.diff(y,t)*sp.diff(eta,t)
rhs = sp.diff(m*sp.diff(y,t)*eta, t) - m*sp.diff(y,t,2)*eta
print("  m·ẏ·η̇ = d/dt(m·ẏ·η) − m·ÿ·η  :", sp.simplify(lhs - rhs) == 0)

boundary = (m*sp.diff(y,t)*eta)            # พจน์ขอบ [m·ẏ·η]
bulk = sp.expand(-m*sp.diff(y,t,2)*eta - m*g*eta)
print("  ⇒ dS/dε|₀ = [m·ẏ·η] จาก 0 ถึง 2  +  ∫ η·(", sp.factor(bulk/eta), ") dt")

print("\n=== ③ พจน์ขอบหายไปก็ต่อเมื่อ η ที่ปลายทั้งสองเป็นศูนย์ ===")
for name, e in [("η = sin(πt/2)  ← ปลายเป็นศูนย์ทั้งคู่", sp.sin(sp.pi*t/2)),
                ("η = t          ← ปลายขวาไม่เป็นศูนย์", t)]:
    b = (m*sp.diff(y,t)*e)
    val = sp.simplify(b.subs(t, t2) - b.subs(t, t1))
    print(f"  {name}:  [m·ẏ·η] = {val}")

print("\n=== ④ บทตั้งหลักของแคลคูลัสการแปรผัน ===")
# ถ้า ∫ η(t)·f(t) dt = 0 สำหรับ *ทุก* η ที่ปลายเป็นศูนย์ ⇒ f ≡ 0
# แสดงด้วยการเลือก η ที่ชี้เป้า: ถ้า f ไม่เป็นศูนย์ที่ไหน เลือก η ให้อินทิกรัลไม่เป็นศูนย์ได้
f = sp.Function('f')
print("  ลองสมมติ f(t) = t − 1 (ไม่ใช่ศูนย์) แล้วหา η ที่ทำให้อินทิกรัลไม่เป็นศูนย์:")
ftest = t - 1
for e in [sp.sin(sp.pi*t/2), (t-1)*sp.sin(sp.pi*t/2)**2]:
    I = sp.integrate(e*ftest, (t, t1, t2))
    print(f"    η = {e}  ⇒  ∫η·f = {sp.simplify(I)} = {float(I):+.4f}")
print("  ⇒ มี η ที่ทำให้อินทิกรัลไม่เป็นศูนย์ ⇒ f ที่ไม่เป็นศูนย์ใช้ไม่ได้ ⇒ f ต้องเป็นศูนย์")

print("\n=== ⑤ สรุป: วงเล็บที่เหลือคือสมการออยเลอร์–ลากรองจ์ ===")
print("  −m·ÿ − m·g = 0  ⇒  ÿ = −g   ← ตรงกับบทที่ 1")
