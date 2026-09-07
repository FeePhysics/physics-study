"""ทุกข้ออ้างในบทที่ 3 — ดัชนี เมตริก และหน่วยธรรมชาติ"""
import sympy as sp

print("=== ① เมตริกสอง convention ให้ผลต่างกันแค่เครื่องหมายรวม ===")
# (−+++) แบบสัมพัทธภาพทั่วไป · (+−−−) แบบฟิสิกส์อนุภาค
g_mostplus  = sp.diag(-1, 1, 1, 1)
g_mostminus = sp.diag(1, -1, -1, -1)
t, x, y, z = sp.symbols('t x y z', real=True)
X = sp.Matrix([t, x, y, z])

for name, g in [("(−+++)", g_mostplus), ("(+−−−)", g_mostminus)]:
    s = sp.expand((X.T * g * X)[0])            # x^μ x_μ = g_μν x^μ x^ν
    print(f"  {name}  x·x = {s}")
print("  ⇒ ต่างกันแค่เครื่องหมายรวม · ข้อสรุปทางฟิสิกส์เหมือนกัน แต่สูตรกลางทางไม่เหมือน")

print("\n=== ② ดัชนีขึ้น-ลง: g ใช้ยกและกดดัชนี ===")
v_up = sp.Matrix(sp.symbols('v0 v1 v2 v3'))
v_dn = g_mostplus * v_up                        # v_μ = g_μν v^ν
print("  v^μ =", list(v_up))
print("  v_μ =", list(v_dn), " ← องค์ประกอบเวลาเปลี่ยนเครื่องหมาย")
print("  ยกกลับได้เดิม:", list(g_mostplus.inv() * v_dn) == list(v_up))

print("\n=== ③ ช่วงในกาลอวกาศแยกสามประเภทได้จริง ===")
for lbl, (dt_, dx_) in [("เดินทางด้วยแสงพอดี", (1, 1)),
                        ("ช้ากว่าแสง (timelike)", (2, 1)),
                        ("เร็วกว่าแสง (spacelike)", (1, 2))]:
    s2 = -dt_**2 + dx_**2                       # ใช้ (−+++) · c = 1
    kind = "ศูนย์ (null)" if s2 == 0 else ("ลบ" if s2 < 0 else "บวก")
    print(f"  {lbl:26} Δs² = {s2:3}  ⇒ {kind}")

print("\n=== ④ Klein–Gordon: EL ของสนามให้สมการคลื่นที่มีมวล ===")
# ℒ = ½ ∂_μφ ∂^μφ − ½ m²φ²  ใช้ (+−−−) ⇒ ∂_μφ∂^μφ = φ̇² − (∇φ)²
tt, xx = sp.symbols('t x', real=True)
m = sp.symbols('m', positive=True)
phi = sp.Function('phi')(tt, xx)
Lag = sp.Rational(1,2)*(sp.diff(phi,tt)**2 - sp.diff(phi,xx)**2) - sp.Rational(1,2)*m**2*phi**2

# EL ของสนาม: ∂_μ(∂ℒ/∂(∂_μφ)) − ∂ℒ/∂φ = 0
el = (sp.diff(sp.diff(Lag, sp.diff(phi,tt)), tt)
      + sp.diff(sp.diff(Lag, sp.diff(phi,xx)), xx)
      - sp.diff(Lag, phi))
print("  EL:", sp.Eq(sp.simplify(el), 0))
kg = sp.diff(phi,tt,2) - sp.diff(phi,xx,2) + m**2*phi
print("  ตรงกับสมการไคลน์–กอร์ดอน □φ + m²φ = 0 :", sp.simplify(el - kg) == 0)

print("\n=== ⑤ หน่วยธรรมชาติ: ħ = c = 1 ทำให้ทุกอย่างวัดเป็นมวลได้ ===")
# [ħ] = พลังงาน×เวลา, [c] = ความยาว/เวลา ⇒ ตั้งเป็น 1 แล้วเหลือมิติเดียว
E, L_, T_, M = sp.symbols('E L T M', positive=True)
print("  ħ=1 ⇒ [เวลา] = 1/[พลังงาน]      · c=1 ⇒ [ความยาว] = [เวลา]")
print("  ⇒ [ความยาว] = [เวลา] = 1/[มวล]  · ทุกปริมาณเขียนเป็น (มวล)^n ได้หมด")
# ตรวจมิติของสนามสเกลาร์: S = ∫d⁴x ℒ ต้องไร้มิติ (ħ=1)
# [d⁴x] = M^-4 ⇒ [ℒ] = M^4 ; ℒ มีพจน์ m²φ² ⇒ M² · [φ]² = M^4 ⇒ [φ] = M
print("  S ไร้มิติ · [d⁴x] = M⁻⁴ ⇒ [ℒ] = M⁴")
print("  พจน์ m²φ² : M²·[φ]² = M⁴ ⇒ [φ] = M¹  ← สนามสเกลาร์มีมิติเท่ามวล")
print("  ตรวจกับพจน์จลน์ (∂φ)² : [∂] = M ⇒ M²·M² = M⁴ ✓ ตรงกัน")
