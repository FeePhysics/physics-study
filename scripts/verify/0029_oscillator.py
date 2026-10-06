"""ตรวจทุกข้ออ้างของบทที่ 17 — สนาม = สปริงจำนวนมาก · ควอนตัม = ขั้นบันได · หนึ่งขั้น = หนึ่งอนุภาค

convention เดียวกับบทที่ 0016: signature (-,+) · 1+1 มิติ · hbar = c = 1
L = 1/2 phi_t^2 - 1/2 phi_x^2 - 1/2 m^2 phi^2  ⇒  phi_tt - phi_xx + m^2 phi = 0
สปริงควอนตัม: H psi = -1/2 psi'' + 1/2 w^2 x^2 psi
"""
import sympy as sp

t, x = sp.symbols('t x', real=True)
k, m, w = sp.symbols('k m w', positive=True)
ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


# ── ① โหมดเดียวของสนาม = สปริงหนึ่งตัว ─────────────────────
q = sp.Function('q')(t)
phi = q*sp.cos(k*x)
L = sp.Rational(1, 2)*sp.diff(phi, t)**2 - sp.Rational(1, 2)*sp.diff(phi, x)**2 - sp.Rational(1, 2)*m**2*phi**2
# สมการ EL ของสนามจาก L จริง (ไม่ได้พิมพ์สมการคลื่นเอง)
pt_, px_, f_ = sp.symbols('pt px f')
Ls = sp.Rational(1, 2)*pt_**2 - sp.Rational(1, 2)*px_**2 - sp.Rational(1, 2)*m**2*f_**2
ff = sp.Function('f')(t, x)
sub = {pt_: sp.diff(ff, t), px_: sp.diff(ff, x), f_: ff}
EL = sp.diff(sp.diff(Ls, pt_).subs(sub), t) + sp.diff(sp.diff(Ls, px_).subs(sub), x) - sp.diff(Ls, f_).subs(sub)
chk("EL จาก L = phi_tt - phi_xx + m^2 phi", EL, sp.diff(ff, t, 2) - sp.diff(ff, x, 2) + m**2*ff)
lhs = (sp.diff(phi, t, 2) - sp.diff(phi, x, 2) + m**2*phi)/sp.cos(k*x)
chk("phi = q cos(kx) ⇒ q'' + (k^2 + m^2) q = 0  [สมการ (1)]", lhs, sp.diff(q, t, 2) + (k**2 + m**2)*q)
# ผลเฉลยคือการแกว่ง cos(w t) ด้วย w^2 = k^2 + m^2
wk = sp.sqrt(k**2 + m**2)
chk("q = cos(w t) แก้สมการ (1) เมื่อ w^2 = k^2 + m^2",
    sp.diff(sp.cos(wk*t), t, 2) + (k**2 + m**2)*sp.cos(wk*t), 0)
# Lagrangian ของโหมด (อินทิเกรตหนึ่งคาบ) = ค่าคงที่ × Lagrangian สปริง
Lmode = sp.integrate(L, (x, 0, 2*sp.pi/k))
qd, qq = sp.symbols('qd qq')
Lmode = Lmode.subs(sp.diff(q, t), qd).subs(q, qq)
chk("∫ L dx หนึ่งคาบ = (pi/k)(1/2 q'^2 - 1/2 w^2 q^2) — Lagrangian สปริง",
    Lmode, (sp.pi/k)*(sp.Rational(1, 2)*qd**2 - sp.Rational(1, 2)*(k**2 + m**2)*qq**2))
# ตัวควบคุมฝั่งลบ: ถ้าลืม -phi_xx จะได้ w^2 = m^2 (ผิด)
bad = (sp.diff(phi, t, 2) + m**2*phi)/sp.cos(k*x)
ok.append(("ตัวควบคุมฝั่งลบ: ลืม phi_xx ⇒ ไม่ได้ k^2 + m^2",
           sp.simplify(bad - (sp.diff(q, t, 2) + (k**2 + m**2)*q)) != 0, bad, "≠"))

# ── ② สปริงควอนตัม: ขั้นบันได E_n = w (n + 1/2) ─────────────
X = sp.symbols('X', real=True)


def H(psi, ww):
    return -sp.Rational(1, 2)*sp.diff(psi, X, 2) + sp.Rational(1, 2)*ww**2*X**2*psi


for n in range(6):
    psi = sp.hermite(n, sp.sqrt(w)*X)*sp.exp(-w*X**2/2)
    chk(f"H psi_{n} = w ({n} + 1/2) psi_{n}  [สมการ (2)]", sp.simplify(H(psi, w)/psi), w*(n + sp.Rational(1, 2)))
psi_bad = sp.exp(-X**2)
ok.append(("ตัวควบคุมฝั่งลบ: e^{-x^2} ไม่ใช่สถานะพลังงานแน่นอนเมื่อ w = 1 (H psi / psi ขึ้นกับ x)",
           sp.diff(sp.simplify(H(psi_bad, 1)/psi_bad), X) != 0, sp.simplify(H(psi_bad, 1)/psi_bad), "ขึ้นกับ x"))

# ข้อ 3: psi0 = e^{-x^2/2}, w = 1
p0 = sp.exp(-X**2/2)
chk("psi0' = -x psi0", sp.diff(p0, X), -X*p0)
chk("psi0'' = (x^2 - 1) psi0", sp.diff(p0, X, 2), (X**2 - 1)*p0)
chk("E0 = 1/2", sp.simplify(H(p0, 1)/p0), sp.Rational(1, 2))
# ข้อ 4: psi1 = x e^{-x^2/2}
p1 = X*sp.exp(-X**2/2)
chk("psi1'' = (x^2 - 3) psi1", sp.diff(p1, X, 2), (X**2 - 3)*p1)
chk("E1 = 3/2", sp.simplify(H(p1, 1)/p1), sp.Rational(3, 2))
chk("E1 - E0 = 1 = w", sp.Rational(3, 2) - sp.Rational(1, 2), 1)
# รูปทั่วไป psi'' = (x^2 + d) psi ⇒ E = -d/2
dd = sp.symbols('d')
chk("psi'' = (x^2 + d) psi ⇒ H psi = (-d/2) psi",
    -sp.Rational(1, 2)*(X**2 + dd) + sp.Rational(1, 2)*X**2, -dd/2)

# ── ③ ตัวเลขในข้อฝึก ───────────────────────────────────────
En = lambda ww, n: ww*(n + sp.Rational(1, 2))
chk("ข้อ 1: d^2/dx^2 cos(3x) = -9 cos(3x)", sp.diff(sp.cos(3*x), x, 2), -9*sp.cos(3*x))
chk("ข้อ 1: k=3, m=4 ⇒ w^2 = 25", 3**2 + 4**2, 25)
chk("ข้อ 1: w = 5", sp.sqrt(25), 5)
chk("ข้อ 2: k=0, m=3 ⇒ w = 3", sp.sqrt(0 + 9), 3)
chk("ข้อ 2: k=4, m=3 ⇒ w = 5", sp.sqrt(16 + 9), 5)
chk("ข้อ 2: k=2, m=0 ⇒ w = 2", sp.sqrt(4 + 0), 2)
chk("ข้อ 5: w=2 ⇒ E0 = 1", En(2, 0), 1)
chk("ข้อ 5: w=2 ⇒ E3 = 7", En(2, 3), 7)
chk("ข้อ 5: E5 - E4 = 2", En(2, 5) - En(2, 4), 2)
chk("ข้อ 6: ขั้นเดียว k=4, m=3 = 5", sp.sqrt(4**2 + 3**2), 5)
chk("ข้อ 6: E^2 = p^2 + m^2 ที่ p=4, m=3 ⇒ E = 5", sp.sqrt(4**2 + 3**2), 5)
chk("ข้อ 6: เหนือพื้น 15 ⇒ n = 3", sp.Rational(15, 5), 3)
chk("ข้อ 7: 2 อนุภาคของ w=3 + 1 อนุภาคของ w=5 = 11", 2*3 + 1*5, 11)
chk("ข้อ 7: จำนวนอนุภาค = 3", 2 + 1, 3)
chk("ข้อ 7 (แก้): ทั้งหมดตามสมการ (3) = 15", En(3, 2) + En(5, 1), 15)
chk("ข้อ 7 (แก้): สุญญากาศ = 4", En(3, 0) + En(5, 0), 4)
chk("ข้อ 7 (แก้): 15 - 4 = 11", (En(3, 2) + En(5, 1)) - (En(3, 0) + En(5, 0)), 11)
chk("ตัวอย่างในบท: ทั้งหมด 18", En(2, 1) + En(6, 2), 18)
chk("ตัวอย่างในบท: สุญญากาศ 4", En(2, 0) + En(6, 0), 4)
chk("ตัวอย่างในบท: เทียบสุญญากาศ 14 = 1*2 + 2*6", (En(2, 1) + En(6, 2)) - (En(2, 0) + En(6, 0)), 1*2 + 2*6)
chk("ข้อ 7: k=0 ⇒ พลังงานหนึ่งขั้น = m = 3", sp.sqrt(0 + 3**2), 3)

bad_ = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
