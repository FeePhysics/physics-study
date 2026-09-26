"""ตรวจทุกข้ออ้างของบทที่ 12 — ให้ alpha ขึ้นกับตำแหน่ง ⇒ ต้องมี A_mu

convention ของบทนี้ (ตำราอื่นอาจกลับเครื่องหมาย — เลยต้อง derive ทุกอย่างที่นี่):
    D_mu = d_mu - i e A_mu        A_mu -> A_mu + (1/e) d_mu alpha
signature g = diag(-1, +1) ใน 1+1 มิติ · ดัชนีล่าง A_mu = (At, Ax)
"""
import sympy as sp

t, x, m, e = sp.symbols('t x m e', real=True, nonzero=True)
eps = sp.symbols('epsilon')                 # ตัวนับอันดับของการแปลง
I = sp.I

phi = sp.Function('phi')(t, x)
phs = sp.Function('phis')(t, x)
al  = sp.Function('alpha')(t, x)            # ⚠️ ขึ้นกับตำแหน่ง — หัวใจของบท
At  = sp.Function('At')(t, x)
Ax  = sp.Function('Ax')(t, x)

ok = []
def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))

def first_order(expr):
    """สัมประสิทธิ์ของ eps^1 — คือ delta(...) ของการแปลง"""
    return sp.expand(sp.diff(expr, eps).subs(eps, 0))

# การแปลงที่อันดับ eps (alpha -> eps*alpha)
phi_T = phi*sp.exp(I*eps*al)
phs_T = phs*sp.exp(-I*eps*al)
At_T  = At + eps*sp.diff(al, t)/e
Ax_T  = Ax + eps*sp.diff(al, x)/e

# ── ① ปัญหา: อนุพันธ์ธรรมดาได้พจน์เกิน ──────────────────────
d_pt = first_order(sp.diff(phi_T, t))
chk("delta(pt) = i*alpha*pt + i*(d_t alpha)*phi",
    d_pt, I*al*sp.diff(phi, t) + I*sp.diff(al, t)*phi)
chk("พจน์เกิน = i*(d_t alpha)*phi",
    d_pt - I*al*sp.diff(phi, t), I*sp.diff(al, t)*phi)
# alpha คงที่ ⇒ พจน์เกินหาย (กลับไปเป็นบทที่ 11)
c = sp.symbols('c', real=True)
chk("alpha คงที่ ⇒ delta(pt) = i*alpha*pt",
    first_order(sp.diff(phi*sp.exp(I*eps*c), t)), I*c*sp.diff(phi, t))

# ── ② ซ่อมด้วยอนุพันธ์โคแวเรียนต์ ───────────────────────────
def Dt(f, A): return sp.diff(f, t) - I*e*A*f
def Dx(f, A): return sp.diff(f, x) - I*e*A*f

# ชิ้นที่ A ขนมา: -i e (delta A_t) phi
chk("-i e (delta At) phi = -i*(d_t alpha)*phi",
    -I*e*(sp.diff(al, t)/e)*phi, -I*sp.diff(al, t)*phi)
chk("พจน์เกิน + ชิ้นของ A = 0",
    I*sp.diff(al, t)*phi + (-I*sp.diff(al, t)*phi), 0)

d_Dt = first_order(Dt(phi_T, At_T))
chk("delta(D_t phi) = i*alpha*(D_t phi) — แปลงตัวเหมือน phi เป๊ะ",
    d_Dt, I*al*Dt(phi, At))
d_Dx = first_order(Dx(phi_T, Ax_T))
chk("delta(D_x phi) = i*alpha*(D_x phi)", d_Dx, I*al*Dx(phi, Ax))

# ── ③ ลากรองจ์ทั้งก้อนไม่เปลี่ยน ─────────────────────────────
def Dt_s(f, A): return sp.diff(f, t) + I*e*A*f      # สังยุคของ D (A จริง)
def Dx_s(f, A): return sp.diff(f, x) + I*e*A*f

def Lag(p, ps, a_t, a_x):
    return (Dt_s(ps, a_t)*Dt(p, a_t) - Dx_s(ps, a_x)*Dx(p, a_x) - m**2*ps*p)

L0 = Lag(phi, phs, At, Ax)
chk("delta L = 0 แม้ alpha ขึ้นกับตำแหน่ง", first_order(Lag(phi_T, phs_T, At_T, Ax_T)), 0)

# ถ้าไม่มี A (อนุพันธ์ธรรมดา) ⇒ พัง
L_noA = lambda p, ps: sp.diff(ps, t)*sp.diff(p, t) - sp.diff(ps, x)*sp.diff(p, x) - m**2*ps*p
broken = first_order(L_noA(phi_T, phs_T))
ok.append(("ไม่มี A ⇒ delta L ไม่เป็นศูนย์ (พังจริง)", sp.simplify(broken) != 0, broken, "≠0"))

# ── ④ F ไม่เปลี่ยน — ใช้อนุพันธ์สลับลำดับได้ (บทที่ 9) ─────────
F_tx = lambda a_t, a_x: sp.diff(a_x, t) - sp.diff(a_t, x)
chk("delta F_tx = 0", first_order(F_tx(At_T, Ax_T)), 0)
chk("d_t(d_x alpha) - d_x(d_t alpha) = 0",
    sp.diff(sp.diff(al, x), t) - sp.diff(sp.diff(al, t), x), 0)

# ── ⑤ กระจาย L: พจน์เชิงเส้นใน A คือ -e A_mu j^mu ─────────────
pt, px   = sp.diff(phi, t), sp.diff(phi, x)
pts, pxs = sp.diff(phs, t), sp.diff(phs, x)
jt = I*(phi*pts - phs*pt)                    # จากบทที่ 11
jx = I*(phs*px - phi*pxs)
a, b = sp.symbols('a b')                      # A_t, A_x ค่าคงที่ ณ จุดหนึ่ง (พอสำหรับพีชคณิต)
Lab = Lag(phi, phs, a, b)
lin = sp.expand(sp.diff(Lab, a).subs({a: 0, b: 0})*a + sp.diff(Lab, b).subs({a: 0, b: 0})*b)
chk("พจน์เชิงเส้นใน A = -e (A_t j^t + A_x j^x)", lin, -e*(a*jt + b*jx))
chk("สัมประสิทธิ์ของ A_t j^t = -e", sp.expand(sp.diff(Lab, a).subs({a: 0, b: 0}) / jt), -e)
quad = sp.expand(Lab - Lab.subs({a: 0, b: 0}) - lin)
chk("พจน์กำลังสอง = e^2 (A_t^2 - A_x^2) phis*phi", quad, e**2*(a**2 - b**2)*phs*phi)
chk("พจน์ไม่มี A = L ของบทที่ 11", sp.expand(Lab.subs({a: 0, b: 0})), pts*pt - pxs*px - m**2*phs*phi)
# A_t^2 - A_x^2 = -A_mu A^mu (signature -+)
chk("A_t^2 - A_x^2 = -(A^mu A_mu)", a**2 - b**2, -(-a**2 + b**2))

# ── ⑥ ตัวเลขในข้อฝึก ───────────────────────────────────────
chk("e=2, dta=6 ⇒ delta At = 3", sp.Integer(6)/2, 3)
chk("e=2, At=3, j^t=5 ⇒ -e At j^t = -30", -2*3*5, -30)
chk("e=2, At=1, phis*phi=3 ⇒ e^2 At^2 |phi|^2 = 12", 2**2*1**2*3, 12)
chk("e=1, At=2, j^t=4, Ax=1, j^x=3 ⇒ -e(At jt + Ax jx) = -11", -1*(2*4 + 1*3), -11)

# ── ⑦ มิติใน 3+1 ──────────────────────────────────────────
d_del, d_A = 1, 1
chk("[e A] = [d] = 1", d_del, 1)
chk("[e] = [d] - [A] = 0 — ประจุไร้มิติใน 4 มิติ", d_del - d_A, 0)
chk("A_mu ใน 3+1 มี 4 ส่วนประกอบ", sp.diag(-1, 1, 1, 1).shape[0], 4)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
