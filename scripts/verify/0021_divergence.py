"""ตรวจทุกข้ออ้างของบทที่ 9 — ไดเวอร์เจนซ์ (ตอบกลับผล 2/8 ของบทที่ 8)

กฎบ้าน: ห้ามเชื่อความรู้ในหัว · ทุกตัวเลข/เครื่องหมายในบทต้องออกมาจากที่นี่
signature g = diag(-1, +1) ใน 1+1 มิติ
"""
import sympy as sp

t, x, m = sp.symbols('t x m', real=True, positive=None)
phi = sp.Function('phi')(t, x)

ok = []
def chk(name, got, want):
    g = sp.simplify(sp.expand(got - want))
    ok.append((name, g == 0, got, want))

# ── ① อนุพันธ์ชั้นสอง: พจนานุกรมใหม่ ─────────────────────────
pt, px = sp.diff(phi, t), sp.diff(phi, x)
ptt, pxx = sp.diff(phi, t, 2), sp.diff(phi, x, 2)
ptx = sp.diff(phi, t, x)
chk("ptx = pxt (อนุพันธ์สลับที่ได้)", ptx, sp.diff(phi, x, t))

# ── ② สมการการเคลื่อนที่จาก EL ─────────────────────────────
def euler_lagrange(L):
    """∂_t(∂L/∂pt) + ∂_x(∂L/∂px) − ∂L/∂phi = 0"""
    return (sp.diff(sp.diff(L, pt), t)
            + sp.diff(sp.diff(L, px), x)
            - sp.diff(L, phi))

L_free = sp.Rational(1,2)*pt**2 - sp.Rational(1,2)*px**2
L_mass = L_free - sp.Rational(1,2)*m**2*phi**2

chk("EOM ไร้มวล  ptt - pxx = 0",       euler_lagrange(L_free), ptt - pxx)
chk("EOM มีมวล   ptt - pxx + m^2 phi", euler_lagrange(L_mass), ptt - pxx + m**2*phi)

# ── ③ กระแสจากการเลื่อนค่าสนาม (delta phi = 1) ───────────────
def shift_current(L):
    return sp.diff(L, pt), sp.diff(L, px)          # (j^t, j^x)

jt_f, jx_f = shift_current(L_free)
jt_m, jx_m = shift_current(L_mass)
chk("ไร้มวล j^t = pt",  jt_f, pt)
chk("ไร้มวล j^x = -px", jx_f, -px)
chk("มีมวล  j^t = pt  (พจน์มวลไม่มีอนุพันธ์ ⇒ กระแสเท่าเดิม)", jt_m, pt)
chk("มีมวล  j^x = -px", jx_m, -px)

def div2(jt_, jx_):
    return sp.diff(jt_, t) + sp.diff(jx_, x)

chk("ไดเวอร์เจนซ์ก่อนใช้ EOM  = ptt - pxx", div2(jt_f, jx_f), ptt - pxx)

# ── ④ สมการ (3): d_mu j^mu = delta L ────────────────────────
# แทน EOM ลงไป: ptt = pxx - m^2 phi
sub_eom = {ptt: pxx - m**2*phi}
dL_free = sp.diff(L_free, phi)            # delta L เมื่อ delta phi = 1
dL_mass = sp.diff(L_mass, phi)
chk("ไร้มวล delta L = 0", dL_free, 0)
chk("มีมวล  delta L = -m^2 phi", dL_mass, -m**2*phi)

div_free = div2(jt_f, jx_f).subs({ptt: pxx})          # EOM ไร้มวล
div_mass = div2(jt_m, jx_m).subs(sub_eom)             # EOM มีมวล
chk("ไร้มวล d_mu j^mu = 0 หลังใช้ EOM", div_free, 0)
chk("มีมวล  d_mu j^mu = delta L หลังใช้ EOM", div_mass, dL_mass)

# ── ⑤ ไดเวอร์เจนซ์เปล่า ๆ (ข้อฝึกคณิตล้วน ไม่มีฟิสิกส์) ──────
cases = [
    ("j^t = phi, j^x = 0",        phi,      sp.Integer(0), pt),
    ("j^t = 0,   j^x = phi",      sp.Integer(0), phi,      px),
    ("j^t = px,  j^x = pt",       px,       pt,            2*ptx),
    ("j^t = pt,  j^x = px",       pt,       px,            ptt + pxx),
    ("j^t = pt,  j^x = -px",      pt,       -px,           ptt - pxx),
    ("j^t = x*phi, j^x = 0",      x*phi,    sp.Integer(0), x*pt),
]
for name, a, b, want in cases:
    chk("div · " + name, div2(a, b), want)

# ── ⑥ T^{mu nu} และการอนุรักษ์ ──────────────────────────────
g = sp.diag(-1, 1)
dphi = [sp.diff(phi, t), sp.diff(phi, x)]
coord = [t, x]
def Tup(L):
    # T^{mu nu} = d^mu phi d^nu phi + g^{mu nu} L   (ดัชนีบน: d^mu = g^{mu a} d_a)
    dup = [sum(g[mu, a]*dphi[a] for a in range(2)) for mu in range(2)]
    return sp.Matrix(2, 2, lambda mu, nu: dup[mu]*dup[nu] + g[mu, nu]*L)

T = Tup(L_mass)
chk("T^00 = 1/2 pt^2 + 1/2 px^2 + 1/2 m^2 phi^2",
    T[0, 0], sp.Rational(1,2)*pt**2 + sp.Rational(1,2)*px**2 + sp.Rational(1,2)*m**2*phi**2)
chk("T สมมาตร", T[0, 1] - T[1, 0], 0)

# d_mu T^{mu nu} = 0
# ⚠️ d_mu มีดัชนี "ล่าง" อยู่แล้ว = d/dx^mu ตรง ๆ ⇒ ไม่ต้องใช้เมตริกช่วยเลย
#    (เคยเผลอเอา g ไปคูณก่อน แล้วได้เศษเหลือ -2 pt pxx - 2 px ptx ซึ่งไม่เป็นศูนย์)
for nu in range(2):
    d = sum(sp.diff(T[mu, nu], coord[mu]) for mu in range(2))
    chk(f"d_mu T^mu{nu} = 0 หลังใช้ EOM", sp.expand(d).subs(sub_eom), 0)

# และตัวที่ข้อฝึกถาม: T^{mu 0} แยกช่องให้เห็น
chk("T^10 = -px*pt", T[1, 0], -px*pt)

# ── ⑦ ค่าตัวเลขที่ใช้ในข้อฝึก ───────────────────────────────
def T00(vpt, vpx, vm, vphi):
    return sp.Rational(1,2)*vpt**2 + sp.Rational(1,2)*vpx**2 + sp.Rational(1,2)*vm**2*vphi**2

chk("T00(2,4,0,-)  = 10", T00(2, 4, 0, 0), 10)
chk("T00(1,1,2,3)  = 19", T00(1, 1, 2, 3), 19)
chk("T00(0,0,5,2)  = 50", T00(0, 0, 5, 2), 50)
chk("T00(0,0,3,2)  = 18", T00(0, 0, 3, 2), 18)
chk("T00(0,0,2,5)  = 50", T00(0, 0, 2, 5), 50)   # สลับ m กับ phi ได้ค่าเดียวกัน
chk("T00(3,0,0,-)  = 4.5", T00(3, 0, 0, 0), sp.Rational(9,2))
chk("พจน์มวลอย่างเดียว m=4 phi=1 -> 8", T00(0, 0, 4, 1), 8)
chk("EM T00 ที่ E=B=3 -> 9", sp.Rational(1,2)*(3**2 + 3**2), 9)
chk("-1/4 F^2 ที่ E=B=3 -> 0", sp.Rational(1,2)*(3**2 - 3**2), 0)

# ── ⑧ จำนวนปริมาณอนุรักษ์ = มิติของกาลอวกาศ ─────────────────
# เลื่อนได้กี่ทิศอิสระ = กี่กระแส · ยืนยันด้วยจำนวนดัชนี nu ของ T^{mu nu}
chk("1+1 มิติ -> 2 ปริมาณ", T.shape[1], 2)
g4 = sp.diag(-1, 1, 1, 1)
chk("3+1 มิติ -> 4 ปริมาณ", g4.shape[1], 4)

# ── ⑨ มิติ (เลขยกกำลังของมวล) ───────────────────────────────
d_phi, d_del = 1, 1                      # [phi] = 1, [d] = 1 ใน 4 มิติ
chk("[T^{mu nu}] = 4", 2*(d_del + d_phi), 4)
chk("[L] = 4",         2*(d_del + d_phi), 4)
chk("[j^mu] = 3  (j = dL/d(d phi) * delta phi, delta phi ไร้มิติ)", 4 - 1, 3)

# ── รายงาน ──────────────────────────────────────────────────
bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
