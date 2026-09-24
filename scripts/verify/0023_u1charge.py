"""ตรวจทุกข้ออ้างของบทที่ 11 — สมมาตรภายใน U(1) ⇒ ประจุ

signature g = diag(-1, +1) ใน 1+1 มิติ (เหมือนบท 8-10)
⚠️ ห้ามพิมพ์สูตร j^mu จากความจำ — ต้อง derive จาก L ทุกครั้ง
   (ตำราเขียนเครื่องหมาย/ตัวประกอบ i ต่างกัน)
"""
import sympy as sp

t, x, m, al = sp.symbols('t x m alpha', real=True)
I = sp.I

# สนามเชิงซ้อนกับสังยุคของมัน — ปฏิบัติเป็น "สองสนามอิสระ" ตามวิธีมาตรฐาน
phi = sp.Function('phi')(t, x)
phs = sp.Function('phis')(t, x)          # phi^*

pt,  px  = sp.diff(phi, t), sp.diff(phi, x)
pts, pxs = sp.diff(phs, t), sp.diff(phs, x)

ok = []
def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))

# ── ① ลากรองจ์ ─────────────────────────────────────────────
# L = -g^{mu nu} d_mu phi* d_nu phi - m^2 phi* phi
L = pts*pt - pxs*px - m**2*phs*phi
chk("L ที่กางออกแล้ว", L, pts*pt - pxs*px - m**2*phs*phi)

# สนามจริง (phi* = phi) ต้องยุบเป็นลากรองจ์ของบทที่ 8-10 คูณสอง
Lreal = L.subs({phs: phi}).doit()
phi_t, phi_x = sp.diff(phi, t), sp.diff(phi, x)
chk("phi* = phi ⇒ L = 2 x (ลากรองจ์สนามจริง)",
    Lreal, 2*(sp.Rational(1,2)*phi_t**2 - sp.Rational(1,2)*phi_x**2 - sp.Rational(1,2)*m**2*phi**2))

# ── ② สมมาตร: phi -> e^{i alpha} phi ──────────────────────
# delta phi = i phi · delta phi* = -i phi*   (จากการกระจาย e^{i alpha} ถึงอันดับหนึ่ง)
e = sp.series(sp.exp(I*al), al, 0, 2).removeO()
chk("กระจาย e^{i alpha} ถึงอันดับหนึ่ง = 1 + i alpha", e, 1 + I*al)
dphi, dphs = I*phi, -I*phs

# delta L ต้องเป็นศูนย์ **เอง** ไม่ต้องใช้สมการการเคลื่อนที่ (off-shell)
dL = (sp.diff(L, phi)*dphi + sp.diff(L, phs)*dphs
      + sp.diff(L, pt)*sp.diff(dphi, t) + sp.diff(L, px)*sp.diff(dphi, x)
      + sp.diff(L, pts)*sp.diff(dphs, t) + sp.diff(L, pxs)*sp.diff(dphs, x))
chk("delta L = 0 แบบ off-shell (ไม่ต้องใช้ EOM)", sp.expand(dL), 0)

# ── ③ กระแส — ต้องบวกสองสนาม (sum_i กลับมา) ───────────────
jt = sp.diff(L, pt)*dphi + sp.diff(L, pts)*dphs
jx = sp.diff(L, px)*dphi + sp.diff(L, pxs)*dphs
chk("dL/d(pt)  = pts", sp.diff(L, pt), pts)
chk("dL/d(pts) = pt",  sp.diff(L, pts), pt)
chk("dL/d(px)  = -pxs", sp.diff(L, px), -pxs)
chk("dL/d(pxs) = -px",  sp.diff(L, pxs), -px)

chk("j^t = i*(phi*pts - phis*pt)", jt, I*(phi*pts - phs*pt))
chk("j^x = i*(phis*px - phi*pxs)", jx, I*(phs*px - phi*pxs))
# พจน์รายตัวที่ข้อฝึกถาม
chk("พจน์แรกของ j^t  = i*phi*pts",  sp.diff(L, pt)*dphi, I*phi*pts)
chk("พจน์ที่สองของ j^t = -i*phis*pt", sp.diff(L, pts)*dphs, -I*phs*pt)

# ── ④ สมการการเคลื่อนที่ ───────────────────────────────────
def EL(field, d_t, d_x):
    return sp.diff(sp.diff(L, d_t), t) + sp.diff(sp.diff(L, d_x), x) - sp.diff(L, field)

ptt,  pxx  = sp.diff(phi, t, 2), sp.diff(phi, x, 2)
ptts, pxxs = sp.diff(phs, t, 2), sp.diff(phs, x, 2)
chk("EL ของ phi*  ⇒ ptt - pxx + m^2 phi", EL(phs, pts, pxs), ptt - pxx + m**2*phi)
chk("EL ของ phi   ⇒ ptts - pxxs + m^2 phis", EL(phi, pt, px), ptts - pxxs + m**2*phs)

# ── ⑤ off-shell / on-shell (กรอบของบทที่ 10) ───────────────
div = sp.expand(sp.diff(jt, t) + sp.diff(jx, x))
chk("off-shell: d_mu j^mu = i*(phi*(EOM ของ phi*) - phis*(EOM ของ phi))",
    div, I*(phi*(ptts - pxxs + m**2*phs) - phs*(ptt - pxx + m**2*phi)))
EOM = {ptt: pxx - m**2*phi, ptts: pxxs - m**2*phs}
chk("on-shell: d_mu j^mu = 0", div.subs(EOM), 0)
# ⚠️ ไม่ใช่ศูนย์เองแบบไม่มีเงื่อนไข — ต้องยืนยันว่า off-shell ไม่เป็นศูนย์
chk("off-shell ไม่เป็นศูนย์: phi=t^2, phis=x^2, m=0 ได้ -2i(t^2 + x^2)",
    div.subs(m, 0).subs({phi: t**2, phs: x**2}).doit(),
    I*(t**2*(-2) - x**2*(2)) )

# ── ⑥ สนามจริงไม่มีประจุ — เป็นศูนย์เอง ไม่ต้องใช้ EOM ─────
chk("phi* = phi ⇒ j^t = 0 เอง", jt.subs({phs: phi}).doit(), 0)
chk("phi* = phi ⇒ j^x = 0 เอง", jx.subs({phs: phi}).doit(), 0)

# ── ⑦ ตัวเลขที่ข้อฝึกใช้ ───────────────────────────────────
def bracket(v_phi, v_pts, v_phis, v_pt):     # ก้อนในวงเล็บของ j^t
    return v_phi*v_pts - v_phis*v_pt

chk("วงเล็บที่ phi=1, pts=4, phis=3, pt=2 ⇒ -2", bracket(1, 4, 3, 2), -2)
chk("วงเล็บที่ phi=3, pts=5, phis=3, pt=5 ⇒ 0 (สนามจริง)", bracket(3, 5, 3, 5), 0)
chk("วงเล็บที่ phi=2, pts=7, phis=5, pt=1 ⇒ 9", bracket(2, 7, 5, 1), 9)
chk("j^t ที่วงเล็บ = -2 คือ -2i", I*(-2), -2*I)

# ── ⑧ มิติ (เลขยกกำลังของมวล ใน 4 มิติ) ────────────────────
chk("[phi] = 1 เท่าสนามจริง", 1, 1)
chk("[j^mu] = 3", (1 + 1) + 1, 3)   # [d phi] + [phi] = 2 + 1
chk("[Q] = [j^0] - 3 = 0 (ประจุไร้มิติ)", 3 - 3, 0)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
