"""ตรวจข้ออ้างของบทที่ 10 — on-shell / off-shell

จุดที่บทที่ 9 ยังไม่ได้ตรวจ: **ค่า off-shell** (ยังไม่แทนสมการการเคลื่อนที่)
บทที่ 9 ตรวจแต่ค่า on-shell ⇒ ตารางสองคอลัมน์ของบทที่ 10 ต้องตรวจคอลัมน์ซ้ายด้วย
"""
import sympy as sp

t, x, m = sp.symbols('t x m', real=True)
phi = sp.Function('phi')(t, x)
pt, px = sp.diff(phi, t), sp.diff(phi, x)
ptt, pxx = sp.diff(phi, t, 2), sp.diff(phi, x, 2)

ok = []
def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))

L_free = sp.Rational(1, 2)*pt**2 - sp.Rational(1, 2)*px**2
L_mass = L_free - sp.Rational(1, 2)*m**2*phi**2
EOM_mass = {ptt: pxx - m**2*phi}          # ptt - pxx + m^2 phi = 0
EOM_free = {ptt: pxx}

def div_j(L):
    """d_mu j^mu ของกระแสเลื่อนค่าสนาม — ยังไม่แทน EOM (= ค่า off-shell)"""
    return sp.diff(sp.diff(L, pt), t) + sp.diff(sp.diff(L, px), x)

# ── ① ตารางสองคอลัมน์ของ d_mu j^mu ─────────────────────────
# ⚠️ หัวใจของบท: คอลัมน์ off-shell **เหมือนกันทั้งสองแถว**
chk("ไร้มวล off-shell = ptt - pxx", div_j(L_free), ptt - pxx)
chk("มีมวล  off-shell = ptt - pxx", div_j(L_mass), ptt - pxx)
chk("off-shell ของสองแถวเท่ากันจริง", div_j(L_free) - div_j(L_mass), 0)
chk("ไร้มวล on-shell = 0", div_j(L_free).subs(EOM_free), 0)
chk("มีมวล  on-shell = -m^2 phi", div_j(L_mass).subs(EOM_mass), -m**2*phi)

# ── ② off-shell ไม่เป็นศูนย์จริง ๆ (ไม่ใช่ศูนย์ที่มองไม่ออก) ──
# ถ้าเผลอเป็นศูนย์อยู่แล้ว ความต่างทั้งบทจะไม่มีอยู่จริง
chk("off-shell ไม่เป็นศูนย์: แทน phi = t^2 ได้ 2", (ptt - pxx).subs(phi, t**2).doit(), 2)
chk("off-shell ไม่เป็นศูนย์: แทน phi = x^3 ได้ -6x", (ptt - pxx).subs(phi, x**3).doit(), -6*x)
# สนามที่เดินตาม EOM ไร้มวลจริง (คลื่นวิ่ง) ⇒ ได้ศูนย์
chk("phi = f(x - t) เป็นคำตอบของ EOM ไร้มวล",
    (ptt - pxx).subs(phi, sp.Function('f')(x - t)).doit(), 0)

# ── ③ T^{mu 0} — ค่า off-shell แบบแยกตัวประกอบ ─────────────
g = sp.diag(-1, 1)
dphi = [pt, px]
coord = [t, x]
dup = [sum(g[mu, a]*dphi[a] for a in range(2)) for mu in range(2)]
T = sp.Matrix(2, 2, lambda mu, nu: dup[mu]*dup[nu] + g[mu, nu]*L_mass)

divT0 = sp.expand(sum(sp.diff(T[mu, 0], coord[mu]) for mu in range(2)))
chk("d_mu T^mu0 off-shell = pt*(ptt - pxx + m^2 phi)", divT0, pt*(ptt - pxx + m**2*phi))
chk("d_mu T^mu0 on-shell = 0", divT0.subs(EOM_mass), 0)

divT1 = sp.expand(sum(sp.diff(T[mu, 1], coord[mu]) for mu in range(2)))
# ⚠️ เดาไว้ว่าเป็น +px*(...) — sympy ตีกลับ ของจริงติดลบ
#    (เครื่องหมายจากเมตริกอีกครั้ง · ห้ามเดาจากรูปของแถวบน)
chk("d_mu T^mu1 off-shell = -px*(ptt - pxx + m^2 phi)", divT1, -px*(ptt - pxx + m**2*phi))
chk("d_mu T^mu1 on-shell = 0", divT1.subs(EOM_mass), 0)

# ⇒ ทั้งสองค่าของ nu มีวงเล็บ "สมการการเคลื่อนที่" เป็นตัวประกอบร่วม
#   ต่างกันแค่ตัวคูณหน้า (pt กับ -px) — ข้ออ้างนี้ต้องตรวจ ไม่ใช่เขียนลอย ๆ
chk("สองแถวต่างกันแค่ตัวคูณหน้า", sp.simplify(divT0*px + divT1*pt), 0)

# ── ④ ตัวเลขที่ข้อฝึกใช้ ───────────────────────────────────
chk("อะไรก็ตามคูณศูนย์ = ศูนย์ (pt = 7)", 7*0, 0)
chk("วงเล็บ = 3 แล้ว pt = 2 ⇒ ได้ 6", 2*3, 6)
chk("off-shell ของ T ที่วงเล็บ = 5 และ pt = 4 ⇒ 20", 4*5, 20)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
