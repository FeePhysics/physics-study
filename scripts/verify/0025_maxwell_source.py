"""ตรวจทุกข้ออ้างของบทที่ 13 — สมการของ A_mu เอง (แมกซ์เวลล์ที่มีแหล่ง)

convention เดียวกับบทที่ 12: D = d - i e A · signature (-,+) · 1+1 มิติ
F = F_tx = d_t Ax - d_x At (มีช่องเดียวใน 1+1)
"""
import sympy as sp

t, x, m, e = sp.symbols('t x m e', real=True, nonzero=True)
I = sp.I
phi = sp.Function('phi')(t, x); phs = sp.Function('phis')(t, x)
At = sp.Function('At')(t, x);   Ax = sp.Function('Ax')(t, x)

ok = []
def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))

# สัญลักษณ์แทนอนุพันธ์ของ A เพื่อหาอนุพันธ์ย่อยของ L เทียบ d_mu A_nu
dtAt, dxAt, dtAx, dxAx = sp.symbols('dtAt dxAt dtAx dxAx')
a_t, a_x = sp.symbols('a_t a_x')           # A เอง (ไม่ใช่อนุพันธ์)
P, Ps, Pt, Pts, Px, Pxs = sp.symbols('P Ps Pt Pts Px Pxs')

Fsym = dtAx - dxAt
# -1/4 F_{mu nu} F^{mu nu} ใน 1+1 (-,+): F^{tx} = g^tt g^xx F_tx = -F_tx
g = sp.diag(-1, 1)
Fdn = sp.Matrix([[0, Fsym], [-Fsym, 0]])
Fup = g*Fdn*g
LF = -sp.Rational(1, 4)*sum(Fdn[i, j]*Fup[i, j] for i in range(2) for j in range(2))
chk("-1/4 F^2 = 1/2 F_tx^2 ใน 1+1", LF, sp.Rational(1, 2)*Fsym**2)
chk("F_tx = 4 ⇒ -1/4 F^2 = 8", (sp.Rational(1, 2)*Fsym**2).subs({dtAx: 4, dxAt: 0}), 8)

Dt, Dts = Pt - I*e*a_t*P, Pts + I*e*a_t*Ps
Dx, Dxs = Px - I*e*a_x*P, Pxs + I*e*a_x*Ps
L = LF + Dts*Dt - Dxs*Dx - m**2*Ps*P

# ── ① อนุพันธ์ย่อยฝั่ง F ────────────────────────────────────
chk("dF/d(dxAt) = -1", sp.diff(Fsym, dxAt), -1)
chk("dL/d(dxAt) = -F", sp.diff(L, dxAt), -Fsym)
chk("dL/d(dtAx) = +F", sp.diff(L, dtAx), Fsym)
chk("dL/d(dtAt) = 0 (ไม่มี d_t At ใน F)", sp.diff(L, dtAt), 0)

# ── ② อนุพันธ์ย่อยฝั่งสสาร ─────────────────────────────────
jt = I*(P*Pts - Ps*Pt)                     # j ของบทที่ 11
jx = I*(Ps*Px - P*Pxs)
Jt = I*(P*Dts - Ps*Dt)                     # j ที่เปลี่ยน d เป็น D
Jx = I*(Ps*Dx - P*Dxs)
chk("J^t = j^t - 2 e At |phi|^2", Jt, jt - 2*e*a_t*Ps*P)
chk("dL/dAt = -e J^t", sp.diff(L, a_t), -e*Jt)
# ⚠️ เดาไว้ว่า +e J^x — sympy ตีกลับ: ของจริง -e J^x เหมือนช่อง t
#    (ทั้งสองช่องคือ dL/dA_mu = -e J^mu ตรงกับรูป -e A_mu J^mu)
chk("dL/dAx = -e J^x", sp.diff(L, a_x), -e*Jx)
chk("ส่วนที่ A = 0: dL/dAt = -e j^t", sp.diff(L, a_t).subs(a_t, 0), -e*jt)
chk("อนุพันธ์ของพจน์คู่ควบ -e At jt เทียบ At = -e jt", sp.diff(-e*a_t*jt, a_t), -e*jt)

# ── ③ สมการ EL ของ A ⇒ แมกซ์เวลล์ในรูป 1+1 ────────────────
# แปลงกลับเป็นฟังก์ชันจริงเพื่อหาอนุพันธ์รวม
F = sp.diff(Ax, t) - sp.diff(At, x)
back = {dtAx: sp.diff(Ax, t), dxAt: sp.diff(At, x), dtAt: sp.diff(At, t), dxAx: sp.diff(Ax, x),
        a_t: At, a_x: Ax, P: phi, Ps: phs,
        Pt: sp.diff(phi, t), Pts: sp.diff(phs, t), Px: sp.diff(phi, x), Pxs: sp.diff(phs, x)}
def B(ex): return ex.subs(back)

EL_At = sp.diff(B(sp.diff(L, dtAt)), t) + sp.diff(B(sp.diff(L, dxAt)), x) - B(sp.diff(L, a_t))
EL_Ax = sp.diff(B(sp.diff(L, dtAx)), t) + sp.diff(B(sp.diff(L, dxAx)), x) - B(sp.diff(L, a_x))
chk("EOM ของ At ⇔ d_x F = e J^t (กฎของเกาส์)", EL_At, -(sp.diff(F, x) - e*B(Jt)))
chk("EOM ของ Ax ⇔ d_t F = -e J^x", EL_Ax, sp.diff(F, t) + e*B(Jx))

# ── ④ ความสอดคล้อง: ด้านซ้ายมีไดเวอร์เจนซ์เป็นศูนย์ off-shell ─────
chk("d_t(d_x F) - d_x(d_t F) = 0 โดยไม่ต้องใช้ EOM ใดเลย",
    sp.diff(sp.diff(F, x), t) - sp.diff(sp.diff(F, t), x), 0)
# ⇒ แทนสองสมการ: e (d_t J^t + d_x J^x) = 0 — ต้องใช้ EOM ของ phi ถึงเป็นศูนย์จริง
div = sp.expand(sp.diff(B(Jt), t) + sp.diff(B(Jx), x))
ptt, pxx = sp.diff(phi, t, 2), sp.diff(phi, x, 2)
ptts, pxxs = sp.diff(phs, t, 2), sp.diff(phs, x, 2)
EOMphi_s = sp.solve(sp.expand(B(sp.diff(L, Ps)) - sp.diff(B(sp.diff(L, Pts)), t) - sp.diff(B(sp.diff(L, Pxs)), x)), ptt)
EOMphs_s = sp.solve(sp.expand(B(sp.diff(L, P)) - sp.diff(B(sp.diff(L, Pt)), t) - sp.diff(B(sp.diff(L, Px)), x)), ptts)
on = sp.simplify(div.subs({ptt: EOMphi_s[0], ptts: EOMphs_s[0]}))
chk("d_mu J^mu = 0 on-shell (แม้มี A)", on, 0)
off = sp.simplify(div.subs(e, 1).subs({At: 0, Ax: 0}).subs({phi: t**2, phs: x**2}).doit())
ok.append(("d_mu J^mu ไม่เป็นศูนย์ off-shell (phi=t^2, phis=x^2, A=0)", off != 0, off, "≠0"))

# ── ⑤ ตัวเลขในข้อฝึก ───────────────────────────────────────
chk("e=1, jt=4 ⇒ -e jt = -4", -1*4, -4)
chk("e=2, jt=3 ⇒ d_x F = 6", 2*3, 6)
chk("e=1, At=2, |phi|^2=3 ⇒ 2 e At |phi|^2 = 12", 2*1*2*3, 12)
chk("jt=5 ⇒ J^t = 5 - 12 = -7", 5 - 12, -7)
chk("d_t J^t = 5 ⇒ d_x J^x = -5", -5, -5)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
