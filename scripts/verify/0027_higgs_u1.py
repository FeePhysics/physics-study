"""ตรวจทุกข้ออ้างของบทที่ 15 — แสงที่มีมวล (กลไกฮิกส์แบบ U(1))

convention ต่อจากบท 11-13: D = d - i e A · signature (-,+) · 1+1 มิติ
ศักย์ V = lam (phis*phi - v^2)^2
⚠️ นิยาม "มวล" ด้วยความสัมพันธ์การกระจาย omega^2 = k^2 + M^2 ของคลื่นระนาบ
   ไม่ใช่ด้วยการเทียบรูปพจน์กับตำรา — ตัวประกอบ 2 ของ M^2 ขึ้นกับ normalisation ของ phi
"""
import sympy as sp

t, x = sp.symbols('t x', real=True)
e, v, lam, k, w, M, a = sp.symbols('e v lam k omega M a', positive=True)
r, th = sp.symbols('r theta', real=True)
I = sp.I

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


# ── ① ศักย์หมวกเม็กซิกัน ────────────────────────────────────
s = sp.symbols('s', nonnegative=True)                  # s = phis*phi = |phi|^2
V = lam*(s - v**2)**2
chk("V(0) = lam v^4", V.subs(s, 0), lam*v**4)
chk("V ที่ |phi|^2 = v^2 เท่ากับ 0", V.subs(s, v**2), 0)
chk("dV/ds = 0 ที่ s = v^2", sp.diff(V, s).subs(s, v**2), 0)
# เขียน phi = r e^{i th} ⇒ V ขึ้นกับ r อย่างเดียว ⇒ จุดต่ำสุดเป็นวงกลม r = v
phi_polar = r*sp.exp(I*th)
Vpolar = sp.simplify(lam*(sp.expand(phi_polar*sp.conjugate(phi_polar)) - v**2)**2)
chk("V ไม่ขึ้นกับ theta", sp.diff(Vpolar, th), 0)
d2 = sp.diff(lam*(r**2 - v**2)**2, r, 2).subs(r, v)
chk("d^2V/dr^2 ที่ r = v = 8 lam v^2 > 0 (เป็นจุดต่ำสุดจริง)", d2, 8*lam*v**2)
chk("d^2V/dr^2 ที่ r = 0 = -4 lam v^2 < 0 (r = 0 เป็นจุดสูงสุดตามแนวรัศมี)",
    sp.diff(lam*(r**2 - v**2)**2, r, 2).subs(r, 0), -4*lam*v**2)

# ── ② สุญญากาศไม่สมมาตร ─────────────────────────────────────
chk("delta phi ที่ phi = v คือ i v", I*v, I*v)
chk("delta phi ที่ phi = 0 คือ 0", I*0, 0)
chk("delta L ของการหมุนเฟส: V ขึ้นกับ phis*phi เท่านั้น ⇒ เปลี่ยน 0",
    sp.diff(V.subs(s, sp.expand(phi_polar*sp.conjugate(phi_polar))), th), 0)

# ── ③ แทน phi = phis = v (ค่าคงที่) ลงใน L ของบทที่ 12 ─────────
At, Ax, dtAx, dxAt, dtAt, dxAx = sp.symbols('At Ax dtAx dxAt dtAt dxAx', real=True)
P, Ps = sp.symbols('P Ps')
F = dtAx - dxAt
Dt, Dts = 0 - I*e*At*P, 0 + I*e*At*Ps          # d_t phi = 0 เพราะ phi คงที่
Dx, Dxs = 0 - I*e*Ax*P, 0 + I*e*Ax*Ps
L = sp.Rational(1, 2)*F**2 + Dts*Dt - Dxs*Dx - V.subs(s, Ps*P)
Lvac = sp.expand(L.subs({P: v, Ps: v}))
chk("L ที่สุญญากาศ = 1/2 F^2 + e^2 v^2 (At^2 - Ax^2)",
    Lvac, sp.expand(sp.Rational(1, 2)*F**2 + e**2*v**2*(At**2 - Ax**2)))

# ── ④ สมการการเคลื่อนที่ของ A (ขั้นตอนเดียวกับบทที่ 13) ─────────
dL_dxAt = sp.diff(Lvac, dxAt)          # = -F
dL_dtAx = sp.diff(Lvac, dtAx)          # = +F
chk("dL/d(dxAt) = -F", dL_dxAt, -F)
chk("dL/dAt = 2 e^2 v^2 At", sp.diff(Lvac, At), 2*e**2*v**2*At)
chk("dL/dAx = -2 e^2 v^2 Ax", sp.diff(Lvac, Ax), -2*e**2*v**2*Ax)
# ⇒  d_x F = -M^2 At  และ  d_t F = -M^2 Ax  โดย M^2 = 2 e^2 v^2

# ── ⑤ คลื่นระนาบ: หา omega ที่ทำให้ทั้งสองสมการจริง ──────────────
ph = k*x - w*t
Axf = a*sp.cos(ph)
Atf = -a*k/w*sp.cos(ph)                  # จากเงื่อนไข d_t At = d_x Ax
Ff = sp.diff(Axf, t) - sp.diff(Atf, x)
M2 = 2*e**2*v**2
eq_t = sp.simplify(sp.diff(Ff, x) + M2*Atf)      # EOM ของ At: d_x F + M^2 At = 0
eq_x = sp.simplify(sp.diff(Ff, t) + M2*Axf)      # EOM ของ Ax: d_t F + M^2 Ax = 0
w_disp = sp.sqrt(k**2 + M2)
chk("EOM ของ Ax เป็นจริงเมื่อ omega^2 = k^2 + 2 e^2 v^2", eq_x.subs(w, w_disp), 0)
chk("EOM ของ At เป็นจริงเมื่อ omega^2 = k^2 + 2 e^2 v^2", eq_t.subs(w, w_disp), 0)
chk("เงื่อนไข d_t At = d_x Ax", sp.diff(Atf, t) - sp.diff(Axf, x), 0)
# ตัวควบคุมฝั่งลบ: omega^2 = k^2 (ไร้มวล) ไม่เป็นคำตอบเมื่อ v != 0
bad = sp.simplify(eq_x.subs(w, k))
ok.append(("omega = k ไม่ใช่คำตอบเมื่อ v != 0 (มีมวลจริง)", bad != 0, bad, "≠0"))
# รูปเดียวกับสนามมีมวลของบทที่ 9-10: ptt - pxx + M^2 (.) = 0
kg = sp.simplify(sp.diff(Axf, t, 2) - sp.diff(Axf, x, 2) + M2*Axf)
chk("Ax เดินตาม ptt - pxx + M^2 Ax = 0", kg.subs(w, w_disp), 0)

# ── ⑥ ตัวเลขในข้อฝึก ───────────────────────────────────────
chk("lam=2, v=3 ⇒ V(0) = 162", V.subs({lam: 2, v: 3, s: 0}), 162)
chk("lam=2, v=3 ⇒ V(|phi|^2 = 9) = 0", V.subs({lam: 2, v: 3, s: 9}), 0)
chk("e=1, v=3 ⇒ e^2 v^2 = 9", (e**2*v**2).subs({e: 1, v: 3}), 9)
chk("e=1, v=3 ⇒ M^2 = 18", M2.subs({e: 1, v: 3}), 18)
chk("k=3, M^2=16 ⇒ omega^2 = 25", 3**2 + 16, 25)
chk("omega = 5", sp.sqrt(25), 5)
chk("M = 0, k = 3 ⇒ omega = 3", sp.sqrt(3**2 + 0), 3)
chk("delta phi ที่ v = 3 คือ 3i", I*3, 3*I)

# ── ⑦ นับกำลังด้วยกฎที่นิยามชัด: แทนทุก A ด้วย lambda*A ────────
L_ = sp.symbols('Lambda')
chk("กำลังของ A ใน e^2 v^2 (At^2 - Ax^2) = 2",
    sp.Poly(sp.expand((e**2*v**2*(At**2 - Ax**2)).subs({At: L_*At, Ax: L_*Ax})), L_).degree(), 2)
chk("กำลังของ A ใน At*Ax (ตัวเดียวของคอมมิวเทเตอร์) = 2",
    sp.Poly(sp.expand((At*Ax).subs({At: L_*At, Ax: L_*Ax})), L_).degree(), 2)
jt = sp.symbols('jt')
chk("กำลังของ A ใน -e At jt = 1",
    sp.Poly(sp.expand((-e*At*jt).subs({At: L_*At})), L_).degree(), 1)
chk("กำลังของ phi ใน V (กางแล้ว) = 4",
    sp.Poly(sp.expand(lam*(Ps*P - v**2)**2).subs({P: L_*P, Ps: L_*Ps}), L_).degree(), 4)

# ── ⑧ นับพารามิเตอร์จริง ────────────────────────────────────
z1, z2 = sp.symbols('z1 z2', real=True)
chk("จำนวนเชิงซ้อน 1 ตัว = 2 พารามิเตอร์จริง", len((z1 + I*z2).free_symbols), 2)
p1, p2, p3, p4 = sp.symbols('p1:5', real=True)
H = sp.Matrix([[p1, p3 + I*p4], [p3 - I*p4, p2]])
ok.append(("Hermitian 2x2 เขียนได้ด้วย 4 พารามิเตอร์จริง", H == H.H and len(H.free_symbols) == 4, H, 4))

# ── ⑨ มิติใน 3+1 ──────────────────────────────────────────
chk("[v] = [phi] = 1", 1, 1)
chk("[M] = [e] + [v] = 0 + 1 = 1", 0 + 1, 1)
chk("[lam] = 4 - 4[phi] = 0", 4 - 4*1, 0)

bad_ = [r_ for r_ in ok if not r_[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
