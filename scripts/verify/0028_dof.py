"""ตรวจทุกข้ออ้างของบทที่ 16 — นับองศาอิสระ: ของหายไปไหนเมื่อแสงมีมวล

convention ต่อจากบท 12-15: D = d - i e A · signature (-,+) · 1+1 มิติ
phi = (v + h) e^{i theta}   (h = การสั่นตามรัศมี · theta = การสั่นตามวงกลม)
⚠️ มวลทุกตัววัดจากความสัมพันธ์การกระจายของคลื่นระนาบ ไม่ได้อ่านจากรูปพจน์
"""
import sympy as sp

t, x = sp.symbols('t x', real=True)
e, v, lam, k = sp.symbols('e v lam k', positive=True)
I = sp.I
ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


# ── ① มิติบน R ─────────────────────────────────────────────
a, b = sp.symbols('a b', real=True)
# C: ฐาน {1, i} — อิสระเชิงเส้นบน R และแผ่ทั่ว
chk("3 + 4i = 3*1 + 4*i (สองสัมประสิทธิ์จริง)", 3*1 + 4*I, 3 + 4*I)
M = sp.Matrix([[1, 0], [0, 1]])   # เวกเตอร์ของ 1 กับ i ในพิกัด (Re, Im)
chk("{1, i} อิสระเชิงเส้นบน R (det ≠ 0) ⇒ dim C = 2", M.det(), 1)

one = sp.eye(2)
sx = sp.Matrix([[0, 1], [1, 0]])
sy = sp.Matrix([[0, -I], [I, 0]])
sz = sp.diag(1, -1)
basis = [one, sx, sy, sz]
c = sp.symbols('c0:4', real=True)
p1, p2, p3, p4 = sp.symbols('p1:5', real=True)
H = sp.Matrix([[p1, p3 + I*p4], [p3 - I*p4, p2]])
combo = sum((ci*Bi for ci, Bi in zip(c, basis)), sp.zeros(2, 2))
sol = sp.solve(list(combo - H), c, dict=True)
ok.append(("Hermitian ทุกตัว = c0*1 + c1*sx + c2*sy + c3*sz ด้วย c จริง (แผ่ทั่ว)",
           len(sol) == 1 and all(sp.simplify(sp.im(sol[0][ci])) == 0 for ci in c), sol, "unique real"))
Z = sp.solve(list(combo), c, dict=True)
ok.append(("ฐานสี่ตัวอิสระเชิงเส้น (มีแต่ c = 0 ที่ให้เมทริกซ์ศูนย์)",
           Z == [{ci: 0 for ci in c}], Z, "c=0"))
chk("เทรซของ sx, sy, sz = 0 ⇒ ไร้เทรซเหลือ 3", sx.trace() + sy.trace() + sz.trace(), 0)
chk("เทรซของ 1 = 2 (ตัวที่ถูกตัด)", one.trace(), 2)

# รูปทั่วไปในบทเรียน [[a, b+ic],[b-ic, d]] เป็น Hermitian จริง และเทรซศูนย์ ⇒ d = -a เหลือ 3 ตัว
ga, gb, gc, gd = sp.symbols('ga gb gc gd', real=True)
Hg = sp.Matrix([[ga, gb + I*gc], [gb - I*gc, gd]])
ok.append(("[[a, b+ic],[b-ic, d]] = ทรานสโพสสังยุคของตัวเอง", sp.simplify(Hg - Hg.H) == sp.zeros(2, 2), Hg.H, Hg))
chk("เทรซศูนย์ ⇒ d = -a", sp.solve(sp.Eq(Hg.trace(), 0), gd)[0], -ga)

# ── ② phi = (v + h) e^{i theta} ─────────────────────────────
h = sp.Function('h')(t, x)
th = sp.Function('theta')(t, x)
phi = (v + h)*sp.exp(I*th)
phs = (v + h)*sp.exp(-I*th)
chk("phis*phi = (v + h)^2 — ไม่มี theta", sp.simplify(phs*phi), (v + h)**2)
V = lam*(phs*phi - v**2)**2
chk("dV/dtheta = 0 ⇒ theta ไม่มีศักย์", sp.diff(sp.simplify(V), th), 0)
hs = sp.symbols('hs')
Vh = sp.expand(lam*((v + hs)**2 - v**2)**2)
chk("สัมประสิทธิ์ h^2 ใน V = 4 lam v^2", Vh.coeff(hs, 2), 4*lam*v**2)
chk("V ไม่มีพจน์ h^1 (v เป็นจุดต่ำสุดจริง)", Vh.coeff(hs, 1), 0)

# จลน์ (ไม่มี A): |d_t phi|^2 = h_t^2 + (v+h)^2 theta_t^2
kin_t = sp.simplify(sp.expand(sp.diff(phs, t)*sp.diff(phi, t)))
chk("|d_t phi|^2 = h_t^2 + (v+h)^2 theta_t^2",
    kin_t, sp.diff(h, t)**2 + (v + h)**2*sp.diff(th, t)**2)

# ── ③ มวลของ h และ theta จากคลื่นระนาบ ──────────────────────
w = sp.symbols('omega', positive=True)
q = sp.symbols('q')


def dispersion(L2, f):
    """L2 = ลากรองจ์กำลังสองของสนาม f · คืนเงื่อนไขของ omega สำหรับ f = cos(kx - wt)"""
    ft, fx = sp.symbols('ft fx')
    Ls = L2.subs({sp.diff(f, t): ft, sp.diff(f, x): fx}).subs(f, q)
    EL = (sp.diff(sp.diff(Ls, ft).subs({ft: sp.diff(f, t), fx: sp.diff(f, x), q: f}), t)
          + sp.diff(sp.diff(Ls, fx).subs({ft: sp.diff(f, t), fx: sp.diff(f, x), q: f}), x)
          - sp.diff(Ls, q).subs({ft: sp.diff(f, t), fx: sp.diff(f, x), q: f}))
    wave = sp.cos(k*x - w*t)
    return sp.simplify(EL.subs(f, wave).doit()/wave)


Lh = sp.diff(h, t)**2 - sp.diff(h, x)**2 - 4*lam*v**2*h**2
Lth = v**2*(sp.diff(th, t)**2 - sp.diff(th, x)**2)
dh = dispersion(Lh, h)
dth = dispersion(Lth, th)
chk("h: เงื่อนไขคือ omega^2 = k^2 + 4 lam v^2",
    sp.solve(dh, w)[0]**2, k**2 + 4*lam*v**2)
chk("theta: เงื่อนไขคือ omega^2 = k^2 (ไร้มวล = โกลด์สโตน)",
    sp.solve(dth, w)[0]**2, k**2)

# ── ④ theta ถูก A กิน: เขียน L ด้วย B = A - theta'/e แล้ว theta หายหมด ──
At = sp.Function('At')(t, x)
Ax = sp.Function('Ax')(t, x)
Dt = sp.diff(phi, t) - I*e*At*phi
Dts = sp.diff(phs, t) + I*e*At*phs
got = sp.simplify(sp.expand(Dts*Dt))
Bt = At - sp.diff(th, t)/e
chk("|D_t phi|^2 = h_t^2 + e^2 (v+h)^2 (A_t - theta_t/e)^2",
    got, sp.diff(h, t)**2 + e**2*(v + h)**2*Bt**2)
# gauge ที่ alpha = -theta ทำให้ phi จริง
alpha = -th
chk("เฟสหลังแปลง = theta + alpha = 0", th + alpha, 0)
chk("A หลังแปลง = A + (1/e) d alpha = A - theta'/e", At + sp.diff(alpha, t)/e, Bt)
# F ของ B เท่ากับ F ของ A (theta ไม่เข้าไปใน F)
Bx = Ax - sp.diff(th, x)/e
chk("F(B) = F(A)", (sp.diff(Bx, t) - sp.diff(Bt, x)) - (sp.diff(Ax, t) - sp.diff(At, x)), 0)

# ── ⑤ แสงไร้มวลใน 1+1 ไม่มีคลื่น ────────────────────────────
# ไม่มีแหล่ง: d_x F = 0 และ d_t F = 0 (บทที่ 13 ที่ J = 0) ⇒ F คงที่
Fconst = sp.Symbol('F0')
chk("F = ค่าคงที่ เป็นไปตาม d_x F = 0", sp.diff(Fconst, x), 0)
chk("F = ค่าคงที่ เป็นไปตาม d_t F = 0", sp.diff(Fconst, t), 0)
wave_F = sp.cos(k*x - w*t)
bad = sp.simplify(sp.diff(wave_F, x))
ok.append(("คลื่น cos(kx - wt) ที่ k ≠ 0 ไม่เป็นไปตาม d_x F = 0 (ตัวควบคุมฝั่งลบ)", bad != 0, bad, "≠0"))

# ── ⑥ นับก่อน/หลัง ใน 1+1 ─────────────────────────────────
before = 2 + 0          # phi เชิงซ้อน (2 คลื่น) + แสงไร้มวล (0)
after = 1 + 1 + 0       # h (1) + แสงมีมวล A_x (1 · บทที่ 15) + theta (0 — ถูกกิน)
chk("ก่อน = 2", before, 2)
chk("หลัง = 2", after, 2)
chk("ก่อน - หลัง = 0 (ของไม่หาย มันย้ายที่)", before - after, 0)

hh, tt = sp.symbols('hh tt', real=True)
chk("หมุนด้วย alpha = -theta ⇒ phi = v + h (จริง)",
    sp.simplify((v + hh)*sp.exp(I*tt)*sp.exp(I*(-tt))), v + hh)

# ── ⑦ ตัวเลขในข้อฝึก ───────────────────────────────────────
chk("lam=2, v=3 ⇒ m_h^2 = 72", (4*lam*v**2).subs({lam: 2, v: 3}), 72)
chk("lam=1, v=1 ⇒ m_h^2 = 4", (4*lam*v**2).subs({lam: 1, v: 1}), 4)
chk("k = 0, m_h^2 = 4 ⇒ omega = 2", sp.sqrt(0 + 4), 2)
chk("theta = 0.3 ⇒ alpha = -0.3", -sp.Rational(3, 10), -sp.Rational(3, 10))
chk("[m_h^2] = [lam] + 2[v] = 0 + 2", 0 + 2*1, 2)

bad_ = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
