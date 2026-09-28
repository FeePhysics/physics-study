"""ตรวจทุกข้ออ้างของบทที่ 14 — เฟสที่เป็นเมทริกซ์ (Yang–Mills)

convention ต่อจากบทที่ 12 (เปลี่ยน e เป็น g):
    D_mu = d_mu - i g A_mu        A_mu เป็นเมทริกซ์ 2x2
    U = 1 + i alpha (อันดับหนึ่ง)  alpha เป็นเมทริกซ์ 2x2 ขึ้นกับตำแหน่ง
    delta A_mu = (1/g) d_mu alpha + i [alpha, A_mu]      <- derive ที่นี่ ไม่ได้จำมา
    F_tx = d_t A_x - d_x A_t - i g [A_t, A_x]
"""
import sympy as sp

t, x, g = sp.symbols('t x g', real=True, nonzero=True)
eps = sp.symbols('epsilon')
I = sp.I


def mat(name):
    return sp.Matrix(2, 2, lambda i, j: sp.Function(f'{name}{i}{j}')(t, x))


def col(name):
    return sp.Matrix(2, 1, lambda i, j: sp.Function(f'{name}{i}')(t, x))


At, Ax, al = mat('At'), mat('Ax'), mat('al')
phi = col('phi')
Id = sp.eye(2)


def comm(A, B):
    return A*B - B*A


ok = []


def chkM(name, got, want):
    d = (got - want).applyfunc(lambda z: sp.simplify(sp.expand(z)))
    ok.append((name, d == sp.zeros(*d.shape), got, want))


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


def first(M):
    return M.applyfunc(lambda z: sp.expand(sp.diff(z, eps).subs(eps, 0)))


def D(f, A, c):
    return sp.diff(f, c) - I*g*A*f


# ── ① การแปลงอันดับหนึ่ง ────────────────────────────────────
U = Id + I*eps*al


def dA(A, c):
    return sp.diff(al, c)/g + I*comm(al, A)


At_T = At + eps*dA(At, t)
Ax_T = Ax + eps*dA(Ax, x)

chkM("D_t(U phi) = U (D_t phi) ถึงอันดับหนึ่ง",
     first(D(U*phi, At_T, t)), first(U*D(phi, At, t)))
chkM("D_x(U phi) = U (D_x phi) ถึงอันดับหนึ่ง",
     first(D(U*phi, Ax_T, x)), first(U*D(phi, Ax, x)))
# ตัวควบคุมฝั่งลบ: ลืมพจน์คอมมิวเทเตอร์ใน delta A ⇒ พัง
At_bad = At + eps*sp.diff(al, t)/g
bad = (first(D(U*phi, At_bad, t)) - first(U*D(phi, At, t))).applyfunc(sp.simplify)
ok.append(("ไม่มี i[alpha,A] ⇒ D ไม่โคแวเรียนต์ (พังจริง)", bad != sp.zeros(2, 1), bad, "≠0"))

# ── ② F มาจากคอมมิวเทเตอร์ของ D ─────────────────────────────
F = sp.diff(Ax, t) - sp.diff(At, x) - I*g*comm(At, Ax)
chkM("[D_t, D_x] phi = -i g F phi",
     D(D(phi, Ax, x), At, t) - D(D(phi, At, t), Ax, x), -I*g*F*phi)

F_T = sp.diff(Ax_T, t) - sp.diff(At_T, x) - I*g*comm(At_T, Ax_T)
chkM("delta F = i [alpha, F] (โคแวเรียนต์ ไม่ใช่ไม่เปลี่ยน)", first(F_T), I*comm(al, F))
chk("delta tr(F^2) = 0", sp.expand(first(F_T*F_T).trace()), 0)

# ── ③ ขีดจำกัดอาบีเลียน ─────────────────────────────────────
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2')
chkM("เมทริกซ์ทแยงคอมมิวต์กัน", comm(sp.diag(a1, a2), sp.diag(b1, b2)), sp.zeros(2, 2))
chkM("[diag(1,2), diag(3,5)] = 0", comm(sp.diag(1, 2), sp.diag(3, 5)), sp.zeros(2, 2))

# ── ④ ตัวเลขในข้อฝึก ───────────────────────────────────────
A1, B1 = sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[0, 0], [1, 0]])
chkM("AB - BA = diag(1,-1)", comm(A1, B1), sp.diag(1, -1))
chk("(AB)_11 = 1", (A1*B1)[0, 0], 1)
chk("(BA)_11 = 0", (B1*A1)[0, 0], 0)
chk("3*5 - 5*3 = 0", 3*5 - 5*3, 0)

sx = sp.Matrix([[0, 1], [1, 0]])
sy = sp.Matrix([[0, -I], [I, 0]])
sz = sp.diag(1, -1)
chk("(sx sy)_11 = i", (sx*sy)[0, 0], I)
chk("(sy sx)_11 = -i", (sy*sx)[0, 0], -I)
chk("[sx, sy]_11 = 2i", comm(sx, sy)[0, 0], 2*I)
chkM("-i [sx, sy] = 2 sz", -I*comm(sx, sy), 2*sz)
Fconst = -I*1*comm(sx, sy)
chk("F_11 ที่ A คงที่ (At=sx, Ax=sy, g=1) = 2", Fconst[0, 0], 2)
chk("F_22 ที่ A คงที่ = -2", Fconst[1, 1], -2)

chkM("[sz, sx] = [[0,2],[-2,0]]", comm(sz, sx), sp.Matrix([[0, 2], [-2, 0]]))
chk("[alpha, A]_12 = 2 (alpha=sz, A=sx)", comm(sz, sx)[0, 1], 2)
chk("[alpha, A]_21 = -2", comm(sz, sx)[1, 0], -2)
chk("g=2, d_t alpha = 6 ⇒ (1/g) d_t alpha = 3", sp.Integer(6)/2, 3)

# ── ⑤ นับกำลังของ A ─────────────────────────────────────────
lam = sp.symbols('lambda')
Fl = sp.diff(lam*Ax, t) - sp.diff(lam*At, x) - I*g*comm(lam*At, lam*Ax)
P = sp.Poly(sp.expand((Fl*Fl).trace()), lam)
chk("กำลังสูงสุดของ A ใน F = 2", sp.Poly(sp.expand(Fl[0, 0]), lam).degree(), 2)
chk("กำลังสูงสุดของ A ใน tr(F^2) = 4", P.degree(), 4)
Fab = sp.diff(lam*Ax, t) - sp.diff(lam*At, x)
chk("อาบีเลียน: กำลังสูงสุดของ A ใน F^2 = 2",
    sp.Poly(sp.expand((Fab*Fab).trace()), lam).degree(), 2)
powers = sorted(set(m[0] for m in P.monoms()))
chk("tr(F^2) มีกำลัง 2, 3, 4 ⇒ พจน์คุยกับตัวเอง 2 แบบ", len([p for p in powers if p > 2]), 2)


# ── ⑥ นับจำนวนสนามแรง ───────────────────────────────────────
def n_params(n):
    # Hermitian n×n: ทแยง n ช่อง (จริง) + นอกทแยง n(n-1)/2 ช่อง (ซับซ้อน = 2 จริง) = n^2 · ไร้เทรซ ตัด 1
    return n + 2*(n*(n - 1)//2) - 1


chk("SU(2): 3", n_params(2), 3)
chk("SU(3): 8", n_params(3), 8)
chk("Hermitian 2×2 มี 4 พารามิเตอร์จริงก่อนตัดเทรซ", 2 + 2*1, 4)

# ── ⑦ มิติ ────────────────────────────────────────────────
chk("[g] = [d] - [A] = 0", 1 - 1, 0)

bad_ = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
