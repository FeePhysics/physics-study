#!/usr/bin/env python
"""ตรวจขั้นพีชคณิตของคู่มืออ่าน Dodelson บทที่ 5 (reference/dodelson-ch05.html)

รัน:  python scripts/verify/dodelson_ch05.py
ทำทุกอย่างใน "real space" (หนังสือทำใน Fourier) แล้วเทียบด้วยกฎ
    k_i -> -i ∂_i  ดังนั้น  k^2 X -> -∇²X  และ  k_i k_j X -> -∂_i∂_j X
นิยาม Ricci ตามหนังสือ (2.31)/(5.8):
    R_{mu nu} = Γ^a_{mu nu,a} - Γ^a_{mu a,nu} + Γ^a_{b a}Γ^b_{mu nu} - Γ^a_{b nu}Γ^b_{mu a}
"""
import sympy as sp

ok = lambda name: print(f"OK  {name}")
eps = sp.symbols('epsilon')
t, x, y, z = sp.symbols('t x y z')
X = [t, x, y, z]
a = sp.Function('a')(t)
H = sp.diff(a, t) / a
add = sp.diff(a, t, 2)


def lin(e):
    """ตัดเหลืออันดับศูนย์ + อันดับหนึ่งใน eps"""
    e0 = e.subs(eps, 0)
    e1 = sp.diff(e, eps).subs(eps, 0)
    return sp.expand(e0 + eps * e1)


def geometry(g):
    ginv = g.inv()
    ginv = ginv.applyfunc(lambda e: lin(sp.series(e, eps, 0, 2).removeO()))
    Gam = [[[lin(sp.Rational(1, 2) * sum(ginv[m, n] * (sp.diff(g[n, al], X[be]) + sp.diff(g[n, be], X[al])
                                                         - sp.diff(g[al, be], X[n])) for n in range(4)))
             for be in range(4)] for al in range(4)] for m in range(4)]

    def Ric(mu, nu):
        r = 0
        for A in range(4):
            r += sp.diff(Gam[A][mu][nu], X[A]) - sp.diff(Gam[A][mu][A], X[nu])
            for B in range(4):
                r += Gam[A][B][A] * Gam[B][mu][nu] - Gam[A][B][nu] * Gam[B][mu][A]
        return lin(r)

    R = sp.Matrix(4, 4, lambda i, j: Ric(i, j))
    Rs = lin(sum(ginv[i, j] * R[i, j] for i in range(4) for j in range(4)))
    return ginv, Gam, R, Rs


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


# ===================================================== scalar perturbations, metric (4.9)
Psi = sp.Function('Psi')(t, x, y, z)
Phi = sp.Function('Phi')(t, x, y, z)
g = sp.diag(-(1 + 2 * eps * Psi), *[a**2 * (1 + 2 * eps * Phi)] * 3)
ginv, Gam, R, Rs = geometry(g)
d = lambda f, *v: sp.diff(f, *v)
lap = lambda f: d(f, x, x) + d(f, y, y) + d(f, z, z)
sp_ = [x, y, z]

# (5.3)-(5.7) Christoffel
assert zero(Gam[0][0][0] - eps * d(Psi, t))
for i in range(3):
    assert zero(Gam[0][0][i + 1] - eps * d(Psi, sp_[i]))
    assert zero(Gam[i + 1][0][0] - eps * d(Psi, sp_[i]) / a**2)
    for j in range(3):
        dij = 1 if i == j else 0
        assert zero(Gam[0][i + 1][j + 1] - dij * a**2 * (H + 2 * eps * H * (Phi - Psi) + eps * d(Phi, t)))
        assert zero(Gam[i + 1][j + 1][0] - dij * (H + eps * d(Phi, t)))
        for k in range(3):
            dik = 1 if i == k else 0
            djk = 1 if j == k else 0
            want = eps * (dij * d(Phi, sp_[k]) + dik * d(Phi, sp_[j]) - djk * d(Phi, sp_[i]))
            assert zero(Gam[i + 1][j + 1][k + 1] - want)
ok("(5.3)-(5.7) Christoffel symbols ของ metric (4.9) อันดับหนึ่ง")

# (5.13) R_00
R00 = -3 * add / a + eps * (lap(Psi) / a**2 - 3 * d(Phi, t, t) + 3 * H * (d(Psi, t) - 2 * d(Phi, t)))
assert zero(R[0, 0] - R00)
ok("(5.13) R_00  [k^2 Psi -> -∇²Psi]")

# (5.14) R_ij
for i in range(3):
    for j in range(3):
        dij = 1 if i == j else 0
        want = dij * ((2 * a**2 * H**2 + a * add) * (1 + 2 * eps * Phi - 2 * eps * Psi)
                      + eps * (a**2 * H * (6 * d(Phi, t) - d(Psi, t)) + a**2 * d(Phi, t, t) - lap(Phi))) \
            - eps * d(Phi + Psi, sp_[i], sp_[j])
        assert zero(R[i + 1, j + 1] - want)
ok("(5.14) R_ij  [k_i k_j X -> -∂_i∂_j X]")

# (5.15)-(5.17) Ricci scalar
Rs_want = 6 * (H**2 + add / a) + eps * (-12 * Psi * (H**2 + add / a) - 2 * lap(Psi) / a**2 + 6 * d(Phi, t, t)
                                        - 6 * H * (d(Psi, t) - 4 * d(Phi, t)) - 4 * lap(Phi) / a**2)
assert zero(Rs - Rs_want)
ok("(5.15)-(5.17) R = 6(H^2 + a''/a) + deltaR")

# (5.19)-(5.21) G^0_0
G00 = lin(ginv[0, 0] * (R[0, 0] - g[0, 0] * Rs / 2))
assert zero(G00 - (-3 * H**2 + eps * (-6 * H * d(Phi, t) + 6 * Psi * H**2 + 2 * lap(Phi) / a**2)))
ok("(5.21) delta G^0_0 = -6H Phi,0 + 6 Psi H^2 - 2k^2 Phi/a^2  (อันดับศูนย์ -3H^2 = Friedmann)")

# (5.28)-(5.29) G^i_j = A delta_ij - ∂_i∂_j(Phi+Psi)/a^2
Gij = sp.Matrix(3, 3, lambda i, j: lin(sum(ginv[i + 1, k + 1] * (R[k + 1, j + 1] - g[k + 1, j + 1] * Rs / 2)
                                           for k in range(3))))
Mrest = sp.Matrix(3, 3, lambda i, j: sp.expand(Gij[i, j] + eps * d(Phi + Psi, sp_[i], sp_[j]) / a**2))
for i in range(3):
    for j in range(3):
        if i != j:
            assert zero(Mrest[i, j])
assert zero(Mrest[0, 0] - Mrest[1, 1]) and zero(Mrest[1, 1] - Mrest[2, 2])
ok("(5.29) G^i_j = A δ_ij + k_i k_j(Phi+Psi)/a^2 — ส่วนที่เหลือเป็นสัดส่วนกับ δ_ij จริง")

# (5.30) projection
k1, k2, k3, A_, Xs = sp.symbols('k1 k2 k3 A X')
kv = [k1, k2, k3]
kk = sp.sqrt(k1**2 + k2**2 + k3**2)
proj = sum((kv[i] * kv[j] / kk**2 - sp.Rational(1, 3) * (1 if i == j else 0))
           * (A_ * (1 if i == j else 0) + kv[i] * kv[j] * Xs / a**2) for i in range(3) for j in range(3))
assert zero(proj - sp.Rational(2, 3) * kk**2 * Xs / a**2)
ok("(5.30) (k̂_i k̂^j - δ/3) G^i_j = 2k^2(Phi+Psi)/(3a^2) และฆ่าพจน์ A δ_ij")

# (5.26)->(5.27) เปลี่ยนเป็น conformal time
eta = sp.symbols('eta')
ac = sp.Function('a')(eta)
Ph = sp.Function('Phi')(eta)
Ps, kk_, src = sp.symbols('Psi k S')
Hc = sp.diff(ac, eta) / ac**2                  # H = (da/dt)/a = a'/a^2
Phi_t = sp.diff(Ph, eta) / ac                  # ∂_t = (1/a)∂_eta
e526 = -3 * Hc * Phi_t + 3 * Ps * Hc**2 - kk_**2 * Ph / ac**2 + 4 * sp.pi * src   # = 0
e527 = kk_**2 * Ph + 3 * sp.diff(ac, eta) / ac * (sp.diff(Ph, eta) - Ps * sp.diff(ac, eta) / ac) \
    - 4 * sp.pi * ac**2 * src
assert zero(e526 * (-ac**2) - e527)
ok("(5.26)->(5.27) คูณ -a^2 และ ∂_t = (1/a)∂_eta")

# (5.24) photon T^0_0 และ (5.32) anisotropic stress
p, mu = sp.symbols('p mu', positive=True)
fb = 1 / (sp.exp(p) - 1)                       # T = 1
rho = (2 * sp.Integral(4 * sp.pi * p**2 * p * fb, (p, 0, sp.oo)) / (2 * sp.pi)**3).evalf(25)
pert = (2 * sp.Integral(4 * sp.pi * p**2 * p * (-p * sp.diff(fb, p)), (p, 0, sp.oo)) / (2 * sp.pi)**3).evalf(25)
assert abs(pert - 4 * rho) < 1e-18
ok("(5.23)-(5.24) ชิ้น -p ∂f/∂p Theta_0 ให้ 4 rho Theta_0")
P2 = sp.legendre(2, mu)
assert zero((mu**2 - sp.Rational(1, 3)) - sp.Rational(2, 3) * P2)
I4 = sp.Integral(p**4 * sp.diff(fb, p), (p, 0, sp.oo)).evalf(25) / (2 * sp.pi**2)
# ∫dmu/2 P2 Theta = -Theta_2 (จาก 4.99, l=2)
stress_over_Theta2 = -2 * I4 * sp.Rational(2, 3) * (-1)
assert abs(stress_over_Theta2 - (-sp.Rational(8, 3) * rho)) < 1e-18
ok("(5.31)-(5.32) mu^2-1/3 = (2/3)P2 และ stress ของ photon = -8 rho Theta_2/3")
# (5.33)
G, rg, rn, Th2, N2, Sum = sp.symbols('G rho_g rho_n Theta2 N2 Sum')
lhs = sp.Rational(2, 3) * kk_**2 * Sum / ac**2          # Sum = Phi + Psi
rhs = 8 * sp.pi * G * (-sp.Rational(8, 3)) * (rg * Th2 + rn * N2)
sol = sp.solve(sp.Eq(lhs, rhs), Sum)[0]
assert zero(kk_**2 * sol + 32 * sp.pi * G * ac**2 * (rg * Th2 + rn * N2))
ok("(5.33) k^2(Phi+Psi) = -32 pi G a^2 [rho_g Theta2 + rho_nu N2]")

# (5.81) จาก (5.27) + (5.84)
Pm_dm, Pr0, Pmv, Pr1, aH, Phd = sp.symbols('rhom_deltam rhor_Thetar0 rhom_vm rhor_Thetar1 aH Phidot')
Phi_s = sp.symbols('Phi_s')
# (5.84): Phidot - aH Psi = 4πGa²/(ik) [rho_m v_m - 4i rho_r Theta_r1]
Phd_sol = aH * Ps + 4 * sp.pi * G * ac**2 / (sp.I * kk_) * (Pmv - 4 * sp.I * Pr1)
k2Phi = 4 * sp.pi * G * ac**2 * (Pm_dm + 4 * Pr0) - 3 * aH * (Phd_sol - Ps * aH)
want = 4 * sp.pi * G * ac**2 * (Pm_dm + 4 * Pr0 + 3 * aH / kk_ * (sp.I * Pmv + 4 * Pr1))
assert zero(k2Phi - want)
ok("(5.81) สมการพีชคณิตของ Phi จาก (5.27) + (5.84)")

# ===================================================== tensor perturbations, metric (5.34)
hp = sp.Function('h_p')(t, z)
hx = sp.Function('h_x')(t, z)
gT = sp.Matrix([[-1, 0, 0, 0],
                [0, a**2 * (1 + eps * hp), a**2 * eps * hx, 0],
                [0, a**2 * eps * hx, a**2 * (1 - eps * hp), 0],
                [0, 0, 0, a**2]])
ginvT, GamT, RT, RsT = geometry(gT)
Hm = sp.Matrix([[hp, hx, 0], [hx, -hp, 0], [0, 0, 0]])     # (5.37)
assert all(zero(GamT[0][0][m]) for m in range(4))                         # (5.35)
for i in range(3):
    for j in range(3):
        gij = gT[i + 1, j + 1]
        assert zero(GamT[0][i + 1][j + 1] - lin(H * gij + a**2 * eps * d(Hm[i, j], t) / 2))   # (5.39)
        assert zero(GamT[i + 1][0][j + 1] - ((1 if i == j else 0) * H + eps * d(Hm[i, j], t) / 2))  # (5.42)
        for k in range(3):
            # (5.43) ใน real space: i k_k H_ij -> ∂_k H_ij
            want = eps * sp.Rational(1, 2) * (d(Hm[i, j], sp_[k]) + d(Hm[i, k], sp_[j]) - d(Hm[j, k], sp_[i]))
            assert zero(GamT[i + 1][j + 1][k + 1] - want)
ok("(5.35)-(5.43) Christoffel ของ tensor perturbation")
assert zero(RT[0, 0] + 3 * add / a)
ok("(5.47) R_00 = -3 a''/a — ไม่มี h ในอันดับหนึ่ง")
for i in range(3):
    for j in range(3):
        gij = gT[i + 1, j + 1]
        want = lin(gij * (add / a + 2 * H**2)) + eps * (sp.Rational(3, 2) * a**2 * H * d(Hm[i, j], t)
                                                       + a**2 * d(Hm[i, j], t, t) / 2 - d(Hm[i, j], z, z) / 2)
        assert zero(RT[i + 1, j + 1] - want)
ok("(5.57) R_ij ของ tensor  [k^2 H_ij -> -∂_z² H_ij]")
assert zero(RsT - 6 * (H**2 + add / a))
ok("(5.58) Ricci scalar ไม่มีพจน์ tensor อันดับหนึ่ง")
GT = lambda i, j: lin(sum(ginvT[i, kq] * (RT[kq, j] - gT[kq, j] * RsT / 2) for kq in range(4)))
E = lambda f: sp.Rational(3, 2) * H * d(f, t) + d(f, t, t) / 2 - d(f, z, z) / (2 * a**2)
assert zero(GT(1, 1) - GT(2, 2) - eps * 2 * E(hp))
assert zero(GT(1, 2) - eps * E(hx))
ok("(5.60)-(5.61) δG^1_1-δG^2_2 = 3H h+,0 + h+,00 + k^2 h+/a^2 · Exercise 9: G^1_2 ให้สมการเดียวกันกับ hx")
assert zero(GT(0, 0) + 3 * H**2)
ok("(5.65) tensor ไม่เข้า G^0_0")
# (5.66) projection ของส่วน tensor: k̂ = z ⇒ ใช้เฉพาะ G^3_3 กับ trace
tr = sum(GT(i, i) for i in (1, 2, 3))
assert zero(sp.diff(GT(3, 3) - tr / 3, eps))           # ส่วนอันดับหนึ่งเป็นศูนย์
ok("(5.66) (k̂_ik̂_j - δ_ij/3) δG^i_j ไม่มีพจน์ tensor")

# (5.62)-(5.63) conformal time
hc = sp.Function('h')(eta)
h_t = sp.diff(hc, eta) / ac
h_tt = sp.diff(h_t, eta) / ac
lhs562 = ac**2 * (3 * Hc * h_t + h_tt + kk_**2 * hc / ac**2)
assert zero(lhs562 - (sp.diff(hc, eta, 2) + 2 * sp.diff(ac, eta) / ac * sp.diff(hc, eta) + kk_**2 * hc))
ok("(5.61)->(5.62) แปลงเป็น conformal time")

# (5.63) ไม่มีการขยายตัว: e^{±ik eta}
for s in (1, -1):
    hh = sp.exp(s * sp.I * kk_ * eta)
    assert zero(sp.diff(hh, eta, 2) + kk_**2 * hh)
# matter era a ∝ eta^2 ⇒ 2a'/a = 4/eta · คำตอบ 3 j1(x)/x (ตัวส่วนใน 5.87)
xx = sp.symbols('x', positive=True)
j1 = sp.sin(xx) / xx**2 - sp.cos(xx) / xx
hm = 3 * j1 / xx
assert zero(sp.diff(hm, xx, 2) + 4 / xx * sp.diff(hm, xx) + hm)
assert sp.limit(hm, xx, 0) == 1
# radiation era a ∝ eta ⇒ 2a'/a = 2/eta · คำตอบ j0 = sin x / x  (นอกหนังสือ: หนังสือไม่ได้เฉลย Ex.12)
hr = sp.sin(xx) / xx
assert zero(sp.diff(hr, xx, 2) + 2 / xx * sp.diff(hr, xx) + hr)
ok("(5.63)-(5.64) e^{±ik eta} · matter era 3j1(k eta)/(k eta) · radiation era j0(k eta)")

# ===================================================== gauge (5.67)-(5.76), Exercise 14
A_f, B_f, ps_f, E_f, xi0, xi = [sp.Function(n)(t, x, y, z) for n in ('A', 'B', 'psi', 'E', 'xi0', 'xi')]


def gmetric(A, B, ps, E, aa):
    gg = sp.zeros(4, 4)
    gg[0, 0] = -(1 + 2 * eps * A)
    for i in range(3):
        gg[0, i + 1] = gg[i + 1, 0] = -aa * eps * d(B, sp_[i])
        for j in range(3):
            gg[i + 1, j + 1] = aa**2 * ((1 if i == j else 0) * (1 + 2 * eps * ps) - 2 * eps * d(E, sp_[i], sp_[j]))
    return gg


# dot = conformal derivative = a ∂_t
dot = lambda f: a * d(f, t)
At = A_f - dot(xi0) / a                      # (5.74)
pst = ps_f - H * xi0                         # (5.75)
Bt = B_f - xi0 / a + dot(xi)
Et = E_f + xi
a_shift = a + eps * xi0 * d(a, t)            # a(t + xi0)
gt = gmetric(At, Bt, pst, Et, a_shift)
xt = [t + eps * xi0] + [sp_[i] + eps * d(xi, sp_[i]) for i in range(3)]   # (5.70)
J = sp.Matrix(4, 4, lambda al, m: d(xt[al], X[m]))
lhs = (J.T * gt * J).applyfunc(lin)
rhs = gmetric(A_f, B_f, ps_f, E_f, a).applyfunc(lin)
assert all(zero(lhs[i, j] - rhs[i, j]) for i in range(4) for j in range(4))
ok("(5.69)-(5.75) กฎแปลง A, B, psi, E ภายใต้ (5.70)  [จุด = d/d eta = a d/dt]")
PhiA = lambda A, B, E: A + d(a * (dot(E) - B), t)       # (1/a)∂_eta = ∂_t
PhiH = lambda ps, B, E: -ps + a * H * (B - dot(E))
assert zero(PhiA(At, Bt, Et) - PhiA(A_f, B_f, E_f))
assert zero(PhiH(pst, Bt, Et) - PhiH(ps_f, B_f, E_f))
ok("(5.76) Phi_A และ Phi_H ไม่เปลี่ยนภายใต้การแปลง gauge")

print("\nผ่านทั้งหมด — dodelson_ch05")
