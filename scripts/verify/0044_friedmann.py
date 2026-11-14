"""บทที่ 29 · สมการ Friedmann ⇒ a ∝ t^(2/3) — ตรวจทุกขั้นและทุกเฉลยของ lessons/0043-friedmann.html

รัน: python scripts/verify/0044_friedmann.py   (ทุกบรรทัดต้องขึ้น ✅ · assert ล้มเมื่อผิด)
ตรวจสองทางเสมอ: (ก) จับคู่กำลังของ t แบบที่บทสอน  (ข) แก้สมการเชิงอนุพันธ์ตรง ๆ ด้วย sympy
ค่าคงตัว (Mpc, G) ใช้ค่าที่หนังสือให้ (Dodelson หน้า 5) · มวล proton เป็นค่านอกหนังสือ
"""
import sys
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


a, x, N, m, t, C, p, K, n = sp.symbols('a x N m t C p K n', positive=True)
G, H0, rho0, rhocr = sp.symbols('G H0 rho0 rho_cr', positive=True)


def power_of(expr, var):
    """กำลังของ var ใน expr ที่เป็นรูป (ค่าคงตัว)·var^k — ใช้ log-derivative: k = var·d(ln expr)/d var"""
    return sp.simplify(var * sp.diff(sp.log(expr), var))


# ── ① (1.2) + (1.3) ───────────────────────────────────────
rho = sp.Symbol('rho', positive=True)
fried = lambda rho_, a_, rho0_: 8 * sp.pi * G / 3 * (rho_ + (rhocr - rho0_) / a_**2)
check("(1.2) ที่ a=1 ในวงเล็บ = ρcr", sp.simplify((rho0 + (rhocr - rho0) / 1) - rhocr) == 0)
check("(1.2) ที่ a=1 + (1.3) ⇒ H² = H0²",
      sp.simplify(fried(rho0, 1, rho0).subs(rhocr, 3 * H0**2 / (8 * sp.pi * G)) - H0**2) == 0)
check("ρ0 = ρcr ⇒ (1.2) เหลือ H² = 8πGρ/3 (flat)",
      sp.simplify(fried(rho, a, rhocr) - 8 * sp.pi * G * rho / 3) == 0)
# (1.3) คือ flat-Friedmann ที่วันนี้ แก้หา ρ0
sol = sp.solve(sp.Eq(H0**2, 8 * sp.pi * G * rho0 / 3), rho0)
check("flat ที่วันนี้ แก้หา ρ0 ได้ 3H0²/(8πG) = (1.3)", sol == [3 * H0**2 / (8 * sp.pi * G)])
# ข้อ 5 (ตัวเลขสมมุติ ρcr = 10, ρ0 = 3)
br = lambda aa: 3 * aa**-3 + (10 - 3) / aa**2     # matter ∝ a^-3 เพื่อดูภาพเท่านั้น
check("ข้อ 5 ขั้น 1: วงเล็บที่ a=1 = 10", sp.nsimplify(3 + (10 - 3) / 1) == 10)
check("ข้อ 5 ขั้น 2: พจน์ความโค้งที่ a=1/2 = 7/(1/4) = 28", sp.Rational(7) / sp.Rational(1, 2)**2 == 28)
check("ข้อ 5 ขั้น 3: flat ⇒ ρ0 = ρcr = 10 ⇒ พจน์ความโค้ง = 0", (10 - 10) == 0)

# ── ② ρ ∝ a^-3 (matter) ───────────────────────────────────
V = (a * x)**3                       # ปริมาตรจริงของกล่องที่มีด้าน comoving x
rho_m = N * m / V
check("ข้อ 1 ขั้น 1: ปริมาตรจริง ∝ a^3", power_of(V, a) == 3)
check("ข้อ 1 ขั้น 2: ρ_m ∝ a^-3", power_of(rho_m, a) == -3)
check("ข้อ 1 ขั้น 3: a=1/2 ⇒ ρ_m = 8 เท่าของวันนี้",
      sp.simplify(rho_m.subs(a, sp.Rational(1, 2)) / rho_m.subs(a, 1)) == 8)

# ── ③ ρ ∝ a^-4 (radiation) · T ∝ 1/a ───────────────────────
lam_emit, hP, c = sp.symbols('lambda_emit h c', positive=True)
lam = lam_emit * a                  # จากบทที่ 28: λ ∝ a (1+z = 1/a)
E1 = hP * c / lam                   # พลังงานต่อโฟตอน (E = hc/λ · ฟิสิกส์ ม.ปลาย)
n_ph = N / V
rho_r = E1 * n_ph
check("ข้อ 2 ขั้น 1: พลังงานต่อโฟตอน ∝ a^-1", power_of(E1, a) == -1)
check("ข้อ 2 ขั้น 2: จำนวนโฟตอนต่อปริมาตร ∝ a^-3", power_of(n_ph, a) == -3)
check("ข้อ 2 ขั้น 3: ρ_r ∝ a^-4", power_of(rho_r, a) == -4)
check("ข้อ 2 ขั้น 4: T = T0/a ที่ a = 1/1000 ⇒ 2725 K (T0 = 2.725 K)",
      abs(2.725 / (1 / 1000) - 2725) < 1e-9)

# ── ④ a ∝ t^p — (ก) จับคู่กำลัง ───────────────────────────
aa = C * t**p
H = sp.simplify(sp.diff(aa, t) / aa)
check("a = C t^p ⇒ H ∝ t^-1", power_of(H, t) == -1)
check("ข้อ 3 ขั้น 1–2: H² ∝ t^-2", power_of(H**2, t) == -2)
check("ข้อ 3 ขั้น 3: ρ ∝ a^-3 = C^-3 t^(-3p) ⇒ กำลังของ t = -3p", sp.simplify(power_of(aa**-3, t) + 3 * p) == 0)
check("ข้อ 3 ขั้น 4: -2 = -3p ⇒ p = 2/3", sp.solve(sp.Eq(-2, -3 * p), p) == [sp.Rational(2, 3)])
check("ข้อ 4: -2 = -4p ⇒ p = 1/2", sp.solve(sp.Eq(-2, -4 * p), p) == [sp.Rational(1, 2)])
check("ข้อ 4 ขั้น 3: p=1/2 ⇒ H t = 1/2", sp.simplify((H * t).subs(p, sp.Rational(1, 2))) == sp.Rational(1, 2))
# กรณีทั่วไป ρ ∝ a^-n ⇒ p = 2/n
check("ทั่วไป ρ ∝ a^-n ⇒ p = 2/n", sp.solve(sp.Eq(-2, -n * p), p) == [2 / n])

# ── ④ (ข) แก้สมการเชิงอนุพันธ์ตรง ๆ: ȧ² = K a^(2-n) (คือ H² = K a^-n) ────────
af = sp.Function('a')
for nn, pp in [(3, sp.Rational(2, 3)), (4, sp.Rational(1, 2))]:
    # ȧ = sqrt(K) a^(1-n/2) — กิ่งที่ขยายตัว
    ode = sp.Eq(af(t).diff(t), sp.sqrt(K) * af(t)**(1 - sp.Rational(nn, 2)))
    sols = sp.dsolve(ode)
    sols = sols if isinstance(sols, list) else [sols]
    # หยิบคำตอบที่เป็นจริงบวก แล้วตั้ง a(0) = 0 (ค่าคงตัวอินทิเกรตเป็น 0)
    good = []
    for s in sols:
        r = s.rhs.subs('C1', 0)
        r = sp.simplify(r)
        if r.subs({t: 1, K: 1}).is_positive:
            good.append(r)
    check(f"dsolve (n={nn}): a(0)=0 ⇒ a ∝ t^{pp}", any(power_of(r, t) == pp for r in good))
    # ตรวจกลับ: แทนคำตอบกำลังลงไปในสมการเดิมได้ศูนย์ (มีค่าคงตัว C ที่ทำให้จริง)
    A = sp.Symbol('A', positive=True)
    trial = A * t**pp
    resid = sp.simplify(sp.diff(trial, t)**2 - K * trial**(2 - nn))
    Asol = sp.solve(resid, A)
    check(f"แทน a = A t^{pp} กลับ: มี A > 0 ที่ทำให้ ȧ² = K a^(2-{nn})", len([s for s in Asol if s.is_positive]) == 1)

# ── ⑤ ตัวเลข (1.5) (1.6) · h = 1 ───────────────────────────
Mpc = 3.0856e24          # cm (หนังสือ หน้า 5)
Gc = 6.67e-8             # cm^3 g^-1 s^-2 (หนังสือ)
H0n = 100e5 / Mpc        # s^-1
check(f"ข้อ 6 ขั้น 1: H0 = {H0n:.4e} s^-1 ≈ 3.24e-18", abs(H0n / 3.24e-18 - 1) < 5e-3)
check(f"ข้อ 6 ขั้น 2: 8πG = {8*3.141592653589793*Gc:.4e} ≈ 1.68e-6", abs(8 * 3.141592653589793 * Gc / 1.68e-6 - 1) < 5e-3)
rc = 3 * H0n**2 / (8 * 3.141592653589793 * Gc)
check(f"ข้อ 6 ขั้น 3: ρcr = {rc:.4e} g/cm^3 ≈ 1.88e-29 (1.6)", abs(rc / 1.88e-29 - 1) < 5e-3)
# ถ้าใช้ค่าที่ปัดแล้วในขั้นก่อน ต้องยังอยู่ใน tolerance 0.5%
check("ใช้ค่าปัด 3.24 และ 1.68 ยังได้ 1.88 ภายใน 0.5%", abs(3 * 3.24**2 / 1.68 / 10 / 1.88 - 1) < 5e-3)
mp_g = 1.6726e-24        # g (นอกหนังสือ)
per_m3 = rc / mp_g * 1e6
check(f"proton ต่อลูกบาศก์เมตร = {per_m3:.2f} ≈ 11", abs(per_m3 - 11.2) < 0.2)

# ── ⑥ a_eq (ตัวเลขสมมุติ ρ_r/ρ_m วันนี้ = 1/1000) ─────────────
ratio = rho_r / rho_m
check("ข้อ 7 ขั้น 1: ρ_r/ρ_m ∝ a^-1", power_of(ratio, a) == -1)
aeq = sp.solve(sp.Eq(sp.Rational(1, 1000) / a, 1), a)
check("ข้อ 7 ขั้น 2: a_eq = 1/1000", aeq == [sp.Rational(1, 1000)])
check("ข้อ 7 ขั้น 3: z_eq = 1/a - 1 = 999", 1 / aeq[0] - 1 == 999)
check("ก่อน a_eq radiation มากกว่า (a = 1/2000 ⇒ อัตราส่วน = 2)",
      sp.Rational(1, 1000) / sp.Rational(1, 2000) == 2)

print(f"\nผ่าน {ok} ข้อ")
