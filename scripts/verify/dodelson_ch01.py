"""คู่มืออ่าน Dodelson บทที่ 1 — ตรวจขั้นพีชคณิต/ตัวเลขที่คู่มือเติมเอง (reference/dodelson-ch01.html)

รัน: python scripts/verify/dodelson_ch01.py   (ทุกบรรทัดต้องขึ้น ✅ · มี assert ล้มเมื่อผิด)
ค่าคงตัวที่ใช้เป็นค่าที่หนังสือให้ไว้ (ปี 2003) เว้นแต่จะเขียนกำกับว่า "นอกหนังสือ"
"""
import sys
import sympy as sp
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")    # Windows console ไม่ใช่ UTF-8 โดยปริยาย

ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


t, C, p = sp.symbols('t C p', positive=True)
a_ = sp.Function('a')

# ① นิยาม (1.1): H = (da/dt)/a · ทดลองกับกฎยกกำลัง a = C t^p
a = C * t**p
H = sp.simplify(sp.diff(a, t) / a)
check("a = C t^p  ⇒  H = p/t  (C หายไป)", sp.simplify(H - p / t) == 0)
check("matter: a ∝ t^(2/3) ⇒ H = 2/(3t)", sp.simplify(H.subs(p, sp.Rational(2, 3)) - sp.Rational(2, 3) / t) == 0)
check("radiation: a ∝ t^(1/2) ⇒ H = 1/(2t)", sp.simplify(H.subs(p, sp.Rational(1, 2)) - 1 / (2 * t)) == 0)
check("matter: H0 t0 = 2/3", sp.simplify((H * t).subs(p, sp.Rational(2, 3))) == sp.Rational(2, 3))

# ② (1.2) ที่ a = 1 (วันนี้) ให้ H = H0 พอดี เพราะพจน์ความโค้งเป็นศูนย์และ ρ0 = ρcr
G, rho0, H0 = sp.symbols('G rho0 H0', positive=True)
rho_cr = 3 * H0**2 / (8 * sp.pi * G)                       # (1.3)
H2_today = 8 * sp.pi * G / 3 * (rho0 + (rho_cr - rho0) / 1**2)
check("(1.2) ที่ a=1: H^2 = H0^2 ไม่ว่า ρ0 เป็นเท่าไร", sp.simplify(H2_today - H0**2) == 0)

# ③ ขั้นที่ข้าม: (1.2) แบน (ρ0 = ρcr) + matter ∝ a^-3 + Λ คงที่  ⇒  (1.10)
aa, OL = sp.symbols('a Omega_Lambda', positive=True)
rho = rho_cr * (OL + (1 - OL) / aa**3)
Hsq = 8 * sp.pi * G / 3 * rho
check("แบน+matter+Λ: H^2 = H0^2 [ΩΛ + (1-ΩΛ)/a^3]",
      sp.simplify(Hsq - H0**2 * (OL + (1 - OL) / aa**3)) == 0)
# dt = da/(aH)
dt_da = 1 / (aa * sp.sqrt(Hsq))
book = 1 / (H0 * aa) * (OL + (1 - OL) / aa**3) ** sp.Rational(-1, 2)
check("dt = da/(aH) ตรงกับรูป (1.10)", sp.simplify(dt_da - book) == 0)
# ΩΛ = 0 ต้องกลับไปได้ 2/3 (ตรงกับ ①) — คำตอบข้อ 2(a) ของหนังสือมีใน Appendix A อยู่แล้ว
age0 = sp.integrate((1 / (aa * sp.sqrt(1 / aa**3))), (aa, 0, 1))
check("ΩΛ=0: H0 t0 = ∫ a^(1/2) da = 2/3 (สอดคล้องกับ ①)", sp.simplify(age0) == sp.Rational(2, 3))

# ④ ตัวเลขใน (1.5) และ (1.6) — ค่าคงตัวจากหนังสือ: Mpc = 3.0856e24 cm, G = 6.67e-8 cgs
Mpc = 3.0856e24            # cm (หนังสือ หน้า 5)
H0_100 = 100e5 / Mpc       # s^-1 ที่ h = 1
yr = 3.156e7               # s  (นอกหนังสือ: 1 ปี ≈ 3.156e7 s)
hbar_eVs = 6.582e-16       # eV s (นอกหนังสือ)
check(f"1/H0 = {1/H0_100/yr/1e10:.3f}e10 h^-1 ปี ≈ 0.98e10  (1.5)", abs(1 / H0_100 / yr / 1e10 - 0.98) < 0.01)
check(f"ħH0 = {H0_100*hbar_eVs:.4e} h eV ≈ 2.133e-33  (1.5)", abs(H0_100 * hbar_eVs / 2.133e-33 - 1) < 2e-3)
Gcgs = 6.67e-8
rcr = 3 * H0_100**2 / (8 * mp.pi * Gcgs)
check(f"ρcr = {float(rcr):.3e} h^2 g/cm^3 ≈ 1.88e-29  (1.6)", abs(rcr / 1.88e-29 - 1) < 5e-3)
nprot = 1.88e-29 * 1e6 / 1.67e-24   # ต่อ m³ (มวล proton นอกหนังสือ)
check(f"ρcr/m_p = {nprot:.1f} h² proton ต่อ m³ ≈ 11", abs(nprot - 11.3) < 0.2)

# ⑤ อายุเอกภพแบน matter-dominated = (2/3)/H0 ที่ h = 0.72 ± 0.08 (ค่าปี 2003)
for h in (0.64, 0.72, 0.80):
    age = 2 / 3 * 0.98e10 / h / 1e9
    check(f"h={h}: (2/3)/H0 = {age:.1f} Gyr อยู่ในช่วง 8–10.5", 8 <= age <= 10.5)

# ⑥ Hubble (1929): ความชัน 1000 km/s ต่อ 2 Mpc · z = v/c
check("ความชัน = 1000/2 = 500 km/s/Mpc", 1000 / 2 == 500)
check(f"z = 1000/c = {1000/2.998e5:.4f} ≈ 0.003", abs(1000 / 2.998e5 - 0.0033) < 2e-4)

# ⑦ redshift (1.7): Lyα 1216 Å จาก QSO ที่ z = 3.572 ไปตกที่ ~5559 Å (รูป 1.9 แกน 5555–5560)
lam = 1216 * (1 + 3.572)
check(f"1216 (1+3.572) = {lam:.0f} Å อยู่ในช่วงแกนรูป 1.9", 5555 <= lam <= 5565)

# ⑧ อุณหภูมิวันนี้ 2.725 K เป็น eV (ใช้ 11605 K = 1 eV ตามหนังสือ หน้า 41)
check(f"2.725 K = {2.725/11605:.3e} eV ≈ 2.35e-4", abs(2.725 / 11605 / 2.35e-4 - 1) < 3e-3)

# ⑨ T = T0/a (1.4) ⇒ ρ_rad ∝ n·E ∝ a^-3 · a^-1 = a^-4
check("n ∝ a^-3, E ∝ a^-1 ⇒ ρ_rad ∝ a^-4", sp.simplify(aa**-3 * aa**-1 - aa**-4) == 0)

# ⑩ มิติของ P(k) ใน (1.9): ถ้า δ(x) ไม่มีมิติ  δ̃(k) = ∫d^3x ... มีมิติ L^3 · δ^3(k) มีมิติ L^3
L = sp.symbols('L', positive=True)
dim_delta_k = L**3
dim_lhs = dim_delta_k**2
dim_dirac = L**3                       # δ^3(k) = 1/[k]^3 = L^3
dim_P = sp.simplify(dim_lhs / dim_dirac)
check("[P(k)] = L^3  ⇒  k^3 P(k) ไม่มีมิติ", dim_P == L**3 and sp.simplify((1 / L)**3 * dim_P) == 1)

# ⑪ (1.8) blackbody: ยอดของ x^3/(e^x - 1) อยู่ที่ x ≈ 2.82 (ใช้ในโจทย์ข้อ 4 — Appendix A มีเฉลย)
xpk = mp.findroot(lambda x: 3 * (mp.e**x - 1) - x * mp.e**x, 2.8)
check(f"ยอดของ x^3/(e^x-1) ที่ x = {float(xpk):.3f}", abs(xpk - 2.821) < 1e-3)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
