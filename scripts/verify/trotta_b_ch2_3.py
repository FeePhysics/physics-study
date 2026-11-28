"""คู่มืออ่าน [B] Trotta, Bayesian Methods in Cosmology §2–§3.4 · ตรวจขั้นที่คู่มือเติมให้ด้วย sympy

รัน: python3 scripts/verify/trotta_b_ch2_3.py
"""
import sympy as sp

ok = 0
def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✓", name)

N = 4  # ใช้ N จริงเล็ก ๆ ให้ sympy กระจายได้ · ผลเป็นรูปทั่วไปตรวจซ้ำด้วยสัญลักษณ์ข้างล่าง
mu = sp.symbols("mu", real=True)
sig = sp.symbols("sigma", positive=True)
x = sp.symbols("x1:%d" % (N + 1), real=True)

# (19)-(21): −2 ln L = χ² + ค่าคงที่ไม่ขึ้นกับ μ
lnL = sum(-sp.log(sp.sqrt(2 * sp.pi) * sig) - (xi - mu) ** 2 / (2 * sig**2) for xi in x)
chi2 = sum((xi - mu) ** 2 / sig**2 for xi in x)
check("(20)-(21) −2 ln L − χ² ไม่ขึ้นกับ μ", sp.diff(sp.expand(-2 * lnL - chi2), mu) == 0)

# (22): ∂lnL/∂μ = 0 ⇒ μ_ML = ค่าเฉลี่ย
sol = sp.solve(sp.diff(lnL, mu), mu)[0]
check("(22) μ_ML = ค่าเฉลี่ยของข้อมูล", sp.simplify(sol - sum(x) / N) == 0)

# (33)-(34): 1/Σ² = −∂²lnL/∂μ² = N/σ²
check("(33)→(34) Σ² = σ²/N", sp.simplify(-1 / sp.diff(lnL, mu, 2) - sig**2 / N) == 0)

# (31)-(32): ln L เป็นพาราโบลาใน μ เป๊ะ ⇒ การประมาณด้วยอนุกรมเทย์เลอร์อันดับสองไม่มีพจน์ที่ตกหล่น
check("(31) อนุพันธ์อันดับสามของ ln L เป็นศูนย์ (Gaussian)", sp.diff(lnL, mu, 3) == 0)

# (23): σ²_ML = (1/N)Σ(x−μ)²
s = sp.symbols("s", positive=True)
lnLs = sum(-sp.log(s) - (xi - mu) ** 2 / (2 * s**2) for xi in x)
sol_s = [r for r in sp.solve(sp.diff(lnLs, s), s) if r.is_positive is not False][0]
check("(23) σ²_ML = (1/N)Σ(x−μ)²", sp.simplify(sol_s**2 - sum((xi - mu) ** 2 for xi in x) / N) == 0)

# (27): binomial ⇒ θ_ML = H/N
th, H, n = sp.symbols("theta H N", positive=True)
lnB = H * sp.log(th) + (n - H) * sp.log(1 - th)
check("(27) θ_ML = H/N", sp.simplify(sp.solve(sp.diff(lnB, th), th)[0] - H / n) == 0)

# (30): Poisson ⇒ λ_ML = n/t
lam, t, k = sp.symbols("lambda t n", positive=True)
lnP = k * sp.log(lam * t) - lam * t
check("(30) λ_ML = n/t", sp.simplify(sp.solve(sp.diff(lnP, lam), lam)[0] - k / t) == 0)

# (68)-(70): ความคลาดไม่เท่ากัน ⇒ ค่าเฉลี่ยถ่วง 1/σ_i² และ Σ⁻² = Σ 1/σ_i²
F = sp.symbols("F", real=True)
si = sp.symbols("s1:%d" % (N + 1), positive=True)
lnLw = sum(-(F - xi) ** 2 / (2 * sj**2) for xi, sj in zip(x, si))
Fml = sp.solve(sp.diff(lnLw, F), F)[0]
w = [1 / sj**2 for sj in si]
check("(68) F_ML = Σ(x_i/σ_i²)/Σ(1/σ_i²)", sp.simplify(Fml - sum(xi * wi for xi, wi in zip(x, w)) / sum(w)) == 0)
check("(70) Σ⁻² = Σ 1/σ_i²", sp.simplify(-sp.diff(lnLw, F, 2) - sum(w)) == 0)

# (72)+(74): prior แบน ⇒ posterior ∝ likelihood · ตรวจกับตัวอย่าง: ระฆัง × ค่าคงที่ แล้วทำให้รวมได้ 1
z = sp.symbols("z", real=True)
Lz = sp.exp(-(z - 3) ** 2 / 2)
post = Lz * sp.Rational(1, 10) / sp.integrate(Lz * sp.Rational(1, 10), (z, -sp.oo, sp.oo))
check("(74) prior แบน: posterior = likelihood ที่ทำให้รวมได้ 1",
      sp.simplify(post - Lz / sp.integrate(Lz, (z, -sp.oo, sp.oo))) == 0)

# ตัวอย่าง 14 ของ [B]: ข้อมูล (39) ⇒ 198.0 < μ < 201.2 ที่ 68.3%
d39 = [sp.Rational(v) for v in "197.2 202.4 201.8 198.8 207.6 191.4 201.4 198.2 195.7 201.2".split()]
m39 = sum(d39) / 10
S39 = 5 / sp.sqrt(10)
check("ตัวอย่าง 14: ช่วง 1σ ปัดแล้ว 198.0 ถึง 201.2",
      round(float(m39 - S39), 1) == 198.0 and round(float(m39 + S39), 1) == 201.2)
print("   ค่าเฉลี่ยของ (39) =", float(m39))

print(f"\nผ่าน {ok} ข้อ")
