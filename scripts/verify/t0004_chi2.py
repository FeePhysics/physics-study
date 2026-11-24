"""T4 · ตรวจทุกตัวเลขในบท t0004-chi2.html ด้วย sympy

รัน: python3 scripts/verify/t0004_chi2.py
"""
import sympy as sp

ok = 0
def check(name, got, want):
    global ok
    assert sp.simplify(got - want) == 0, f"{name}: ได้ {got} ต้องการ {want}"
    ok += 1
    print("✓", name, "=", want)

d, t, s, th, a, b = sp.symbols("d t sigma theta a b", real=True)
sig = sp.symbols("sigma", positive=True)

# ── (1) Gaussian likelihood ⇒ −2 ln L = ((d − t)/σ)² + ค่าคงที่ (ไม่ขึ้นกับ t) ──
L = sp.exp(-(d - t) ** 2 / (2 * sig ** 2)) / sp.sqrt(2 * sp.pi * sig ** 2)
rest = sp.expand(sp.expand_log(-2 * sp.log(L), force=True) - ((d - t) / sig) ** 2)
assert not rest.has(t) and not rest.has(d)
ok += 1; print("✓ −2 ln L − pull² ไม่ขึ้นกับทฤษฎี:", sp.simplify(rest))

# ── (3) ผลคูณของ likelihood อิสระ ⇒ χ² บวกกัน ──
L1, L2 = sp.symbols("L1 L2", positive=True)
check("−2 ln(L1 L2) = −2 ln L1 − 2 ln L2", sp.expand_log(-2 * sp.log(L1 * L2)), -2 * sp.log(L1) - 2 * sp.log(L2))

# ── ของเล่น (ภาพ warmup ใช้ข้อมูลเดียวกัน) ──
xs = [10, 20, 30, 40]; ys = [700, 1500, 2100, 3600]; ss = [100, 100, 100, 600]
k = sp.symbols("k")
chi2 = sum(((sp.Integer(y) - k * x) / sp.Integer(s_)) ** 2 for x, y, s_ in zip(xs, ys, ss))
kbest = sp.solve(sp.diff(chi2, k), k)[0]
check("ของเล่น H0 ที่ดีที่สุด", kbest, 72)
check("ของเล่น χ² ต่ำสุด", chi2.subs(k, kbest), sp.Rational(11, 5))
check("warmup เส้น A (H0 = 70) χ²", chi2.subs(k, 70), 1 + sp.Rational(16, 9))
check("warmup เส้น B (H0 = 90) χ²", chi2.subs(k, 90), 4 + 9 + 36 + 0)
# พาราโบลา: χ²(k) = χ²min + ((k − k̂)/σ_k)² · σ_k = 1/√(Σx²/σ²)
sk = 1 / sp.sqrt(sum(sp.Integer(x) ** 2 / sp.Integer(s_) ** 2 for x, s_ in zip(xs, ss)))
check("χ² เป็นพาราโบลาเป๊ะ", sp.expand(chi2 - (sp.Rational(11, 5) + ((k - 72) / sk) ** 2)), 0)
print("   σ ของ H0 ในของเล่น ≈", float(sk))
# (4) ห่าง nσ ⇒ Δχ² = n²
for n in (1, 2):
    check(f"Δχ² ที่ห่าง {n}σ", chi2.subs(k, 72 + n * sk) - chi2.subs(k, 72), n ** 2)

# ── ข้อ 1 · จุดเดียว 14 ± 2 · ทฤษฎี 10 ──
p1 = sp.Rational(14 - 10, 2)
check("ข้อ1 ขั้น1 ห่างกี่ σ", p1, 2)
check("ข้อ1 ขั้น2 χ²", p1 ** 2, 4)
check("ข้อ1 ขั้น3 L ทฤษฎี 10 ÷ L ทฤษฎี 14", sp.exp(-p1 ** 2 / 2) / sp.exp(0), sp.exp(-2))
print("   e^-2 ≈", float(sp.exp(-2)))

# ── ข้อ 2 · d = 10, 20, 30 · v = 750, 1300, 2250 · σ = 50, 100, 150 ──
X2 = [10, 20, 30]; Y2 = [750, 1300, 2250]; S2 = [50, 100, 150]
pulls = lambda H: [sp.Rational(y - H * x, s_) for x, y, s_ in zip(X2, Y2, S2)]
check("ข้อ2 ขั้น1 จุดแรกห่างกี่ σ (H0=70)", pulls(70)[0], 1)
check("ข้อ2 ขั้น2 จุดที่สอง (H0=70)", pulls(70)[1], -1)
check("ข้อ2 ขั้น3 χ² ที่ H0=70", sum(p ** 2 for p in pulls(70)), 3)
check("ข้อ2 ขั้น4 χ² ที่ H0=75", sum(p ** 2 for p in pulls(75)), 4)

# ── ข้อ 3 · (45) รวมสองชุดข้อมูล ──
H = [66, 70, 74]; SN = [2, 5, 10]; LE = [9, 4, 3]
tot = [a_ + b_ for a_, b_ in zip(SN, LE)]
check("ข้อ3 SN ชอบ", H[SN.index(min(SN))], 66)
check("ข้อ3 lensing ชอบ", H[LE.index(min(LE))], 74)
check("ข้อ3 χ²tot ที่ 70", tot[1], 9)
check("ข้อ3 รวมกันชอบ", H[tot.index(min(tot))], 70)

# ── ข้อ 4 · ตาราง χ²(H0) = 3 + ((H0 − 70)/2)² ──
tab = {h: 3 + sp.Rational(h - 70, 2) ** 2 for h in (66, 68, 70, 72, 74)}
assert [tab[h] for h in (66, 68, 70, 72, 74)] == [7, 4, 3, 4, 7]
ok += 1; print("✓ ข้อ4 ตาราง 7 4 3 4 7")
check("ข้อ4 ช่วง 1σ: Δχ²=1 ที่ ±2", sp.solve(sp.Eq(sp.Rational(1, 4) * (th - 70) ** 2, 1), th)[1] - 70, 2)
check("ข้อ4 ขอบล่าง 2σ", sp.solve(sp.Eq(sp.Rational(1, 4) * (th - 70) ** 2, 4), th)[0], 66)

print(f"\nผ่าน {ok} ข้อ")
