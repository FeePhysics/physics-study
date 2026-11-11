"""ตรวจทุกข้ออ้างของบทที่ 24 — Bayes กับพารามิเตอร์ต่อเนื่อง (Trotta 2008 §2.2-2.3)

prior N(mu, Sigma^2) × likelihood N(m | theta, sigma^2) ⇒ posterior N(mubar, tau^2)
  1/tau^2 = 1/Sigma^2 + 1/sigma^2           (precision บวกกัน)
  mubar   = (mu/Sigma^2 + m/sigma^2) * tau^2 (ค่าเฉลี่ยถ่วงด้วย precision)
และต้องเท่ากับ Trotta (9)-(10) ที่เขียนด้วย (sigma/Sigma)^2
"""
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(got - want) == 0, got, want))


th, mu, m, mbar = sp.symbols('theta mu m mbar', real=True)
Sig, sig, N = sp.symbols('Sigma sigma N', positive=True)

prior = sp.exp(-(th - mu)**2 / (2*Sig**2)) / sp.sqrt(2*sp.pi*Sig**2)
like = sp.exp(-(th - m)**2 / (2*sig**2))          # Trotta (7) · L0 ไม่สำคัญ

prec = 1/Sig**2 + 1/sig**2
tau2 = 1/prec
mub = (mu/Sig**2 + m/sig**2) * tau2

# ── สมการ (2)(3): เติมกำลังสองจริง ────────────────────────────
# log ของผลคูณ = ผลบวกของเลขชี้กำลัง (ตัดค่าคงตัวที่ไม่มี θ ทิ้งได้)
neg2log = sp.expand(-2*(-(th - mu)**2/(2*Sig**2) - (th - m)**2/(2*sig**2)))
ok.append(("prior×likelihood = exp(เลขชี้กำลังข้างบน) × ค่าคงตัว",
           sp.simplify(prior*like/sp.exp(-neg2log/2)).has(th) is False, None, None))
c2 = neg2log.coeff(th, 2)
c1 = neg2log.coeff(th, 1)
chk("สัมประสิทธิ์ θ² ของ -2ln(prior×likelihood) = 1/Σ² + 1/σ²", c2, prec)
chk("สัมประสิทธิ์ θ¹ = -2(μ/Σ² + m/σ²)", c1, -2*(mu/Sig**2 + m/sig**2))
# Gaussian (θ-μ̄)²/τ² มี θ² = 1/τ² และ θ¹ = -2μ̄/τ²  ⇒ อ่าน τ², μ̄ ออกมา
chk("τ² = 1/(สัมประสิทธิ์ θ²)", 1/c2, tau2)
chk("μ̄ = -(สัมประสิทธิ์ θ¹)/(2·สัมประสิทธิ์ θ²)", -c1/(2*c2), mub)

# posterior ที่ normalize แล้วเป็น Gaussian N(μ̄, τ²) จริง (อินทิเกรตตรง)
post_un = prior*like
Z = sp.simplify(sp.integrate(post_un, (th, -sp.oo, sp.oo)))
post = post_un / Z
gauss = sp.exp(-(th - mub)**2/(2*tau2)) / sp.sqrt(2*sp.pi*tau2)
vals = {mu: 0, Sig: 2, m: 5, sig: 1}
for tv in [-1, 0, 3, sp.Rational(7, 2), 6]:
    chk(f"posterior = N(μ̄,τ²) ที่ θ={tv} (Alice)", sp.N(post.subs(vals).subs(th, tv), 30),
        sp.N(gauss.subs(vals).subs(th, tv), 30))

# ── ตรงกับ Trotta (9)(10) ────────────────────────────────────
r = (sig/Sig)**2
chk("Trotta (9): (m + (σ/Σ)²μ)/(1+(σ/Σ)²) = μ̄", (m + r*mu)/(1 + r), mub)
chk("Trotta (10): σ²/(1+(σ/Σ)²) = τ²", sig**2/(1 + r), tau2)

# ── ข้อ 1: สองสมมติฐาน (Trotta (5)(6)) ──────────────────────
pA, pB = sp.Rational(1, 2), sp.Rational(1, 2)
LA, LB = sp.Rational(1, 5), sp.Rational(3, 5)
ev = LA*pA + LB*pB
chk("ข้อ 1 ขั้น 1: evidence p(d) = 0.4", ev, sp.Rational(2, 5))
chk("ข้อ 1 ขั้น 2: p(B|d) = 0.75", LB*pB/ev, sp.Rational(3, 4))
ok.append(("ตัวควบคุม: ลืมหาร evidence ได้ 0.3 ≠ 0.75", LB*pB != LB*pB/ev, LB*pB, None))

# ── ข้อ 2: Alice μ=0 Σ²=4 · ข้อมูล m=5 σ²=1 ─────────────────
A = {mu: 0, Sig: 2, m: 5, sig: 1}
chk("ข้อ 2 ขั้น 1: 1/Σ² = 1/4", (1/Sig**2).subs(A), sp.Rational(1, 4))
chk("ข้อ 2 ขั้น 2: 1/τ² = 5/4", prec.subs(A), sp.Rational(5, 4))
chk("ข้อ 2 ขั้น 3: τ² = 4/5", tau2.subs(A), sp.Rational(4, 5))
chk("ข้อ 2 ขั้น 4: μ̄ = 4", mub.subs(A), 4)
chk("ข้อ 2 ขั้น 5: (σ/Σ)² = 1/4", r.subs(A), sp.Rational(1, 4))
chk("ข้อ 2 ขั้น 6: Trotta (9) ได้ 4", ((m + r*mu)/(1 + r)).subs(A), 4)
ok.append(("ตัวควบคุม: บวกความแปรปรวนแทน precision ได้ τ²=5 ≠ 4/5",
           (Sig**2 + sig**2).subs(A) != tau2.subs(A), (Sig**2 + sig**2).subs(A), None))
ok.append(("ตัวควบคุม: เฉลี่ยธรรมดา (0+5)/2 = 5/2 ≠ 4",
           sp.Rational(5, 2) != mub.subs(A), sp.Rational(5, 2), None))

# ── ข้อ 3: Bob μ=10 Σ²=1 · ข้อมูลเดิม ─────────────────────────
B = {mu: 10, Sig: 1, m: 5, sig: 1}
chk("ข้อ 3 ขั้น 1: 1/τ² = 2", prec.subs(B), 2)
chk("ข้อ 3 ขั้น 2: τ² = 1/2", tau2.subs(B), sp.Rational(1, 2))
chk("ข้อ 3 ขั้น 3: μ̄ = 15/2", mub.subs(B), sp.Rational(15, 2))
chk("ระยะห่าง Alice–Bob หลังข้อมูลตัวเดียว = 7/2", mub.subs(B) - mub.subs(A), sp.Rational(7, 2))

# ── ข้อมูล N ตัว: likelihood ของ θ เป็น Gaussian ความแปรปรวน σ²/N รอบค่าเฉลี่ย ──
xs = sp.symbols('x1:5', real=True)        # ทดสอบ N = 4 ด้วยตัวแปรจริง
n = len(xs)
s = sum(-(th - xi)**2/(2*sig**2) for xi in xs)
xbar = sum(xs)/n
rest = sp.expand(s - (-(th - xbar)**2/(2*sig**2/n)))
chk("Σ-(θ-xᵢ)²/2σ² = -(θ-x̄)²/(2σ²/N) + (ไม่มี θ)  [N=4]", sp.diff(rest, th), 0)

# ── ข้อ 4: Bob กับข้อมูล N = 99 ตัว ค่าเฉลี่ย 5 ─────────────────
B99 = {mu: 10, Sig: 1, m: 5, sig: 1/sp.sqrt(99)}
chk("ข้อ 4 ขั้น 1: precision ของข้อมูล N/σ² = 99", (1/sig**2).subs(B99), 99)
chk("ข้อ 4 ขั้น 2: 1/τ² = 100", prec.subs(B99), 100)
chk("ข้อ 4 ขั้น 3: μ̄ = 101/20 = 5.05", mub.subs(B99), sp.Rational(101, 20))
A99 = {mu: 0, Sig: 2, m: 5, sig: 1/sp.sqrt(99)}
gap = mub.subs(B99) - mub.subs(A99)
ok.append(("Alice–Bob ห่างกัน < 0.07 ที่ N=99 (จาก 3.5 ที่ N=1)", gap < sp.Rational(7, 100), sp.N(gap), None))
ok.append(("ตัวควบคุม: ลืมคูณ N ได้ μ̄ = 15/2 ≠ 5.05", mub.subs(B) != mub.subs(B99), mub.subs(B), None))

# ── ตัวอย่างในเนื้อบท: μ=0 Σ²=9 · m=10 σ²=1 ─────────────────
E = {mu: 0, Sig: 3, m: 10, sig: 1}
chk("ตัวอย่าง: 1/τ² = 10/9", prec.subs(E), sp.Rational(10, 9))
chk("ตัวอย่าง: τ² = 0.9", tau2.subs(E), sp.Rational(9, 10))
chk("ตัวอย่าง: μ̄ = 9", mub.subs(E), 9)
ok.append(("ข้อ 4: Alice ที่ N=99 ปัดได้ 4.99", round(float(mub.subs(A99)), 2) == 4.99, sp.N(mub.subs(A99)), None))

bad = 0
for name, good, got, want in ok:
    print(("✅" if good else "❌"), name, "" if good else f"got={got} want={want}")
    bad += not good
print(f"\n{len(ok) - bad}/{len(ok)}")
raise SystemExit(1 if bad else 0)
