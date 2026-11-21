"""Python · matplotlib — รูปประกอบ T0 (Laplace) และ T1 (Gaussian) · ฝังเป็น inline SVG ในบท
ตัวเลขในรูปตั้งใจให้ "ไม่ซ้ำ" กับข้อฝึก ไม่งั้นรูปจะเฉลยคำตอบ · ทุกค่า assert
รัน: python3 scripts/verify/t0001_figs.py  → เขียน SVG ทับระหว่าง <!--FIG:x--> ... <!--/FIG:x--> ในบท
"""
import io, re, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from fractions import Fraction as F
plt.rcParams.update({"font.family": ["Loma", "DejaVu Sans"], "font.size": 17,
                     "axes.unicode_minus": False, "svg.fonttype": "none", "svg.hashsalt": "t0001"})
g = lambda x, c, v: np.exp(-(x-c)**2/(2*v))/np.sqrt(2*np.pi*v)
figs = {}

def to_svg(fig, label):
    buf = io.StringIO(); fig.savefig(buf, format="svg", facecolor="white"); plt.close(fig)
    s = buf.getvalue()
    s = s[s.index("<svg"):]                                  # ตัด XML header / DOCTYPE
    s = re.sub(r"<metadata>.*?</metadata>", "", s, flags=re.S)
    s = re.sub(r'(<svg[^>]*?) width="[^"]*" height="[^"]*"', r'\1', s, count=1)  # ให้ CSS คุมขนาด
    s = s.replace("<svg ", f'<svg role="img" aria-label="{label}" ', 1)
    assert "$" not in s                                       # กฎ workspace: ห้าม $ ใน svg
    return s

# รูปกล่องสองใบเดิม → แทนด้วยตาราง 100 รอบแบบเลื่อนได้ (assets/bayes.js) ใน T0

# ---------- รูป 2: ตัวอย่างในเนื้อบท μ=0 Σ²=9 · m=10 σ²=1 → μ̄=9 τ²=0.9 ----------
mu, S2, m, s2 = 0.0, 9.0, 10.0, 1.0
tau2 = 1/(1/S2 + 1/s2); mubar = tau2*(mu/S2 + m/s2)
assert abs(tau2 - 0.9) < 1e-12 and abs(mubar - 9) < 1e-12
th = np.linspace(-10, 16, 801); dth = th[1] - th[0]
prior = g(th, mu, S2); like = np.exp(-(th - m)**2/(2*s2))
prod = prior*like; post = prod/(prod.sum()*dth)
assert np.max(np.abs(post - g(th, mubar, tau2))) < 1e-6     # normalize จริงแล้วตรงกับ N(9, 0.9)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
ax = axs[0]
ax.plot(th, prior/prior.max(), lw=3, c="#1f77b4", label="prior: เชื่อว่าราว 0 ± 3")
ax.plot(th, like, lw=3, c="#d62728", ls="--", label="likelihood: วัดได้ 10 ± 1")
ax.plot(th, prod/prod.max()*0.6, lw=3, c="#555555", label="prior × likelihood (ย่อสเกล)")
ax.set_title("① คูณกันทีละค่า θ"); ax.set_xlabel("θ"); ax.set_yticks([]); ax.set_ylim(0, 1.45)
ax.legend(loc="upper left", fontsize=14, frameon=False)
ax = axs[1]
ax.plot(th, prior, lw=2, c="#1f77b4", alpha=.6, label="prior")
ax.plot(th, post, lw=3.5, c="#2ca02c", label="posterior: กลาง 9")
ax.axvline(m, c="#d62728", ls="--", lw=1.5, label="ค่าที่วัดได้ 10")
ax.axvline(5, c="k", ls=":", lw=1.5, label="ค่าเฉลี่ยธรรมดา 5 (ไม่ใช่)")
ax.set_title("② หารด้วยพื้นที่ใต้กราฟ (evidence)"); ax.set_xlabel("θ"); ax.set_yticks([]); ax.set_ylim(0, .62)
ax.legend(loc="upper left", fontsize=14, frameon=False)
fig.tight_layout()
figs["product"] = to_svg(fig, "prior กว้างรอบ 0 คูณ likelihood แคบรอบ 10 ได้ posterior แคบรอบ 9 ไม่ใช่ 5")

# ---------- รูป 3: Alice (0, 9) กับ Bob (20, 1) · ข้อมูลเฉลี่ย 10 · N = 0, 1, 99 ----------
def post_of(mu, S2, mbar, s2, N):
    t2 = 1/(1/S2 + N/s2); return t2*(mu/S2 + N*mbar/s2), t2
res = {}
th = None
fig, axs = plt.subplots(1, 3, figsize=(12, 4.4))
for ax, N in zip(axs, [0, 1, 99]):
    for name, (mu_, S2_), c in [("Alice", (0, 9), "#1f77b4"), ("Bob", (20, 1), "#ff7f0e")]:
        mb, t2 = (mu_, S2_) if N == 0 else post_of(mu_, S2_, 10, 1, N)
        res[(name, N)] = mb
        # จุดถี่เฉพาะรอบยอด ไม่งั้น SVG บวมเป็น MB
        x = np.unique(np.concatenate([np.linspace(-6, 24, 200), np.linspace(mb - 6*t2**.5, mb + 6*t2**.5, 200)]))
        y = g(x, mb, t2); ax.plot(x, y/y.max(), lw=3, c=c)
        ax.fill_between(x, y/y.max(), color=c, alpha=.15)
    ax.axvline(10, c="#d62728", ls="--", lw=1)
    head = "ก่อนวัด" if N == 0 else f"หลังวัด {N} ครั้ง"
    ax.set_title(f"{head}\nAlice {res[('Alice', N)]:.2f} · Bob {res[('Bob', N)]:.2f}")
    ax.set_xlabel("θ"); ax.set_yticks([]); ax.set_ylim(0, 1.1); ax.set_xlim(-10, 25)
assert abs(res[("Alice", 1)] - 9) < 1e-12 and abs(res[("Bob", 1)] - 15) < 1e-12
assert abs(res[("Bob", 99)] - 10.1) < 1e-12 and abs(res[("Alice", 99)] - 990/(1/9 + 99)) < 1e-12
fig.tight_layout()
figs["converge"] = to_svg(fig, "Alice เริ่มที่ 0 Bob เริ่มที่ 20 หลังวัด 1 ครั้งได้ 9 กับ 15 หลังวัด 99 ครั้งได้ราว 9.99 กับ 10.1")


# ---------- รูป 0: เดิมพันของ Laplace (ดาวเสาร์ชนะ · ดาวพฤหัสแพ้) ----------
# Laplace: เสาร์ 1/3512 · พฤหัส 1/1070 ของมวลดวงอาทิตย์ · พนันว่าคลาดไม่ถึง 1/100
# ค่าปัจจุบัน (IAU refined): ดวงอาทิตย์/เสาร์ = 3497.898 · ดวงอาทิตย์/พฤหัส = 1047.3486
bets = [("ดาวเสาร์\n11,000 : 1", 3512, 3497.898), ("ดาวพฤหัส\n1,000,000 : 1", 1070, 1047.3486)]
dev = {}
fig, ax = plt.subplots(figsize=(7.5, 3.4))
for i, (name, L, M) in enumerate(bets):
    d = (1/M - 1/L)/(1/L)*100                 # ค่าจริงห่างจากค่าของ Laplace กี่ %
    dev[L] = d; y = 1 - i
    ax.add_patch(Rectangle((-1, y - .18), 2, .36, fc="#c9e7c9", ec="#2ca02c"))
    ax.plot([0], [y], "o", ms=10, c="k")
    ax.plot([d], [y], "*", ms=22, c="#d62728" if abs(d) > 1 else "#2ca02c", mec="k")
    ax.text(-3.5, y, name, ha="left", va="center", fontsize=17)
    ax.text(d, y + .3, ("แพ้" if abs(d) > 1 else "ชนะ") + f"  (ค่าจริง {d:+.1f}%)", ha="center", fontsize=15)
assert abs(dev[3512] - 0.40) < 0.01 and abs(dev[1070] - 2.16) < 0.01
ax.set_xlim(-3.6, 3.4); ax.set_ylim(-.35, 1.6); ax.set_yticks([])
ax.set_xlabel("ห่างจากค่าของ Laplace (%)")
for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
fig.tight_layout()
figs["laplace"] = to_svg(fig, "Laplace พนันว่ามวลดาวเสาร์คลาดไม่ถึง 1 เปอร์เซ็นต์ ค่าจริงห่าง 0.4 เปอร์เซ็นต์ ชนะ ดาวพฤหัสค่าจริงห่าง 2.2 เปอร์เซ็นต์ แพ้")

# ---------- ฝังลงบท ----------
WHERE = {"laplace": "lessons/t0000-bayes-boxes.html",
         "product": "lessons/t0001-bayes-gauss.html", "converge": "lessons/t0001-bayes-gauss.html"}
for k, s_ in figs.items():
    path = WHERE[k]; html = open(path, encoding="utf-8").read()
    pat = re.compile(rf"(<!--FIG:{k}-->).*?(<!--/FIG:{k}-->)", re.S)
    assert pat.search(html), f"ไม่เจอ marker FIG:{k} ใน {path}"
    html = pat.sub(lambda mt: mt.group(1) + "\n" + s_ + "\n" + mt.group(2), html)
    open(path, "w", encoding="utf-8").write(html)
print("ok", {k: len(v)//1024 for k, v in figs.items()}, "KB", res)
