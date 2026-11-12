"""Python · matplotlib — รูปประกอบ T1 (ฝังเป็น inline SVG ในบท)
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
LESSON = "lessons/t0001-bayes-gauss.html"
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

# ---------- รูป 1: กล่องสองใบ (prior 1/2 · แดง 0.3 กับ 0.7) ----------
pA = pB = F(1, 2); LA, LB = F(3, 10), F(7, 10)
ev = LA*pA + LB*pB
assert ev == F(1, 2) and LA*pA == F(3, 20) and LB*pB == F(7, 20) and LB*pB/ev == F(7, 10)
fig, axs = plt.subplots(1, 3, figsize=(11, 5.0))
def box(ax, title, grey=False):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=17)
    ax.add_patch(Rectangle((0, 0), .5, 1, fc="#eeeeee" if grey else "#cfe0f5", ec="k"))
    ax.add_patch(Rectangle((.5, 0), .5, 1, fc="#eeeeee" if grey else "#f5dcc8", ec="k"))
box(axs[0], "① ก่อนหยิบ = prior\nโยนเหรียญเลือกกล่อง")
axs[0].text(.25, .5, "กล่อง A\n½", ha="center", va="center", fontsize=19)
axs[0].text(.75, .5, "กล่อง B\n½", ha="center", va="center", fontsize=19)
for i, (grey, a, b, t) in enumerate([(False, "0.3", "0.7", "② ถ้าหยิบจากกล่องนี้ ได้แดงกี่ส่วน\n= likelihood"),
                                     (True, "0.15", "0.35", "③ เห็นว่าได้แดง ตัดส่วนอื่นทิ้ง\n= posterior")], start=1):
    box(axs[i], t, grey)
    axs[i].add_patch(Rectangle((0, 0), .5, float(LA), fc="#d62728", ec="k"))
    axs[i].add_patch(Rectangle((.5, 0), .5, float(LB), fc="#d62728", ec="k"))
    axs[i].text(.25, float(LA)/2, a, ha="center", va="center", color="w", fontsize=18)
    axs[i].text(.75, float(LB)/2, b, ha="center", va="center", color="w", fontsize=18)
axs[1].text(.25, .65, "A", ha="center", fontsize=19); axs[1].text(.75, .85, "B", ha="center", fontsize=19)
axs[2].text(.5, .80, "แดงทั้งหมด 0.15 + 0.35 = 0.5\n(evidence)\nB ในส่วนแดง 0.35 / 0.5 = 0.7",
            ha="center", va="center", fontsize=15, bbox=dict(fc="w", ec="k"))
fig.tight_layout()
figs["boxes"] = to_svg(fig, "สี่เหลี่ยมพื้นที่ 1 แบ่งครึ่งเป็นกล่อง A กับ B แต่ละครึ่งมีส่วนแดง 0.3 และ 0.7 พอเห็นว่าได้แดงก็เหลือแค่ส่วนแดง 0.15 กับ 0.35 ได้ B เท่ากับ 0.7")

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

# ---------- ฝังลงบท ----------
html = open(LESSON, encoding="utf-8").read()
for k, s in figs.items():
    pat = re.compile(rf"(<!--FIG:{k}-->).*?(<!--/FIG:{k}-->)", re.S)
    assert pat.search(html), f"ไม่เจอ marker FIG:{k} ในบท"
    html = pat.sub(lambda mt: mt.group(1) + "\n" + s + "\n" + mt.group(2), html)
open(LESSON, "w", encoding="utf-8").write(html)
print("ok", {k: len(v)//1024 for k, v in figs.items()}, "KB", res)
