# Python · matplotlib — รูปประกอบ Bayes (T1) · ตัวเลขทุกตัวคำนวณจริงและ assert
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family": ["Loma", "DejaVu Sans"], "font.size": 13, "axes.unicode_minus": False})
from fractions import Fraction as F

# ---------- รูป 1: กล่องสองใบ (แบบนับพื้นที่) ----------
pA, pB = F(1,2), F(1,2)          # prior: โยนเหรียญเลือกกล่อง
LA, LB = F(1,5), F(3,5)          # likelihood: โอกาสได้ลูกแดงจากแต่ละกล่อง
ev = LA*pA + LB*pB               # evidence
assert ev == F(2,5) and LB*pB/ev == F(3,4)

fig, axs = plt.subplots(1, 3, figsize=(15, 5.8))
def square(ax, title):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=14)
# (ก) prior
ax = axs[0]; square(ax, "① ก่อนหยิบ (prior)\nโยนเหรียญเลือกกล่อง")
ax.add_patch(Rectangle((0,0),.5,1,fc="#cfe0f5",ec="k")); ax.add_patch(Rectangle((.5,0),.5,1,fc="#f5dcc8",ec="k"))
ax.text(.25,.5,"กล่อง A\n½",ha="center",va="center",fontsize=16); ax.text(.75,.5,"กล่อง B\n½",ha="center",va="center",fontsize=16)
# (ข) likelihood
ax = axs[1]; square(ax, "② แต่ละกล่องมีลูกแดงกี่ส่วน\n(likelihood)")
ax.add_patch(Rectangle((0,0),.5,1,fc="#cfe0f5",ec="k")); ax.add_patch(Rectangle((.5,0),.5,1,fc="#f5dcc8",ec="k"))
ax.add_patch(Rectangle((0,0),.5,.2,fc="#d62728",ec="k")); ax.add_patch(Rectangle((.5,0),.5,.6,fc="#d62728",ec="k"))
ax.text(.25,.1,"0.2",ha="center",va="center",color="w",fontsize=15); ax.text(.75,.3,"0.6",ha="center",va="center",color="w",fontsize=15)
ax.text(.25,.6,"A",ha="center",fontsize=16); ax.text(.75,.8,"B",ha="center",fontsize=16)
# (ค) posterior
ax = axs[2]; square(ax, "③ ได้ลูกแดงแล้ว ⇒ เหลือแค่ส่วนสีแดง\n(posterior)")
ax.add_patch(Rectangle((0,0),.5,1,fc="#eeeeee",ec="#bbbbbb")); ax.add_patch(Rectangle((.5,0),.5,1,fc="#eeeeee",ec="#bbbbbb"))
ax.add_patch(Rectangle((0,0),.5,.2,fc="#d62728",ec="k")); ax.add_patch(Rectangle((.5,0),.5,.6,fc="#d62728",ec="k"))
ax.text(.25,.1,"0.1",ha="center",va="center",color="w",fontsize=15); ax.text(.75,.3,"0.3",ha="center",va="center",color="w",fontsize=15)
ax.text(.5,.82,"แดงทั้งหมด = 0.1 + 0.3 = 0.4  (evidence)\nP(B | แดง) = 0.3 / 0.4 = 0.75",ha="center",fontsize=13,
        bbox=dict(fc="w",ec="k"))
fig.tight_layout(); fig.savefig("fig1-boxes.png", dpi=110); plt.close(fig)

# ---------- รูป 2: prior × likelihood → posterior (Alice) ----------
mu, S2, m, s2 = 0.0, 4.0, 5.0, 1.0
tau2 = 1/(1/S2 + 1/s2); mubar = tau2*(mu/S2 + m/s2)
assert abs(tau2-0.8)<1e-12 and abs(mubar-4)<1e-12
th = np.linspace(-6, 10, 2001); dth = th[1]-th[0]
g = lambda x, c, v: np.exp(-(x-c)**2/(2*v))/np.sqrt(2*np.pi*v)
prior = g(th, mu, S2); like = np.exp(-(th-m)**2/(2*s2))   # likelihood ไม่ normalize ตั้งใจ
prod = prior*like; Z = prod.sum()*dth; post = prod/Z
assert np.max(np.abs(post - g(th, mubar, tau2))) < 1e-6      # ตรวจว่า posterior = N(4, 0.8) จริง
fig, axs = plt.subplots(1, 2, figsize=(15, 5.2))
ax = axs[0]
ax.plot(th, prior, lw=3, c="#1f77b4", label="prior  p(θ):  เชื่อว่าราว ๆ 0 ± 2")
ax.plot(th, like*prior.max(), lw=3, c="#d62728", ls="--", label="likelihood  L(θ):  วัดได้ 5 ± 1  (ย่อสเกลให้เห็น)")
ax.plot(th, prod/prod.max()*prior.max()*0.6, lw=3, c="#7f7f7f", label="prior × likelihood  (ยังไม่ normalize)")
ax.set_xlabel("θ"); ax.set_yticks([]); ax.legend(loc="upper left", fontsize=11); ax.set_title("① คูณกันทีละจุด θ")
ax = axs[1]
ax.plot(th, prior, lw=2, c="#1f77b4", alpha=.5, label="prior  N(0, 4)")
ax.plot(th, post, lw=3.5, c="#2ca02c", label="posterior  N(4, 0.8)")
ax.axvline(m, c="#d62728", ls="--", lw=1.5, label="ค่าที่วัดได้ m = 5")
ax.axvline(2.5, c="k", ls=":", lw=1.5, label="ค่าเฉลี่ยธรรมดา 2.5 (ไม่ใช่!)")
ax.annotate("", xy=(4,0.53), xytext=(0,0.53), arrowprops=dict(arrowstyle="->", lw=2))
ax.text(-5.5,0.56,"ถูกดึงจาก 0 ไปทางการวัด เพราะการวัดแม่นกว่า prior",ha="left",fontsize=11)
ax.set_xlabel("θ"); ax.set_yticks([]); ax.set_ylim(0,0.62); ax.set_xlim(-6,13); ax.legend(loc="center right", fontsize=11)
ax.set_title("② หารด้วยพื้นที่ใต้กราฟ (evidence) → posterior")
fig.tight_layout(); fig.savefig("fig2-alice.png", dpi=110); plt.close(fig)

th = np.linspace(-6, 14, 4001)
# ---------- รูป 3: Alice vs Bob, N = 1 กับ N = 99 ----------
def post_of(mu, S2, mbar, s2, N):
    t2 = 1/(1/S2 + N/s2); return t2*(mu/S2 + N*mbar/s2), t2
res = {}
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), sharey=False)
for ax, N in zip(axs, [0, 1, 99]):
    for name, (mu_, S2_), c in [("Alice", (0,4), "#1f77b4"), ("Bob", (10,1), "#ff7f0e")]:
        mb, t2 = (mu_, S2_) if N == 0 else post_of(mu_, S2_, 5, 1, N)
        res[(name,N)] = mb
        ax.plot(th, g(th, mb, t2), lw=3, c=c, label=f"{name}: กลาง {mb:.2f}")
        ax.fill_between(th, g(th, mb, t2), color=c, alpha=.15)
    ax.axvline(5, c="#d62728", ls="--", lw=1)
    ax.set_title("ก่อนวัด (prior)" if N==0 else f"หลังวัด N = {N} ครั้ง"); ax.set_xlabel("θ"); ax.set_yticks([]); ax.legend(fontsize=11)
assert abs(res[("Alice",1)]-4)<1e-12 and abs(res[("Bob",1)]-7.5)<1e-12
assert abs(res[("Bob",99)]-5.05)<1e-12 and abs(res[("Alice",99)]-495/99.25)<1e-12
fig.tight_layout(); fig.savefig("fig3-converge.png", dpi=110); plt.close(fig)
print("ok", float(ev), float(LB*pB/ev), tau2, mubar, res)
