"""รูปประกอบ T3 (SVG เขียนมือ · ไม่มี $ ข้างใน · ข้อความใช้ currentColor ให้ dark mode อ่านได้)

รัน: python3 scripts/verify/t0003_figs.py  → เขียน SVG ทับระหว่าง <!--FIG:x--> ... <!--/FIG:x--> ในบท

  ex3       ระฆังสองลูกของข้อ 3 (แสงรวม 12 ± 3 · พื้นหลัง 5 ± 4) — ไม่บอกคำตอบ
  peaks     ซ้าย: posterior ใกล้ Gaussian ยอดภูเขากับยอดเงาอยู่ที่ s เดียวกัน
            ขวา: ตารางใน ⓪ หมุนให้ s เป็นแนวนอน ช่องสูงสุด (s=1) ≠ ผลรวมคอลัมน์สูงสุด (s=2)
  triangle  triangle plot ของตัวอย่างกาแล็กซี 10 ± 1 · พื้นหลัง 3 ± 2 (แบบเดียวกับ Cobaya รูป 5)
"""
import math
import re
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
LESSON = ROOT / "lessons" / "t0003-marginal.html"
BLUE, RED, GREEN = "#3b6fb6", "#d62728", "#2ca02c"
FONT = 'font-family="Sarabun, sans-serif" fill="currentColor"'


def gauss(x, m, s):
    return math.exp(-((x - m) ** 2) / (2 * s * s))


def path(pts):
    return "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)


def svg(w, h, label, body):
    return (f'<figure><svg role="img" aria-label="{label}" viewBox="0 0 {w} {h}" '
            f'xmlns="http://www.w3.org/2000/svg" style="max-width:{w}px">{body}</svg></figure>')


# ── ex3: ระฆังสองลูกบนแกนเดียวกัน ────────────────────────────────────────────
def fig_ex3():
    W, L, R = 400, 20, 385
    X = lambda v: L + (v - 0) / 20 * (R - L)
    out = []
    rows = [(12, 3, BLUE, "แสงรวมที่กล้องวัดได้ T", "12 ± 3", 90),
            (5, 4, "#888", "พื้นหลังจากท้องฟ้าข้าง ๆ b", "5 ± 4", 190)]
    for m, s, col, name, val, base in rows:
        pts = [(X(v), base - 48 * gauss(v, m, s)) for v in np.linspace(0, 20, 201)]
        out.append(f'<path d="{path(pts)}" fill="none" stroke="{col}" stroke-width="3"/>')
        out.append(f'<line x1="{X(m):.1f}" y1="{base}" x2="{X(m):.1f}" y2="{base-48}" stroke="{col}" stroke-dasharray="3 3"/>')
        out.append(f'<text x="{L}" y="{base-62}" font-size="15" {FONT}>{name} = {val}</text>')
        out.append(f'<line x1="{L}" y1="{base}" x2="{R}" y2="{base}" stroke="currentColor" opacity=".5"/>')
    for v in (0, 5, 10, 15, 20):
        out.append(f'<text x="{X(v):.1f}" y="208" font-size="13" text-anchor="middle" {FONT}>{v}</text>')
    return svg(W, 215, "ระฆังสองลูก แสงรวมที่วัดได้ 12 บวกลบ 3 และพื้นหลัง 5 บวกลบ 4", "".join(out))


# ── peaks: Gaussian ยอดตรงกัน vs ตาราง ⓪ ยอดไม่ตรงกัน ──────────────────────
def fig_peaks():
    out = []
    # ซ้าย: ภูเขาของตัวอย่าง ② (T 10±1, b 3±2) · s ∈ [1,13] · b ∈ [-3,9]
    L, R, TOP, BOT = 15, 185, 40, 170
    X = lambda s: L + (s - 1) / 12 * (R - L)
    Y = lambda b: BOT - (b + 3) / 12 * (BOT - TOP)
    cov = np.array([[5.0, -4.0], [-4.0, 4.0]])          # (s, b) · var s = 1+4 · cov = −var b
    vals, vecs = np.linalg.eigh(cov)
    for k, op in ((math.sqrt(6.17), .25), (math.sqrt(2.30), .55)):   # 95% · 68% ของ 2 มิติ
        pts = []
        for t in np.linspace(0, 2 * math.pi, 120):
            v = vecs @ (k * np.sqrt(vals) * np.array([math.cos(t), math.sin(t)]))
            pts.append((X(7 + v[0]), Y(3 + v[1])))
        out.append(f'<path d="{path(pts)}Z" fill="{BLUE}" fill-opacity="{op}" stroke="none"/>')
    out.append(f'<circle cx="{X(7):.1f}" cy="{Y(3):.1f}" r="4" fill="currentColor"/>')
    out.append(f'<text x="{X(7)+8:.1f}" y="{Y(3)-6:.1f}" font-size="13" {FONT}>ยอดภูเขา</text>')
    base = 250
    pts = [(X(s), base - 55 * gauss(s, 7, math.sqrt(5))) for s in np.linspace(1, 13, 121)]
    out.append(f'<path d="{path(pts)}" fill="none" stroke="{GREEN}" stroke-width="3"/>')
    out.append(f'<line x1="{X(7):.1f}" y1="{Y(3):.1f}" x2="{X(7):.1f}" y2="{base}" stroke="currentColor" stroke-dasharray="4 3" opacity=".7"/>')
    out.append(f'<line x1="{L}" y1="{base}" x2="{R}" y2="{base}" stroke="currentColor" opacity=".5"/>')
    out.append(f'<text x="{X(7):.1f}" y="{base+17}" font-size="13" text-anchor="middle" {FONT}>s = 7</text>')
    out.append(f'<text x="{L}" y="20" font-size="15" font-weight="600" {FONT}>ใกล้ Gaussian</text>')
    out.append(f'<text x="{L}" y="{base+40}" font-size="13" {FONT}>ยอดภูเขากับยอดเงา</text>')
    out.append(f'<text x="{L}" y="{base+58}" font-size="13" {FONT}>อยู่ที่ s เดียวกัน</text>')

    # ขวา: ตาราง ⓪ · คอลัมน์ = s · แถว = b (b = 2 อยู่บน) · ผลรวมคอลัมน์ = แท่งข้างล่าง
    tab = {1: [0, 0, 28], 2: [14, 14, 14], 3: [10, 10, 10]}    # tab[s][b]
    x0, cw, ch = 245, 45, 40
    for i, s in enumerate((1, 2, 3)):
        for j, b in enumerate((2, 1, 0)):
            v = tab[s][b]
            x, y = x0 + i * cw, TOP + j * ch
            out.append(f'<rect x="{x}" y="{y}" width="{cw-3}" height="{ch-3}" fill="{BLUE}" fill-opacity="{0.08 + v/28*0.8:.2f}"/>')
            out.append(f'<text x="{x+cw/2-1.5}" y="{y+ch/2+4}" font-size="14" text-anchor="middle" {FONT}>{v}</text>')
        out.append(f'<text x="{x0+i*cw+cw/2-1.5}" y="{TOP+3*ch+14}" font-size="12" text-anchor="middle" {FONT}>s={s}</text>')
    for j, b in enumerate((2, 1, 0)):
        out.append(f'<text x="{x0-5}" y="{TOP+j*ch+ch/2+4}" font-size="12" text-anchor="end" {FONT}>b={b}</text>')
    out.append(f'<rect x="{x0-1}" y="{TOP-1}" width="{cw-1}" height="{ch-1}" fill="none" stroke="{RED}" stroke-width="2.5"/>')
    for i, s in enumerate((1, 2, 3)):
        tot = sum(tab[s]); h = tot / 42 * 55
        x = x0 + i * cw
        col = GREEN if s == 2 else "#999"
        out.append(f'<rect x="{x+6}" y="{base-h:.1f}" width="{cw-15}" height="{h:.1f}" fill="{col}"/>')
        out.append(f'<text x="{x+cw/2-1.5}" y="{base-h-4:.1f}" font-size="13" text-anchor="middle" {FONT}>{tot}</text>')
    out.append(f'<line x1="{x0}" y1="{base}" x2="{x0+3*cw}" y2="{base}" stroke="currentColor" opacity=".5"/>')
    out.append(f'<text x="{x0-25}" y="20" font-size="15" font-weight="600" {FONT}>ไม่ใช่ Gaussian (ตาราง ⓪)</text>')
    out.append(f'<text x="{x0-25}" y="{base+40}" font-size="13" {FONT}>ช่องสูงสุดอยู่ที่ s = 1</text>')
    out.append(f'<text x="{x0-25}" y="{base+58}" font-size="13" {FONT}>แต่เงาสูงสุดอยู่ที่ s = 2</text>')
    return svg(400, 318, "ซ้าย posterior รูประฆัง ยอดภูเขากับยอดเงาอยู่ที่ s เท่ากับ 7 ทั้งคู่ "
               "ขวา ตารางจากหัวข้อศูนย์ ช่องสูงสุดอยู่ที่ s เท่ากับ 1 แต่ผลรวมคอลัมน์สูงสุดอยู่ที่ s เท่ากับ 2",
               "".join(out))


# ── triangle plot ของตัวอย่างกาแล็กซี ─────────────────────────────────────────
def fig_triangle():
    out = []
    P, G = 150, 22                        # ขนาดช่อง · ช่องว่าง
    ox, oy = 50, 15
    sr, br = (0, 14), (-4, 10)
    cov = np.array([[5.0, -4.0], [-4.0, 4.0]]); mean = (7, 3)
    vals, vecs = np.linalg.eigh(cov)
    def box(cx, cy): return f'<rect x="{cx}" y="{cy}" width="{P}" height="{P}" fill="none" stroke="currentColor" opacity=".45"/>'
    # ช่องบนซ้าย p(s)
    X = lambda s, cx: cx + (s - sr[0]) / (sr[1] - sr[0]) * P
    out.append(box(ox, oy))
    pts = [(X(s, ox), oy + P - (P - 15) * gauss(s, 7, math.sqrt(5))) for s in np.linspace(*sr, 141)]
    out.append(f'<path d="{path(pts)}" fill="none" stroke="{GREEN}" stroke-width="2.5"/>')
    # ช่องล่างซ้าย p(s, b)
    cy = oy + P + G
    Yb = lambda b: cy + P - (b - br[0]) / (br[1] - br[0]) * P
    out.append(box(ox, cy))
    for k, op in ((math.sqrt(6.17), .3), (math.sqrt(2.30), .65)):
        pts = []
        for t in np.linspace(0, 2 * math.pi, 120):
            v = vecs @ (k * np.sqrt(vals) * np.array([math.cos(t), math.sin(t)]))
            pts.append((X(mean[0] + v[0], ox), Yb(mean[1] + v[1])))
        out.append(f'<path d="{path(pts)}Z" fill="{BLUE}" fill-opacity="{op}"/>')
    # ช่องล่างขวา p(b) — แกนนอนเป็น b
    cx2 = ox + P + G
    Xb = lambda b: cx2 + (b - br[0]) / (br[1] - br[0]) * P
    out.append(box(cx2, cy))
    pts = [(Xb(b), cy + P - (P - 15) * gauss(b, 3, 2)) for b in np.linspace(*br, 141)]
    out.append(f'<path d="{path(pts)}" fill="none" stroke="{GREEN}" stroke-width="2.5"/>')
    # ป้ายแกน
    for v in (0, 5, 10):
        out.append(f'<text x="{X(v, ox):.1f}" y="{cy+P+15}" font-size="12" text-anchor="middle" {FONT}>{v}</text>')
    for v in (0, 4, 8):
        out.append(f'<text x="{ox-5}" y="{Yb(v)+4:.1f}" font-size="12" text-anchor="end" {FONT}>{v}</text>')
        out.append(f'<text x="{Xb(v):.1f}" y="{cy+P+15}" font-size="12" text-anchor="middle" {FONT}>{v}</text>')
    out.append(f'<text x="{ox+P/2}" y="{cy+P+33}" font-size="14" text-anchor="middle" {FONT}>กาแล็กซี s</text>')
    out.append(f'<text x="{cx2+P/2}" y="{cy+P+33}" font-size="14" text-anchor="middle" {FONT}>พื้นหลัง b</text>')
    out.append(f'<text transform="translate(16,{cy+P/2}) rotate(-90)" font-size="14" text-anchor="middle" {FONT}>พื้นหลัง b</text>')
    # คำอธิบายในมุมขวาบนที่ว่าง
    tx = cx2 - 5
    for i, line in enumerate(["บนซ้าย: เงาบนแกน s", "ล่างซ้าย: ภูเขามองจากบน",
                              "  (สีเข้ม = ส่วนที่สูง)", "ล่างขวา: เงาบนแกน b"]):
        out.append(f'<text x="{tx}" y="{oy+40+i*22}" font-size="13" {FONT}>{line}</text>')
    return svg(400, 395, "triangle plot ของตัวอย่างกาแล็กซี ช่องบนซ้ายคือเงาบนแกน s ช่องล่างซ้ายคือภูเขาสองมิติที่เอียงลง "
               "ช่องล่างขวาคือเงาบนแกน b", "".join(out))


FIGS = {"ex3": fig_ex3, "peaks": fig_peaks, "triangle": fig_triangle}

if __name__ == "__main__":
    html = LESSON.read_text(encoding="utf-8")
    for k, f in FIGS.items():
        pat = re.compile(rf"(<!--FIG:{k}-->).*?(<!--/FIG:{k}-->)", re.S)
        assert pat.search(html), f"ไม่เจอ marker FIG:{k}"
        s = f()
        assert "$" not in s
        html = pat.sub(lambda m: m.group(1) + s + m.group(2), html)
    LESSON.write_text(html, encoding="utf-8")
    print("ฝัง", ", ".join(FIGS))
