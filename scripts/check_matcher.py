#!/usr/bin/env python3
"""ประตูที่ 4 — ตัวตรวจคำตอบรับ/ไม่รับอะไรบ้าง

ประตูอีกสามตัวดูว่า "หน้าเปิดได้ไหม" ไม่ได้ดูว่า **กดตรวจแล้วได้ผลถูกไหม**
ซึ่งเป็นจุดที่พังเงียบที่สุด: คำตอบถูกแต่ถูกนับผิด ผู้เรียนไม่มีทางรู้ว่า
ตัวเองถูกหรือตัวตรวจพัง (เคยเจอจริง — num() ไม่รับ "+1" ในบท 0017)

สคริปต์นี้จำลอง norm()/num()/judge() ของ assets/quiz.js ในภาษา Python
แล้วยิงสองฝั่งกับทุก data-answer ในทุกบท:

  ฝั่งรับ  — คีย์เอง · alt ทุกตัว · รูปแบบที่คนพิมพ์จริง (เว้นวรรค · ขีดยูนิโคด
            · เครื่องหมายบวกนำหน้า · ช่องว่างรอบ *) ต้องได้ 'right'
  ฝั่งปฏิเสธ — ค่าที่ผิดแบบใกล้เคียงที่สุด (เครื่องหมายกลับ · คลาดไปหนึ่ง)
            ต้องได้ 'wrong' · ตัวตรวจที่รับทุกอย่างคือตัวตรวจที่ไม่ได้ตรวจอะไร

⚠️ ห้ามแก้ที่คีย์เพื่อให้ประตูผ่าน — ถ้าฝั่งรับตก แปลว่า quiz.js ต้องแก้
"""
import re, sys, html
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


# ── จำลอง quiz.js ─────────────────────────────────────────
SUP = str.maketrans("\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079",
                    "0123456789")


def norm(v: str) -> str:
    """ต้องตรงกับ norm() ใน assets/quiz.js เป๊ะ — นี่คือเหตุผลทั้งหมดที่ไฟล์นี้มีอยู่"""
    s = str(v).strip().lower()
    s = re.sub(r"[−–—]", "-", s)
    s = s.translate(SUP)                    # ⁰¹²³… → เลขธรรมดา
    s = s.replace("**", "^")
    s = re.sub(r"[\s*·]", "", s)
    # ⚠️ ตัด ^ ทิ้ง — ก๊อป $m^2$ จากหน้าที่เรนเดอร์แล้วได้ "m 2" ไม่มี ^
    #    ไม่ตัดจะปฏิเสธคำตอบที่ถูก (เกิดจริง 2026-09-24)
    s = s.replace("^", "")
    return re.sub(r"[()]", "", s)


def num(v: str):
    s = norm(v)
    m = re.match(r"^([+-]?\d*\.?\d+)/([+-]?\d*\.?\d+)$", s)
    if m:
        return float(m.group(1)) / float(m.group(2))
    return float(s) if re.match(r"^[+-]?\d*\.?\d+$", s) else None


def judge(keys, typed: str) -> str:
    if not typed.strip():
        return "empty"
    got = num(typed)
    for k in keys:
        want = num(k)
        if got is not None and want is not None:
            tol = max(1e-9, abs(want) * 5e-3)
            if abs(got - want) <= tol:
                return "right"
        elif norm(typed) == norm(k):
            return "right"
    return "wrong"


# ── สร้างเคสจากคีย์ ───────────────────────────────────────
def accept_forms(keys):
    """รูปแบบที่ผู้เรียนพิมพ์จริงแล้วต้องนับว่าถูก"""
    out = []
    for k in keys:
        out += [k, " " + k + " ", k.upper()]
        n = num(k)
        if n is not None and n > 0 and not k.startswith("+"):
            out.append("+" + k)
        if k.startswith("-"):
            out.append("−" + k[1:])          # ขีดยูนิโคด
        # ⚠️ คีย์ที่เขียน ** อยู่แล้วห้ามแปลงเรื่อง * — จะกลายเป็นสตริงที่ไม่มีใครพิมพ์
        if "*" in k and "**" not in k:
            out.append(k.replace("*", " * "))
            out.append(k.replace("*", ""))        # norm() ตัด * ทิ้งอยู่แล้ว
        if "^" in k:
            out.append(k.replace("^", "**"))
    return out


def reject_forms(keys):
    """ค่าที่ผิดแบบใกล้เคียงที่สุด — ต้องไม่ถูกรับ"""
    out = set()
    main = keys[0]
    n = num(main)
    if n is not None:
        out.add(str(n + 1))
        if abs(n) > 1e-12:
            out.add(str(-n))                      # เครื่องหมายกลับ
    else:
        out.add(main[1:] if main.startswith("-") else "-" + main)
    # ตัดตัวที่ดันไปตรงกับ alt ที่ตั้งใจรับ (ไม่มีในบทไหนตอนนี้ แต่กันไว้)
    return [r for r in out if not any(norm(r) == norm(k) for k in keys)]


# ── เคสที่เขียนมือ — ความเข้าใจผิดที่เจอจริงของบทนั้น ─────
# ตัวสร้างอัตโนมัติจับได้แค่ "เครื่องหมายกลับ" กับ "คลาดไปหนึ่ง"
# ส่วนนี้ไว้ใส่ค่าที่ผู้เรียนตอบผิดมาจริง ๆ แล้วต้องไม่ถูกรับ
HAND = [
    # (ไฟล์, คีย์, สิ่งที่พิมพ์, ผลที่ต้องได้)
    ("0019-current-T.html", "-px", "px", "wrong"),          # ลืมยกดัชนี (g มีลบ)
    ("0019-current-T.html", "-m^2*phi", "m^2*phi", "wrong"),  # ลืมลบจาก EL
    ("0019-current-T.html", "9", "0", "wrong"),             # นึกว่า T00 = -1/4 F^2
    ("0019-current-T.html", "9", "18", "wrong"),            # ลืม 1/2
    ("0019-current-T.html", "4", "0", "wrong"),             # [T] ไม่ใช่ไร้มิติ
    ("0016-fields.html", "-m^2*phi", "m^2*phi", "wrong"),

    # ── บทที่ 9 · ความเข้าใจผิดที่เกิดจริงในผล 2/8 ของบทที่ 8 ──────────
    ("0020-divergence.html", "ptt-pxx", "pt-px", "wrong"),    # ลืมอนุพันธ์ตัวนอก (คำตอบจริงของแบ)
    ("0020-divergence.html", "ptt-pxx", "ptt+pxx", "wrong"),  # เอาลบของ j^x ไปรวมกับตรงกลาง
    ("0020-divergence.html", "ptt+pxx", "ptt-pxx", "wrong"),  # กลับกัน
    ("0020-divergence.html", "-m^2*phi", "0", "wrong"),       # คำท่อง "อนุรักษ์ = ศูนย์"
    ("0020-divergence.html", "-m^2*phi", "m^2*phi", "wrong"), # ลืมลบตอนย้ายข้าง
    ("0020-divergence.html", "2*ptx", "ptx", "wrong"),        # ลืมว่าสองก้อนเป็นตัวเดียวกัน
    ("0020-divergence.html", "50", "10", "wrong"),            # ลืมยกกำลังสอง m (คำตอบจริงของแบ)
    ("0020-divergence.html", "19", "20", "wrong"),            # ลืม 1/2 ของพจน์จลน์ (คำตอบจริงของแบ)
    ("0020-divergence.html", "1", "2", "wrong"),              # ลืม 1/2 ทั้งสองพจน์
    ("0020-divergence.html", "18", "36", "wrong"),            # ลืม 1/2 ของพจน์มวล
    ("0020-divergence.html", "8", "16", "wrong"),
    ("0020-divergence.html", "-4*phi", "4*phi", "wrong"),
    ("0020-divergence.html", "pt", "phi", "wrong"),           # ไม่ได้อนุพันธ์

    # ── บทที่ 10 · ความเข้าใจผิดที่เกิดจริงในผล 5/8 ของบทที่ 9 ─────────
    # สองเคสแรกคือ "สลับช่อง off-shell กับ on-shell" ซึ่งเป็นทั้งหมดของบทที่ 10
    ("0021-onshell.html", "0", "ptt-pxx", "wrong"),           # เอาค่า off-shell ใส่ช่อง on-shell
    ("0021-onshell.html", "ptt-pxx", "0", "wrong"),           # เอาค่า on-shell ใส่ช่อง off-shell
    ("0021-onshell.html", "ptt-pxx", "Ptt-pxx", "right"),     # ⚠️ ตัวใหญ่ต้องรับ (พิมพ์มาจริง)
    ("0021-onshell.html", "-px", "px", "wrong"),              # เดาเครื่องหมายจากแถวบน
    ("0021-onshell.html", "-m^2*phi", "m^2*phi", "wrong"),
    ("0021-onshell.html", "2", "0", "wrong"),                 # นึกว่า off-shell เป็นศูนย์เสมอ
    ("0021-onshell.html", "6", "5", "wrong"),
    ("0021-onshell.html", "20", "0", "wrong"),                # นึกว่า off-shell เป็นศูนย์เสมอ

    # ── ⚠️ ฝั่งรับ: ก๊อปสมการที่เรนเดอร์แล้วมาวาง ────────────────────
    # $m^2$ บนหน้าเว็บ ก๊อปมาได้ "m 2" (ไม่มี ^) หรือ "m²" (ตัวยกยูนิโคด)
    # ทั้งสองแบบเป็นคำตอบที่ถูก ⇒ ตัวตรวจต้องรับ
    # (เกิดจริง 2026-09-24: ตอบ "pt · ( ptt − pxx + m 2 phi )" ถูกแต่ถูกปฏิเสธ)
    ("0021-onshell.html", "-m^2*phi", "-m 2 phi", "right"),
    ("0021-onshell.html", "-m^2*phi", "-m²*phi", "right"),
    ("0021-onshell.html", "-m^2*phi", "−m 2 phi", "right"),
    ("0020-divergence.html", "-m^2*phi", "-m 2 phi", "right"),
    ("0019-current-T.html", "-m^2*phi", "-m²phi", "right"),
    ("0016-fields.html", "m^2", "m²", "right"),
    # …แต่ยังต้องปฏิเสธเครื่องหมายที่ผิดเหมือนเดิม
    ("0021-onshell.html", "-m^2*phi", "m 2 phi", "wrong"),

    # ── บทที่ 11 · สนามเชิงซ้อน — ความพลาดที่เดาได้จากโครงของบท ──────
    ("0022-u1charge.html", "-i*phis", "i*phis", "wrong"),     # ลืมสังยุคของ i
    ("0022-u1charge.html", "i*phi", "phi", "wrong"),          # ลืม i
    ("0022-u1charge.html", "i*phi", "i*phis", "wrong"),       # สลับสนามกับสังยุค
    ("0022-u1charge.html", "pts", "pt", "wrong"),             # อนุพันธ์เทียบผิดตัว
    ("0022-u1charge.html", "-pxs", "pxs", "wrong"),           # ลืมลบจากเมตริก
    ("0022-u1charge.html", "-2i", "2i", "wrong"),
    ("0022-u1charge.html", "-2", "2", "wrong"),               # ลบผิดทาง (6-4 แทน 4-6)
    ("0022-u1charge.html", "-2i", "-2", "wrong"),             # ลืม i
    ("0022-u1charge.html", "0", "30", "wrong"),               # บวกแทนลบในข้อสนามจริง
    ("0022-u1charge.html", "12", "16", "wrong"),              # 2*7+1*2
    # ฝั่งรับ: เขียนสลับที่ได้ (การคูณสลับที่)
    ("0022-u1charge.html", "i*phi*pts", "i*pts*phi", "right"),
    ("0022-u1charge.html", "i*phi*pts", "phi*pts*i", "right"),
    ("0022-u1charge.html", "-i*phis*pt", "−i*phis*pt", "right"),
    ("0022-u1charge.html", "-2i", "-2*i", "right"),

    # ── บทที่ 12 · gauge — ความพลาดที่เดาได้จากโครงของบท ────────────
    ("0023-gauge.html", "i*dta*phi", "i*alpha*pt", "wrong"),  # ตอบพจน์ปกติแทนพจน์เกิน
    ("0023-gauge.html", "i*dta*phi", "i*alpha*phi", "wrong"), # ลืมว่าอนุพันธ์ไปโดน alpha
    ("0023-gauge.html", "-i*dta*phi", "i*dta*phi", "wrong"),  # เครื่องหมายของชิ้น A
    ("0023-gauge.html", "i*alpha", "i*alpha*phi", "wrong"),   # ตอบทั้ง delta phi แทนตัวคูณ
    ("0023-gauge.html", "dta/e", "e*dta", "wrong"),           # คูณแทนหาร
    ("0023-gauge.html", "dta/e", "dta", "wrong"),             # ลืม 1/e
    ("0023-gauge.html", "-e", "e", "wrong"),                  # เครื่องหมายของการคู่ควบ
    ("0023-gauge.html", "-e", "-i*e", "wrong"),               # ลืมว่า i ไปรวมเป็น j
    ("0023-gauge.html", "3", "12", "wrong"),                  # e*dta แทน dta/e
    ("0023-gauge.html", "-30", "30", "wrong"),
    ("0023-gauge.html", "-11", "-8", "wrong"),                # คิดแค่ช่อง t ลืมช่อง x
    ("0023-gauge.html", "12", "6", "wrong"),                  # ลืมยกกำลัง e
    # ฝั่งรับ
    ("0023-gauge.html", "i*alpha*pt", "i*α*pt", "right"),
    ("0023-gauge.html", "dta/e", "(1/e)*dta", "right"),
    ("0023-gauge.html", "pts*pt", "pt*pts", "right"),
]


def run_hand(bad):
    for fname, key, typed, want in HAND:
        f = ROOT / "lessons" / fname
        src = f.read_text(encoding="utf-8")
        hit = None
        for tag in re.findall(r"<(?:div|li)[^>]*data-answer=[^>]*>", src):
            ks = keys_of(tag)
            if ks and ks[0] == key:
                hit = ks
                break
        if hit is None:
            bad.append(f"{fname}: เคสมือชี้คีย์ {key!r} ซึ่งไม่มีในไฟล์แล้ว")
            continue
        got = judge(hit, typed)
        if got != want:
            bad.append(f"{fname}: คีย์ {key!r} + {typed!r} ได้ {got!r} (ควร {want!r})")


def keys_of(tag: str):
    a = re.search(r'data-answer="([^"]*)"', tag)
    if not a:
        return None
    ks = [html.unescape(a.group(1))]
    alt = re.search(r'data-alt="([^"]*)"', tag)
    if alt:
        ks += [x for x in html.unescape(alt.group(1)).split("|") if x]
    return ks


def main():
    files = sorted((ROOT / "lessons").glob("*.html"))
    bad, n_keys, n_cases = [], 0, 0

    for f in files:
        src = f.read_text(encoding="utf-8")
        # ทุก element ที่มี data-answer (ทั้ง .quiz แบบพิมพ์เอง และ <li> ของขั้น)
        for tag in re.findall(r"<(?:div|li)[^>]*data-answer=[^>]*>", src):
            if "data-answer" not in tag:
                continue
            ks = keys_of(tag)
            # ปรนัยใช้ data-answer เป็นดัชนีตัวเลือก ไม่ได้ผ่าน judge()
            if "data-open" not in tag and "<li" not in tag and "data-steps" not in src[: src.index(tag)] :
                pass
            if re.match(r"^<div", tag) and "data-open" not in tag:
                continue                          # ปรนัย — ข้าม
            n_keys += 1
            for t in accept_forms(ks):
                n_cases += 1
                if judge(ks, t) != "right":
                    bad.append(f"{f.name}: คีย์ {ks[0]!r} ไม่รับ {t!r} (ควรถูก)")
            for t in reject_forms(ks):
                n_cases += 1
                if judge(ks, t) == "right":
                    bad.append(f"{f.name}: คีย์ {ks[0]!r} ดันรับ {t!r} (ควรผิด)")

    run_hand(bad)

    if bad:
        print("❌ ตัวตรวจคำตอบมีปัญหา")
        for b in bad:
            print("   " + b)
        sys.exit(1)
    print(f"✅ ผ่าน — {n_keys} คีย์ · {n_cases} เคส ({len(files)} ไฟล์)")


if __name__ == "__main__":
    main()
