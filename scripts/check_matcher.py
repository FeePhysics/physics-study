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
def norm(v: str) -> str:
    s = str(v).strip().lower()
    s = re.sub(r"[−–—]", "-", s)
    s = s.replace("**", "^")
    s = re.sub(r"[\s*·]", "", s)
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
