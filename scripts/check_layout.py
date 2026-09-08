#!/usr/bin/env python3
"""ประตูตรวจการแสดงผลจริง — เปิดทุกบทเรียนในเบราว์เซอร์จริง

`check_lessons.py` อ่านไฟล์อย่างเดียว จับ layout ไม่ได้ · ตัวนี้เปิด Chromium
จริงแล้ววัด:

  1. หน้าต้องไม่ล้นแนวนอนที่ 375px และ 1200px
     (ตารางกว้างต้องเลื่อนในกล่องตัวเอง ไม่ใช่ดันทั้งหน้า — เกิดจริงในบทที่ 1)
  2. ไม่มี JavaScript error
  3. quiz ทุกอันคลิกแล้วทำงาน: ปุ่มล็อก · เฉลยแสดง · ติดคลาสถูก/ผิดตรงตาม data-answer

⚠️ ตรวจการเรนเดอร์สมการไม่ได้ — CDN ถูก proxy บล็อกในคลาวด์ (ดู CLAUDE.md)
   ต้องเปิดไฟล์บนเครื่องยืนยันอีกชั้น

    pip install playwright        # ครั้งเดียว · อย่ารัน playwright install
    python3 scripts/check_layout.py
"""
import pathlib, sys

try:
    from playwright.sync_api import sync_playwright
except ImportError:
    sys.exit("ต้องติดตั้งก่อน: pip install playwright")

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHROME = "/opt/pw-browsers/chromium"
problems = []

files = sorted((ROOT / "lessons").glob("*.html")) + sorted((ROOT / "reference").glob("*.html"))
if not files:
    sys.exit("ยังไม่มีบทเรียนให้ตรวจ")

with sync_playwright() as p:
    browser = p.chromium.launch(executable_path=CHROME)
    for f in files:
        rel = f.relative_to(ROOT)
        for width in (375, 1200):
            errs = []
            page = browser.new_page(viewport={"width": width, "height": 900})
            page.on("pageerror", lambda e: errs.append(str(e)))
            page.goto(f.resolve().as_uri())
            page.wait_for_timeout(1200)

            if page.evaluate("document.documentElement.scrollWidth"
                             " > document.documentElement.clientWidth + 1"):
                problems.append(f"{rel} @{width}px: ล้นแนวนอน")

            if width == 1200:                       # ทดสอบ quiz ครั้งเดียวพอ
                for i in range(page.locator(".quiz").count()):
                    quiz = page.locator(".quiz").nth(i)

                    if quiz.get_attribute("data-open") is not None:
                        # แบบพิมพ์เอง: ยิงเฉลยของตัวเองเข้าไป ต้องได้ right
                        # (จับกรณีเฉลยที่ normalizer อ่านไม่ตรงกับตัวเอง — เกิดได้จริง
                        #  เมื่อเฉลยมีวงเล็บหรือเครื่องหมายคูณที่ถูก normalize ทิ้ง)
                        key = quiz.get_attribute("data-answer")
                        quiz.locator(".ans input").fill(key)
                        quiz.locator(".ans button").click()
                        page.wait_for_timeout(120)
                        st = quiz.evaluate("q => ({res: q.dataset.result,"
                                           " marked: q.querySelector('.ans input').classList.contains('right'),"
                                           " fb: !!q.querySelector('.fb.show'),"
                                           " locked: q.querySelector('.ans button').disabled})")
                        if not (st["res"] == "right" and st["marked"] and st["fb"] and st["locked"]):
                            problems.append(f"{rel}: quiz[data-open] #{i+1} เฉลยของตัวเองไม่ผ่าน — {st}")
                        continue

                    ans = int(quiz.get_attribute("data-answer"))
                    n = quiz.locator("button").count()
                    pick = (ans + 1) % n            # จงใจตอบผิด
                    quiz.locator("button").nth(pick).click()
                    page.wait_for_timeout(120)
                    st = quiz.evaluate("""q => ({
                        res: q.dataset.result,
                        right: q.querySelectorAll('button')[+q.dataset.answer].classList.contains('right'),
                        fb: !!q.querySelector('.fb.show'),
                        locked: [...q.querySelectorAll('button')].every(b => b.disabled)})""")
                    if not (st["res"] == "wrong" and st["right"] and st["fb"] and st["locked"]):
                        problems.append(f"{rel}: quiz #{i+1} ไม่ทำงานถูก — {st}")

            for e in errs:
                problems.append(f"{rel} @{width}px: JS error — {e[:90]}")
            page.close()
    browser.close()

if problems:
    print(f"❌ เจอ {len(problems)} จุด\n")
    for x in problems:
        print("  •", x)
    sys.exit(1)

print(f"✅ ผ่าน — เปิดจริง {len(files)} ไฟล์ × 2 ความกว้าง")
