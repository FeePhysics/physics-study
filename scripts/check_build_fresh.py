#!/usr/bin/env python3
"""ประตู: ไฟล์ใน build/ ต้องตรงกับ lessons/ ที่เป็นต้นฉบับ

ทำไมต้องมี — ตัวที่ถูก **publish** คือ `build/artifact-*.html` ไม่ใช่ `lessons/*.html`
ส่วนประตูอีกสองตัวอ่าน `lessons/` ทั้งคู่ ⇒ **build/ ค้างเก่าได้โดยไม่มีอะไรฟ้อง**
แล้วทั้งสองประตูจะผ่านฉลุยในขณะที่ของที่ส่งให้ผู้เรียนเป็นฉบับก่อนแก้

เกิดจริงมาแล้ว: รัน `build_standalone.py` เปล่า ๆ (ไม่ใส่ `--all`) มันพิมพ์ docstring
แล้ว exit ซึ่ง **หน้าตาเหมือนรันสำเร็จ** · ประตู layout ผ่าน · เกือบ publish ของเก่า

วิธีตรวจ: สร้างใหม่ลงที่ชั่วคราวแล้วเทียบทีละไบต์ (ตัวสร้างเป็นฟังก์ชันบริสุทธิ์
⇒ ต่างเมื่อไหร่แปลว่า build/ ค้าง)

    python3 scripts/check_build_fresh.py
"""
import pathlib, shutil, subprocess, sys, tempfile, filecmp

ROOT = pathlib.Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"
srcs = sorted((ROOT/"lessons").glob("*.html")) + sorted((ROOT/"reference").glob("*.html"))

if not BUILD.exists():
    sys.exit("ยังไม่เคยสร้าง build/ — รัน: python3 scripts/build_standalone.py --all")

stale, missing = [], []
with tempfile.TemporaryDirectory() as tmp:
    tmp = pathlib.Path(tmp)
    # สร้างใหม่ลงที่ชั่วคราว โดยสลับ build/ ออกไปก่อน แล้วคืนที่เดิมเสมอ
    keep = ROOT / "build.checking"
    BUILD.rename(keep)
    try:
        subprocess.run([sys.executable, "scripts/build_standalone.py", "--all"],
                       cwd=ROOT, check=True, capture_output=True)
        fresh = BUILD
        for s in srcs:
            for name in (s.name, "artifact-" + s.name):
                a, b = keep / name, fresh / name
                if not a.exists():
                    missing.append(name)
                elif not filecmp.cmp(a, b, shallow=False):
                    stale.append(name)
        shutil.rmtree(BUILD)          # ของสดถูกเขียนทับตอนท้ายอยู่แล้ว
    finally:
        keep.rename(BUILD)

if missing or stale:
    for n in missing: print(f"✗ ไม่มีใน build/: {n}")
    for n in stale:   print(f"✗ ค้างเก่า:      {n}")
    sys.exit("\nbuild/ ไม่ตรงกับต้นฉบับ — รัน: python3 scripts/build_standalone.py --all")

print(f"✅ ผ่าน — build/ ตรงกับต้นฉบับครบ {len(srcs)} ไฟล์")
