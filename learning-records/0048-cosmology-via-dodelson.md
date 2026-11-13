# 0048 · อาจารย์สั่งอ่าน Dodelson ⇒ cosmology (ช่วง G2) ขึ้นมาก่อน · สร้างบทที่ 28

**ที่มา (2026-11-13):** แบแจ้งว่าอาจารย์ให้อ่าน S. Dodelson, *Modern Cosmology* (2003) "เพื่อเข้าใจพื้นฐานให้แน่น"
และบอกเองว่า**ไม่ชำนาญภาษาอังกฤษ** — ขอให้แปลเป็นไทย

## การตัดสินใจ

- ไม่ทำฉบับแปลทั้งเล่ม (หนังสือมีลิขสิทธิ์ · repo อยู่บน GitHub) ⇒ ทำ **คู่มืออ่าน** รายบทแทน
  (`reference/dodelson-chNN.html` · สเปกที่ `reference/DODELSON-FORMAT.md` · แผนที่ `reference/dodelson-index.html`)
  มีตารางศัพท์+สำนวนอังกฤษทุกบท เพราะอุปสรรคที่แบบอกคือภาษา ไม่ใช่คณิต
- MISSION ไม่เปลี่ยน (อ่าน paper hep-th/gr-qc) · ที่เปลี่ยนคือ**ลำดับ**: G2 มาก่อน C–F ที่เหลือ
  และ GR ปูเฉพาะที่ FRW ต้องใช้ ไปพร้อมกับ Dodelson บทที่ 2 (ทางเลือกที่ ROADMAP เปิดไว้อยู่แล้ว)
- บทที่ 27 (time ordering) ยังรอผลรอบแรก — ไม่ทิ้ง กลับมาเก็บได้

## บทที่ 28 (`lessons/0041-expanding-universe.html`) — Dodelson §1.1–1.2

① scale factor · ② comoving $x$ กับ physical $d = ax$ สมการ (1) · ③ $H = \dot a/a$ สมการ (2) · ตัวอย่าง $a = t^{2/3} \Rightarrow H = 2/(3t)$, $Ht = 2/3$ ·
④ $v = \dot d = \dot a x = Hd$ สมการ (3) · Hubble diagram ความชัน $= H$ · ⑤ $1 + z = \lambda_{obs}/\lambda_{emit} = 1/a$ สมการ (4)

ทุกสมการที่ตั้งเลขมีข้อฝึกที่ผู้เรียน*ผลิต*ค่าเอง (6 ข้อ 19 ขั้น) · ข้อ 6 ทวนการนับมิติของบทที่ 9

ยังไม่แตะโดยตั้งใจ: สมการ Friedmann · $\rho \propto a^{-3}, a^{-4}$ · critical density · เหตุที่แสงยืดตามกริด

`scripts/verify/0042_expanding_universe.py` 21/21 (sympy: $\dot d = Hd$ สำหรับ $a(t)$ ใด ๆ)

ของที่แก้ระหว่างทาง: `quiz.js` tolerance 0.5% ไม่มีเพดาน ⇒ คีย์ ≥ 200 รับค่าที่คลาดไปหนึ่ง (500 รับ 501) — ใส่เพดาน 0.5
· `check_layout.py` ผูก path Chromium ของคลาวด์ ⇒ รันบนเครื่องแบไม่ได้ — ให้หา Chrome/Edge ในเครื่องด้วย

## ผลรอบแรก

(รอผล)
