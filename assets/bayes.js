/* bayes.js — ของเล่น Bayes แบบกล่องสองใบ (JavaScript ล้วน ไม่มี library)
   ใช้ซ้ำได้ทุกบทในเส้น thesis · แค่วาง markup ไม่ต้องเขียน JS เพิ่ม

   ── เกมกล่องลับ ─────────────────────────────────────────
     <div class="bayes-game" data-ra="0.3" data-rb="0.7" data-pb="0.5"></div>
     ra / rb = สัดส่วนลูกแดงในกล่อง A / B · pb = สัดส่วนกล่อง B ในถุง (คอมพิวเตอร์สุ่มหยิบกล่องจากถุง 1 ใบ)

   ── ตาราง 100 รอบ ──────────────────────────────────────
     <div class="bayes-grid" data-ra="0.3" data-rb="0.7" data-pb="0.5"></div>
     มีแถบเลื่อนปรับทั้งสามค่า (ทีละ 0.1 ⇒ ทุกจำนวนรอบเป็นจำนวนเต็มพอดี ไม่ต้องปัด)

   ⚠️ build_standalone.py ฝังไฟล์นี้ลงฉบับ build/ ให้เอง (ทุกสคริปต์ที่ชี้ไป ../assets/*.js) */
(function () {
  "use strict";
  const pct = x => Math.round(x * 100) + "%";
  const el = (tag, cls, txt) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (txt != null) e.textContent = txt;
    return e;
  };
  const num = (root, k, d) => {
    const v = parseFloat(root.dataset[k]);
    return isNaN(v) ? d : v;
  };

  /* ---------- เกมกล่องลับ ---------- */
  function game(root) {
    const ra = num(root, "ra", 0.3), rb = num(root, "rb", 0.7), pb = num(root, "pb", 0.5);
    let secret, drawn, pB, wins = 0, plays = 0, totalDraws = 0, over;

    const boxes = el("div", "bg-boxes");
    [["A", ra], ["B", rb]].forEach(([n, r]) => {
      const b = el("div", "bg-box bg-box-" + n);
      b.append(el("strong", null, "กล่อง " + n));
      const dots = el("div", "bg-dots");
      for (let i = 0; i < 10; i++) dots.append(el("span", i < Math.round(r * 10) ? "ball red" : "ball white"));
      b.append(dots, el("small", null, "แดง " + Math.round(r * 10) + " ใน 10"));
      boxes.append(b);
    });
    const prompt = el("p", "bg-prompt");
    const row = el("div", "bg-drawn");
    const meterWrap = el("div", "bg-meter-wrap");
    const meter = el("div", "bg-meter");
    const mA = el("div", "bg-mA"), mB = el("div", "bg-mB");
    meter.append(mA, mB);
    const toggle = el("label", "bg-toggle");
    const cb = el("input"); cb.type = "checkbox";
    toggle.append(cb, document.createTextNode(" เปิดเครื่องวัดความเชื่อ (Bayes คำนวณให้หลังหยิบแต่ละลูก)"));
    meterWrap.append(toggle, meter);
    const btns = el("div", "bg-btns");
    const bDraw = el("button", "bg-draw", "🎲 หยิบลูก (ใส่คืนทุกครั้ง)");
    const gA = el("button", null, "ทายว่า A"), gB = el("button", null, "ทายว่า B");
    const bNew = el("button", "bg-new", "↻ เกมใหม่");
    [bDraw, gA, gB, bNew].forEach(b => b.type = "button");
    btns.append(bDraw, gA, gB, bNew);
    const result = el("p", "bg-result");
    const score = el("p", "bg-score");
    root.append(boxes, prompt, row, btns, meterWrap, result, score);

    function paint() {
      mA.style.width = pct(1 - pB); mB.style.width = pct(pB);
      mA.textContent = "A " + pct(1 - pB); mB.textContent = "B " + pct(pB);
      meter.hidden = !cb.checked;
    }
    function reset() {
      secret = Math.random() < pb ? "B" : "A";
      drawn = 0; pB = pb; over = false;
      row.textContent = ""; result.textContent = "";
      prompt.textContent = "คอมพิวเตอร์สุ่มหยิบกล่องจากถุงไว้แล้ว 1 ใบ แต่ไม่บอกว่าใบไหน · หยิบลูกดูก่อน แล้วค่อยทาย";
      bDraw.disabled = gA.disabled = gB.disabled = false;
      boxes.classList.remove("reveal-A", "reveal-B");
      paint();
    }
    bDraw.onclick = () => {
      if (over) return;
      const r = secret === "A" ? ra : rb;
      const red = Math.random() < r;
      drawn++;
      row.append(el("span", red ? "ball red" : "ball white"));
      // Bayes หนึ่งขั้น: posterior ของลูกก่อน = prior ของลูกนี้
      const lB = red ? rb : 1 - rb, lA = red ? ra : 1 - ra;
      pB = lB * pB / (lB * pB + lA * (1 - pB));
      prompt.textContent = "หยิบไปแล้ว " + drawn + " ลูก · พร้อมทายหรือยัง";
      paint();
    };
    function guess(g) {
      if (over) return;
      over = true; plays++; totalDraws += drawn;
      const ok = g === secret;
      if (ok) wins++;
      boxes.classList.add("reveal-" + secret);
      result.textContent = (ok ? "✅ ถูก" : "❌ ผิด") + " — กล่องลับคือ " + secret +
        " · ตอนทาย Bayes ให้ B " + pct(pB) + " (หยิบ " + drawn + " ลูก)";
      score.textContent = "ทายถูก " + wins + " จาก " + plays + " เกม · หยิบเฉลี่ย " +
        (totalDraws / plays).toFixed(1) + " ลูกต่อเกม";
      bDraw.disabled = gA.disabled = gB.disabled = true;
    }
    gA.onclick = () => guess("A");
    gB.onclick = () => guess("B");
    bNew.onclick = reset;
    cb.onchange = paint;
    reset();
  }

  /* ---------- ตาราง 100 รอบ ---------- */
  function grid(root) {
    const st = { pb: num(root, "pb", 0.5), ra: num(root, "ra", 0.3), rb: num(root, "rb", 0.7), seen: false };
    const ctrls = el("div", "bgr-ctrls");
    const mk = (key, label) => {
      const w = el("label", "bgr-slider");
      const s = el("input"); s.type = "range"; s.min = 0; s.max = 10; s.step = 1;
      s.value = Math.round(st[key] * 10);
      const out = el("output");
      s.oninput = () => { st[key] = s.value / 10; draw(); };
      w.append(el("span", null, label), s, out);
      ctrls.append(w);
      return out;
    };
    const oPB = mk("pb", "กล่อง B ในถุง"), oRA = mk("ra", "ลูกแดงในกล่อง A"), oRB = mk("rb", "ลูกแดงในกล่อง B");
    const bSeen = el("button", "bgr-seen"); bSeen.type = "button";
    bSeen.onclick = () => { st.seen = !st.seen; draw(); };
    const field = el("div", "bgr-field");
    const out = el("div", "bgr-out");
    root.append(ctrls, bSeen, field, out);

    function draw() {
      // ทีละ 0.1 ⇒ nB เป็นพหุคูณของ 10 ⇒ nB·rb เป็นจำนวนเต็มพอดี
      const nB = Math.round(st.pb * 100), nA = 100 - nB;
      const redA = Math.round(nA * st.ra), redB = Math.round(nB * st.rb), red = redA + redB;
      // แสดงเป็น "กี่ใบ/กี่ลูก จาก 10" — นับได้ด้วยตา ไม่ต้องแปลทศนิยม
      oPB.textContent = Math.round(st.pb * 10) + " / 10"; oRA.textContent = Math.round(st.ra * 10) + " / 10";
      oRB.textContent = Math.round(st.rb * 10) + " / 10";
      bSeen.textContent = st.seen ? "↺ ดูทั้ง 100 รอบอีกครั้ง" : "👁 สมมติว่าหยิบได้แดง — ตัดรอบที่ไม่แดงทิ้ง";
      field.textContent = "";
      [["A", nA, redA], ["B", nB, redB]].forEach(([n, tot, r]) => {
        const g = el("div", "bgr-group bgr-" + n);
        g.append(el("div", "bgr-label", "กล่อง " + n + " · " + tot + " รอบ"));
        const dots = el("div", "bgr-dots");
        for (let i = 0; i < tot; i++) {
          const isRed = i < r;
          dots.append(el("span", "ball " + (isRed ? "red" : "white") + (st.seen && !isRed ? " gone" : "")));
        }
        g.append(dots);
        if (tot) field.append(g);
      });
      out.textContent = "";
      const line = t => out.append(el("p", null, t));
      line("แดงจากกล่อง A: " + nA + " × " + st.ra.toFixed(1) + " = " + redA + " รอบ");
      line("แดงจากกล่อง B: " + nB + " × " + st.rb.toFixed(1) + " = " + redB + " รอบ");
      line("แดงทั้งหมด: " + redA + " + " + redB + " = " + red + " รอบ");
      const p = el("p", "bgr-answer");
      p.textContent = red
        ? "ในรอบที่ได้แดง เป็นกล่อง B " + redB + " / " + red + " = " + (redB / red).toFixed(3).replace(/\.?0+$/, "")
        : "ไม่มีรอบไหนได้แดงเลย — หยิบได้แดงเป็นไปไม่ได้";
      out.append(p);
    }
    draw();
  }

  function init() {
    document.querySelectorAll(".bayes-game").forEach(game);
    document.querySelectorAll(".bayes-grid").forEach(grid);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
