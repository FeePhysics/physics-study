/* bayes.js — ของเล่น Bayes แบบกล่องสองใบ (JavaScript ล้วน ไม่มี library)
   ใช้ซ้ำได้ทุกบทในเส้น thesis · แค่วาง markup ไม่ต้องเขียน JS เพิ่ม

   ── เกมกล่องลับ ─────────────────────────────────────────
     <div class="bayes-game" data-ra="0.3" data-rb="0.7" data-pb="0.5"></div>
     ra / rb = สัดส่วนลูกแดงในกล่อง A / B · pb = สัดส่วนกล่อง B ในถุง (คอมพิวเตอร์สุ่มหยิบกล่องจากถุง 1 ใบ)

   ── ตาราง 100 รอบ ──────────────────────────────────────
     <div class="bayes-grid" data-ra="0.3" data-rb="0.7" data-pb="0.5"></div>
     มีแถบเลื่อนปรับทั้งสามค่า (ทีละ 0.1 ⇒ ทุกจำนวนรอบเป็นจำนวนเต็มพอดี ไม่ต้องปัด)

   ── ระฆัง Gaussian / Alice กับ Bob / ยืดแกน ── ดู comment ที่ bell() · duel() · stretch() ข้างล่าง

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

  /* ---------- ระฆัง Gaussian: prior × likelihood → posterior ----------
     <div class="bayes-bell" data-mu="0" data-sw="3" data-m="10" data-s="1" data-n="1" data-hide="n"></div>
     mu/sw = จุดกลาง/ความกว้างของ prior (Σ) · m/s = ค่าที่วัดได้/ความคลาดของเครื่อง (σ)
     n = จำนวนครั้งที่วัด (m คือค่าเฉลี่ย) · data-hide = ชื่อแถบที่ซ่อน คั่นด้วยช่องว่าง */
  const COL = { prior: "#3b6fb6", like: "#d62728", post: "#2ca02c", bob: "#c8692a" };
  const gauss = (x, c, w) => Math.exp(-(x - c) * (x - c) / (2 * w * w)) / (w * Math.sqrt(2 * Math.PI));
  // สูตร (3)(4): precision บวกกัน · ค่ากลางถ่วงด้วย precision
  const update = (mu, S, m, s, n) => {
    const P = 1 / (S * S), D = n / (s * s), T = P + D;
    return { P, D, T, c: (mu * P + m * D) / T, w: Math.sqrt(1 / T) };
  };
  function canvasFor(root, h) {
    const cv = el("canvas", "bb-canvas");
    const ctx = cv.getContext("2d");
    const size = () => {
      const W = Math.max(260, root.clientWidth - 32), r = window.devicePixelRatio || 1;
      cv.width = W * r; cv.height = h * r; cv.style.width = W + "px"; cv.style.height = h + "px";
      ctx.setTransform(r, 0, 0, r, 0, 0);
      return W;
    };
    return { cv, ctx, size };
  }
  function axes(ctx, W, H, x0, x1) {
    const ink = getComputedStyle(document.body).color || "#555";
    ctx.strokeStyle = ink; ctx.fillStyle = ink; ctx.globalAlpha = .6; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, H - 18); ctx.lineTo(W, H - 18); ctx.stroke();
    ctx.font = "12px sans-serif"; ctx.textAlign = "center";
    const step = (x1 - x0) > 30 ? 10 : 5;
    for (let t = Math.ceil(x0 / step) * step; t <= x1; t += step) {
      const X = (t - x0) / (x1 - x0) * W;
      ctx.fillText(String(t), X, H - 3);
      ctx.beginPath(); ctx.moveTo(X, H - 18); ctx.lineTo(X, H - 14); ctx.stroke();
    }
    ctx.textAlign = "right"; ctx.fillText("θ", W - 2, H - 22);
    ctx.globalAlpha = 1;
  }
  function curve(ctx, W, H, x0, x1, f, ymax, col, dash, fill) {
    ctx.beginPath();
    for (let i = 0; i <= 400; i++) {
      const x = x0 + (x1 - x0) * i / 400, X = W * i / 400, Y = (H - 18) - f(x) / ymax * (H - 30);
      i ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y);
    }
    ctx.strokeStyle = col; ctx.lineWidth = 3; ctx.setLineDash(dash ? [7, 5] : []); ctx.stroke(); ctx.setLineDash([]);
    if (fill) {
      ctx.lineTo(W, H - 18); ctx.lineTo(0, H - 18); ctx.closePath();
      ctx.globalAlpha = .15; ctx.fillStyle = col; ctx.fill(); ctx.globalAlpha = 1;
    }
  }

  function bell(root) {
    const st = { mu: num(root, "mu", 0), sw: num(root, "sw", 3), m: num(root, "m", 10), s: num(root, "s", 1),
                 n: num(root, "n", 1), boxes: false };
    const hide = (root.dataset.hide || "").split(/\s+/);
    const x0 = num(root, "xmin", -6), x1 = num(root, "xmax", 16);
    const ctrls = el("div", "bgr-ctrls");
    const mk = (key, label, lo, hi, step) => {
      const w = el("label", "bgr-slider");
      const s = el("input"); s.type = "range"; s.min = lo; s.max = hi; s.step = step; s.value = st[key];
      const out = el("output");
      s.oninput = () => { st[key] = parseFloat(s.value); draw(); };
      w.append(el("span", null, label), s, out);
      if (!hide.includes(key)) ctrls.append(w);
      return out;
    };
    const o = {
      mu: mk("mu", "prior: เชื่อว่าราว ๆ", -5, 15, 0.5), sw: mk("sw", "prior: คลาดได้ ±", 0.5, 6, 0.5),
      m: mk("m", "วัดได้", -5, 15, 0.5), s: mk("s", "เครื่องคลาดได้ ±", 0.5, 6, 0.5),
      n: mk("n", "วัดกี่ครั้ง", 1, 100, 1)
    };
    const tg = el("label", "bg-toggle"); const cb = el("input"); cb.type = "checkbox";
    tg.append(cb, document.createTextNode(" แสดงเป็นกล่อง (θ = จำนวนเต็ม แต่ละค่าคือกล่องหนึ่งใบแบบ T0)"));
    cb.onchange = () => { st.boxes = cb.checked; draw(); };
    const legend = el("p", "bb-legend");
    [["prior", COL.prior, false], ["likelihood", COL.like, true], ["posterior", COL.post, false]].forEach(([t, c, d]) => {
      const k = el("span", "bb-key"); k.style.borderTop = `3px ${d ? "dashed" : "solid"} ${c}`; legend.append(k, document.createTextNode(" " + t + "  "));
    });
    const { cv, ctx, size } = canvasFor(root, 220);
    const out = el("div", "bgr-out");
    root.append(ctrls, tg, legend, cv, out);

    function draw() {
      const W = size(), H = 220;
      ctx.clearRect(0, 0, W, H);
      for (const k in o) o[k].textContent = String(st[k]);
      const r = update(st.mu, st.sw, st.m, st.s, st.n);
      const sd = st.s / Math.sqrt(st.n);                       // ข้อมูล n ตัว = วัดครั้งเดียวที่แม่นขึ้น √n เท่า
      const fP = x => gauss(x, st.mu, st.sw), fL = x => gauss(x, st.m, sd), fQ = x => gauss(x, r.c, r.w);
      const ymax = Math.max(fP(st.mu), fL(st.m), fQ(r.c)) * 1.08;
      axes(ctx, W, H, x0, x1);
      if (st.boxes) {
        // T0 ทุกขั้น: กล่องละ prior × likelihood แล้วหารด้วยผลรวม ⇒ แท่งต้องทับเส้นเขียวพอดี
        const xs = []; for (let t = Math.ceil(x0); t <= x1; t++) xs.push(t);
        const prod = xs.map(t => fP(t) * fL(t)), Z = prod.reduce((a, b) => a + b, 0);
        const bw = W / (x1 - x0) * 0.6;
        ctx.fillStyle = COL.post; ctx.globalAlpha = .45;
        xs.forEach((t, i) => {
          const hh = Z ? prod[i] / Z / ymax * (H - 30) : 0, X = (t - x0) / (x1 - x0) * W;
          ctx.fillRect(X - bw / 2, (H - 18) - hh, bw, hh);
        });
        ctx.globalAlpha = 1;
      }
      curve(ctx, W, H, x0, x1, fP, ymax, COL.prior, false, false);
      curve(ctx, W, H, x0, x1, fL, ymax, COL.like, true, false);
      curve(ctx, W, H, x0, x1, fQ, ymax, COL.post, false, true);
      out.textContent = "";
      const line = (t, cls) => out.append(el("p", cls || null, t));
      const f = x => (Math.round(x * 100) / 100).toString();
      line("precision ของ prior = 1/" + f(st.sw * st.sw) + " = " + f(r.P));
      line("precision ของข้อมูล = " + st.n + "/" + f(st.s * st.s) + " = " + f(r.D));
      line("precision ของ posterior = " + f(r.P) + " + " + f(r.D) + " = " + f(r.T) + "  ⇒ คลาดได้ ±" + f(r.w));
      line("posterior อยู่ราว ๆ " + f(r.c) + " ± " + f(r.w), "bgr-answer");
    }
    window.addEventListener("resize", draw);
    draw();
  }

  /* ---------- Alice กับ Bob: วัดเพิ่มทีละครั้ง ----------
     <div class="bayes-duel" data-truth="10"></div>
     Alice prior 0 ± 3 · Bob prior 20 ± 1 · เครื่องคลาด ±1 · ค่าจริงซ่อนอยู่ */
  function duel(root) {
    const truth = num(root, "truth", 10), s = 1, A = [0, 3], B = [20, 1], x0 = -6, x1 = 26;
    let data = [];
    const randn = () => { // Box–Muller
      let u = 0; while (!u) u = Math.random();
      return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random());
    };
    const btns = el("div", "bg-btns");
    const b1 = el("button", "bg-draw", "📏 วัดเพิ่ม 1 ครั้ง"), b10 = el("button", null, "+10 ครั้ง"),
          b0 = el("button", "bg-new", "↻ เริ่มใหม่");
    [b1, b10, b0].forEach(b => b.type = "button");
    btns.append(b1, b10, b0);
    const legend = el("p", "bb-legend");
    [["Alice (เริ่ม 0 ± 3)", COL.prior], ["Bob (เริ่ม 20 ± 1)", COL.bob]].forEach(([t, c]) => {
      const k = el("span", "bb-key"); k.style.borderTop = `3px solid ${c}`; legend.append(k, document.createTextNode(" " + t + "  "));
    });
    const { cv, ctx, size } = canvasFor(root, 200);
    const out = el("div", "bgr-out");
    root.append(btns, legend, cv, out);
    function draw() {
      const W = size(), H = 200, n = data.length;
      const mean = n ? data.reduce((a, b) => a + b, 0) / n : 0;
      const a = n ? update(A[0], A[1], mean, s, n) : { c: A[0], w: A[1] };
      const b = n ? update(B[0], B[1], mean, s, n) : { c: B[0], w: B[1] };
      const fa = x => gauss(x, a.c, a.w), fb = x => gauss(x, b.c, b.w);
      const ymax = Math.max(fa(a.c), fb(b.c)) * 1.08;
      ctx.clearRect(0, 0, W, H);
      axes(ctx, W, H, x0, x1);
      curve(ctx, W, H, x0, x1, fa, ymax, COL.prior, false, true);
      curve(ctx, W, H, x0, x1, fb, ymax, COL.bob, false, true);
      out.textContent = "";
      const f = x => (Math.round(x * 100) / 100).toString();
      out.append(el("p", null, "วัดไปแล้ว " + n + " ครั้ง" + (n ? " · ค่าเฉลี่ยของข้อมูล " + f(mean) : "")));
      out.append(el("p", null, "Alice: " + f(a.c) + " ± " + f(a.w) + " · Bob: " + f(b.c) + " ± " + f(b.w)));
      out.append(el("p", "bgr-answer", "สองคนห่างกัน " + f(Math.abs(a.c - b.c))));
    }
    const add = k => { for (let i = 0; i < k; i++) data.push(truth + s * randn()); draw(); };
    b1.onclick = () => add(1); b10.onclick = () => add(10); b0.onclick = () => { data = []; draw(); };
    window.addEventListener("resize", draw);
    draw();
  }

  /* ---------- ยืดแกน: สุ่มด้านยาว L แล้วดูพื้นที่ A = L² ----------
     <div class="bayes-stretch" data-n="2000"></div>
     ปุ่มเลือกว่าจะสุ่ม "แบนบน L" หรือ "แบนบน A" · จุดเดียวกันโผล่บนสองแกน มีเส้นเชื่อม
     histogram 8 ช่องบนแต่ละแกน ⇒ เห็นว่าแบนบนแกนหนึ่ง ไม่แบนบนอีกแกน */
  function stretch(root) {
    const N = num(root, "n", 2000), L0 = 1, L1 = 3, A0 = 1, A1 = 9, H = 260;
    let mode = "L", pts = [];
    const btns = el("div", "bg-btns");
    const bL = el("button", null, "สุ่มแบบแบนบนด้านยาว L"), bA = el("button", null, "สุ่มแบบแบนบนพื้นที่ A"),
          bR = el("button", "bg-new", "↻ สุ่มใหม่");
    [bL, bA, bR].forEach(b => b.type = "button");
    btns.append(bL, bA, bR);
    const { cv, ctx, size } = canvasFor(root, H);
    const out = el("div", "bgr-out");
    root.append(btns, cv, out);
    function roll() {
      pts = [];
      for (let i = 0; i < N; i++) {
        const L = mode === "L" ? L0 + (L1 - L0) * Math.random() : Math.sqrt(A0 + (A1 - A0) * Math.random());
        pts.push(L);
      }
      draw();
    }
    function hist(vals, lo, hi, k) {
      const h = new Array(k).fill(0);
      vals.forEach(v => { h[Math.min(k - 1, Math.floor((v - lo) / (hi - lo) * k))]++; });
      return h;
    }
    function draw() {
      const W = size(), pad = 16, w = W - 2 * pad;
      const ink = getComputedStyle(document.body).color || "#555";
      ctx.clearRect(0, 0, W, H);
      const yL = 95, yA = 185, bh = 70;
      const xL = L => pad + (L - L0) / (L1 - L0) * w, xA = A => pad + (A - A0) / (A1 - A0) * w;
      // histogram 8 ช่อง: L ขึ้นด้านบนแกน · A ลงด้านล่างแกน
      const hL = hist(pts, L0, L1, 8), hA = hist(pts.map(L => L * L), A0, A1, 8);
      const mx = Math.max(...hL, ...hA, 1);
      ctx.globalAlpha = .35;
      ctx.fillStyle = COL.prior;
      hL.forEach((c, i) => { const h = c / mx * bh; ctx.fillRect(pad + i * w / 8 + 1, yL - h, w / 8 - 2, h); });
      ctx.fillStyle = COL.bob;
      hA.forEach((c, i) => { const h = c / mx * bh; ctx.fillRect(pad + i * w / 8 + 1, yA, w / 8 - 2, h); });
      // เส้นเชื่อม: จุดเดียวกันบนสองแกน (วาดแค่ 60 เส้นแรกพอให้เห็นการยืด)
      ctx.globalAlpha = .25; ctx.strokeStyle = ink; ctx.lineWidth = 1;
      pts.slice(0, 60).forEach(L => { ctx.beginPath(); ctx.moveTo(xL(L), yL); ctx.lineTo(xA(L * L), yA); ctx.stroke(); });
      ctx.globalAlpha = 1;
      // แกน + ขีดเลข
      ctx.strokeStyle = ink; ctx.fillStyle = ink; ctx.font = "12px sans-serif"; ctx.textAlign = "center";
      [[yL, L0, L1, xL, "ด้านยาว L", [1, 1.5, 2, 2.5, 3]], [yA, A0, A1, xA, "พื้นที่ A = L²", [1, 3, 5, 7, 9]]].forEach(([y, lo, hi, X, name, ticks]) => {
        ctx.beginPath(); ctx.moveTo(pad, y); ctx.lineTo(W - pad, y); ctx.stroke();
        ticks.forEach(t => { ctx.beginPath(); ctx.moveTo(X(t), y - 3); ctx.lineTo(X(t), y + 3); ctx.stroke(); ctx.fillText(String(t), X(t), y === yL ? y + 15 : y - 6); });
        // ป้ายชื่อแกน: L อยู่มุมซ้ายบน · A อยู่ชิดขวาเหนือแกน (มุมซ้ายล่างทับแท่ง)
        if (y === yL) { ctx.textAlign = "left"; ctx.fillText(name, pad, 14); }
        else { ctx.textAlign = "right"; ctx.fillText(name, W - pad, y - 20); }
        ctx.textAlign = "center";
      });
      // เส้นแบ่ง L = 2 ⇔ A = 4
      ctx.strokeStyle = COL.like; ctx.setLineDash([5, 4]); ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(xL(2), yL - bh); ctx.lineTo(xL(2), yL); ctx.lineTo(xA(4), yA); ctx.lineTo(xA(4), yA + bh); ctx.stroke();
      ctx.setLineDash([]);
      const k = pts.filter(L => L < 2).length;
      out.textContent = "";
      out.append(el("p", null, "กำลังสุ่มแบบ: แบนบน" + (mode === "L" ? "ด้านยาว L" : "พื้นที่ A") + " · " + N + " จุด"));
      out.append(el("p", "bgr-answer", "ด้านยาวน้อยกว่า 2 (เส้นประแดง): " + k + " จาก " + N + " = " + Math.round(k / N * 100) + "%"));
      bL.classList.toggle("bg-draw", mode === "L"); bA.classList.toggle("bg-draw", mode === "A");
    }
    bL.onclick = () => { mode = "L"; roll(); };
    bA.onclick = () => { mode = "A"; roll(); };
    bR.onclick = roll;
    window.addEventListener("resize", draw);
    roll();
  }

  function init() {
    document.querySelectorAll(".bayes-game").forEach(game);
    document.querySelectorAll(".bayes-grid").forEach(grid);
    document.querySelectorAll(".bayes-bell").forEach(bell);
    document.querySelectorAll(".bayes-duel").forEach(duel);
    document.querySelectorAll(".bayes-stretch").forEach(stretch);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
