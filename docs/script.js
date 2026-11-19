/* Vasicek Bond Risk Simulator: simulation and page behaviour. */
"use strict";

/* ---------- Reference values from the notebook (update these if the notebook output changes) ---------- */
const PRESETS = {
    calm: { kappa: 0.7655150217941049, theta: 0.08615227029073727, sigma: 0.014098880214804786, r0: 0.10082215839866182 },
    crisis: { kappa: 0.21423885415882352, theta: 0.1393768140848096, sigma: 0.05094228272509795, r0: 0.10082215839866182 }
};
const REF = { zeroPrice: 63.84, tolerance: 0.05, minPaths: 50000 };   // notebook: 10 replications x 100,000 paths, weekly step

/* ---------- Fixed settings (same as the notebook) ---------- */
const T_BOND = 5;            // years to maturity
const FACE = 100;            // LKR
const STEPS_PER_YEAR = 52;   // weekly Euler-Maruyama step
const COUPONS_PER_YEAR = 2;
const HORIZON = 1;           // years, for Value-at-Risk
const ALPHA = 0.95;
const Z_ALPHA = 1.6448536269514722;   // 95th percentile of the standard normal
const BATCH = 250;           // paths added per animation frame
const STORED_PATHS = 400;    // full paths kept for the chart and its percentile bands
const DRAWN_PATHS = 60;      // individual paths drawn
const REVEAL_FRAMES = 80;    // frames taken to sweep the paths across the chart

/* ---------- Random numbers: seeded uniform generator and Box-Muller normal variates ---------- */
function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
function makeNormal(uniform) {
    let spare = null;
    return function () {
        if (spare !== null) { const s = spare; spare = null; return s; }
        let u = 0;
        while (u === 0) u = uniform();
        const v = uniform();
        const m = Math.sqrt(-2 * Math.log(u));
        spare = m * Math.sin(2 * Math.PI * v);
        return m * Math.cos(2 * Math.PI * v);
    };
}

/* ---------- Exact Vasicek results ---------- */
// Zero-coupon price P = A(tau) * exp(-B(tau) * r)
function vasicekAB(tau, p) {
    const B = (1 - Math.exp(-p.kappa * tau)) / p.kappa;
    const lnA = (p.theta - p.sigma * p.sigma / (2 * p.kappa * p.kappa)) * (B - tau) - p.sigma * p.sigma * B * B / (4 * p.kappa);
    return { A: Math.exp(lnA), B: B };
}
function vasicekZero(r, tau, p) {
    const ab = vasicekAB(tau, p);
    return ab.A * Math.exp(-ab.B * r);
}
// Mean and standard deviation of the rate after t years
function rateMean(t, p) { return p.theta + (p.r0 - p.theta) * Math.exp(-p.kappa * t); }
function rateSd(t, p) { return p.sigma * Math.sqrt((1 - Math.exp(-2 * p.kappa * t)) / (2 * p.kappa)); }

/* ---------- Bond cash flows ---------- */
function cashflows(couponRate) {
    const n = T_BOND * COUPONS_PER_YEAR, times = [], amounts = [];
    for (let j = 1; j <= n; j++) {
        times.push(j / COUPONS_PER_YEAR);
        amounts.push(FACE * couponRate / COUPONS_PER_YEAR + (j === n ? FACE : 0));
    }
    return { times: times, amounts: amounts };
}
function yieldToMaturity(price, cf) {
    let lo = -0.5, hi = 2;                       // bisection on the semi-annual yield
    for (let i = 0; i < 80; i++) {
        const y = (lo + hi) / 2;
        let pv = 0;
        for (let j = 0; j < cf.times.length; j++) pv += cf.amounts[j] / Math.pow(1 + y / 2, cf.times[j] * 2);
        if (pv > price) lo = y; else hi = y;
    }
    return (lo + hi) / 2;
}
// Value of the position after one year if the rate then is r1 (remaining payments + coupons already received)
function valueAtHorizon(r1, cf, p) {
    let v = 0;
    for (let j = 0; j < cf.times.length; j++) {
        v += cf.times[j] <= HORIZON + 1e-9 ? cf.amounts[j] : cf.amounts[j] * vasicekZero(r1, cf.times[j] - HORIZON, p);
    }
    return v;
}
// Formula values shown before any simulation has run
function formulaPreview(p, cf) {
    const df = cf.times.map(function (t) { return vasicekZero(p.r0, t, p); });
    let price = 0;
    for (let j = 0; j < df.length; j++) price += cf.amounts[j] * df[j];
    const badRate = rateMean(HORIZON, p) + Z_ALPHA * rateSd(HORIZON, p);     // 95th percentile of the one-year rate
    return { price: price, df: df, valueAtRisk: price - valueAtHorizon(badRate, cf, p) };
}

/* ---------- Simulation ---------- */
function newRun(p, couponRate, nTarget, seed) {
    const cf = cashflows(couponRate);
    return {
        p: p, cf: cf, nTarget: nTarget, done: 0,
        normal: makeNormal(mulberry32(seed)),
        sumZ: 0, sumZ2: 0, sumC: 0, sumC2: 0,
        sumD: new Float64Array(cf.times.length),   // running sum of discount factors at each payment date
        r1: new Float64Array(nTarget),             // rate at the one-year horizon, per path
        stored: [], fan: null, conv: [], risk: null, reveal: 0
    };
}
// Simulate `count` more paths with the Euler-Maruyama scheme.
function simulateBatch(run, count) {
    const p = run.p, cf = run.cf;
    const dt = 1 / STEPS_PER_YEAR, sqrtDt = Math.sqrt(dt);
    const nSteps = T_BOND * STEPS_PER_YEAR;
    const stepsPerCoupon = STEPS_PER_YEAR / COUPONS_PER_YEAR;
    const horizonStep = HORIZON * STEPS_PER_YEAR;
    for (let i = 0; i < count && run.done < run.nTarget; i++) {
        const keep = run.stored.length < STORED_PATHS ? new Float32Array(nSteps + 1) : null;
        let r = p.r0, integral = 0, pvCoupon = 0, discount = 1;
        if (keep) keep[0] = r;
        for (let n = 1; n <= nSteps; n++) {
            const rNew = r + p.kappa * (p.theta - r) * dt + p.sigma * sqrtDt * run.normal();   // Euler-Maruyama step
            integral += 0.5 * (r + rNew) * dt;                                                 // trapezoidal rule
            r = rNew;
            if (keep) keep[n] = r;
            if (n === horizonStep) run.r1[run.done] = r;
            if (n % stepsPerCoupon === 0) {
                const j = n / stepsPerCoupon - 1;
                discount = Math.exp(-integral);
                run.sumD[j] += discount;
                pvCoupon += cf.amounts[j] * discount;
            }
        }
        const pvZero = FACE * discount;
        run.sumZ += pvZero; run.sumZ2 += pvZero * pvZero;
        run.sumC += pvCoupon; run.sumC2 += pvCoupon * pvCoupon;
        if (keep) run.stored.push(keep);
        run.done++;
    }
    const z = estimate(run.sumZ, run.sumZ2, run.done);
    run.conv.push([run.done, z.mean, z.half]);
    if (!run.fan && (run.stored.length >= STORED_PATHS || run.done >= run.nTarget)) run.fan = percentileFan(run.stored);
}
// Mean and 95% confidence half-width from running sums.
function estimate(sum, sum2, n) {
    if (n < 2) return { mean: NaN, half: NaN };
    const mean = sum / n;
    const variance = Math.max(0, (sum2 - n * mean * mean) / (n - 1));
    return { mean: mean, half: 1.96 * Math.sqrt(variance / n) };
}
// Percentiles of the stored paths at every fourth step, for the shaded bands.
function percentileFan(paths) {
    const fan = { step: [], q05: [], q25: [], q50: [], q75: [], q95: [] };
    const col = new Float64Array(paths.length);
    for (let n = 0; n <= T_BOND * STEPS_PER_YEAR; n += 4) {
        for (let i = 0; i < paths.length; i++) col[i] = paths[i][n];
        const s = Float64Array.from(col).sort();
        const q = function (f) { return s[Math.round(f * (s.length - 1))]; };
        fan.step.push(n); fan.q05.push(q(0.05)); fan.q25.push(q(0.25)); fan.q50.push(q(0.5)); fan.q75.push(q(0.75)); fan.q95.push(q(0.95));
    }
    return fan;
}
// One-year loss distribution: revalue the position at each path's one-year rate.
function computeRisk(run) {
    const n = run.done, cf = run.cf, p = run.p;
    if (n < 200) return null;
    const v0 = run.sumC / n;
    let received = 0; const A = [], B = [], c = [];
    for (let j = 0; j < cf.times.length; j++) {
        if (cf.times[j] <= HORIZON + 1e-9) received += cf.amounts[j];
        else { const ab = vasicekAB(cf.times[j] - HORIZON, p); A.push(ab.A); B.push(ab.B); c.push(cf.amounts[j]); }
    }
    const loss = new Float64Array(n);
    let nLoss = 0;
    for (let i = 0; i < n; i++) {
        let v1 = received;
        for (let j = 0; j < c.length; j++) v1 += c[j] * A[j] * Math.exp(-B[j] * run.r1[i]);
        loss[i] = v0 - v1;
        if (loss[i] > 0) nLoss++;
    }
    const sorted = loss.sort();
    const pos = ALPHA * (n - 1), k = Math.floor(pos);
    const valueAtRisk = sorted[k] + (pos - k) * (sorted[Math.min(k + 1, n - 1)] - sorted[k]);
    let tail = 0, nTail = 0;
    for (let i = n - 1; i >= 0 && sorted[i] >= valueAtRisk; i--) { tail += sorted[i]; nTail++; }
    return { valueAtRisk: valueAtRisk, shortfall: tail / nTail, pLoss: nLoss / n, sorted: sorted };
}

/* ---------- Page wiring ---------- */
(function () {
    const $ = function (id) { return document.getElementById(id); };
    const pct = function (v) { return (v * 100).toFixed(2) + "%"; };
    const sliders = {
        kappa: { el: $("kappa"), out: $("kappaOut"), scale: 1, fmt: function (v) { return v.toFixed(3) + " / yr"; } },
        theta: { el: $("theta"), out: $("thetaOut"), scale: 100, fmt: pct },
        sigma: { el: $("sigma"), out: $("sigmaOut"), scale: 100, fmt: pct },
        r0: { el: $("r0"), out: $("r0Out"), scale: 100, fmt: pct },
        coupon: { el: $("coupon"), out: $("couponOut"), scale: 100, fmt: pct }
    };
    const values = { kappa: 0, theta: 0, sigma: 0, r0: 0, coupon: 0.10 };   // model units (rates as decimals)
    const reducedMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    let activePreset = null, run = null, running = false, frame = 0, rafId = 0;
    let colors = {}, hoverT = null, drag = null, frozenRange = null;
    const hit = { paths: null, conv: null, loss: null, flow: null };

    function modelParams() { return { kappa: values.kappa, theta: values.theta, sigma: values.sigma, r0: values.r0 }; }
    function showValues() {
        for (const key in sliders) {
            const s = sliders[key];
            s.el.value = values[key] * s.scale;
            s.out.textContent = s.fmt(values[key]);
            s.el.style.setProperty("--fill", ((s.el.value - s.el.min) / (s.el.max - s.el.min) * 100) + "%");
        }
        $("presetCalm").setAttribute("aria-pressed", String(activePreset === "calm"));
        $("presetCrisis").setAttribute("aria-pressed", String(activePreset === "crisis"));
        $("presetNote").textContent = activePreset ? "Values fitted to Central Bank of Sri Lanka data." : "Custom values.";
    }
    function loadPreset(name) {
        Object.assign(values, PRESETS[name]);
        activePreset = name;
        showValues(); resetRun();
    }
    function setValue(key, v) {
        values[key] = v;
        if (key !== "coupon") activePreset = null;
        showValues(); resetRun();
    }
    for (const key in sliders) {
        sliders[key].el.addEventListener("input", function () { setValue(key, parseFloat(sliders[key].el.value) / sliders[key].scale); });
    }
    $("paths").addEventListener("change", resetRun);
    $("seed").addEventListener("change", resetRun);
    $("presetCalm").addEventListener("click", function () { loadPreset("calm"); });
    $("presetCrisis").addEventListener("click", function () { loadPreset("crisis"); });

    function resetRun() {
        cancelAnimationFrame(rafId);
        running = false; run = null; frame = 0;
        $("start").disabled = false; $("pause").disabled = true;
        render();
    }
    function start() {
        if (!run) {
            const seed = Math.max(0, Math.floor(Number($("seed").value) || 0));
            run = newRun(modelParams(), values.coupon, parseInt($("paths").value, 10), seed);
            run.range = rateRange(run.p);
            if (reducedMotion) run.reveal = 1;
        }
        if (run.done >= run.nTarget) return;
        running = true;
        $("start").disabled = true; $("pause").disabled = false;
        rafId = requestAnimationFrame(tick);
    }
    function pause() {
        running = false;
        cancelAnimationFrame(rafId);
        $("start").disabled = false; $("pause").disabled = true;
        if (run) { run.risk = computeRisk(run); render(); }
    }
    function tick() {
        if (!running) return;
        simulateBatch(run, BATCH);
        frame++;
        run.reveal = Math.min(1, run.reveal + 1 / REVEAL_FRAMES);
        const finished = run.done >= run.nTarget;
        if (finished) run.reveal = 1;
        if (finished || frame % 8 === 0) run.risk = computeRisk(run);
        render();
        if (finished) { running = false; $("start").disabled = true; $("pause").disabled = true; }
        else rafId = requestAnimationFrame(tick);
    }
    $("start").addEventListener("click", start);
    $("pause").addEventListener("click", pause);
    $("reset").addEventListener("click", resetRun);

    /* ----- Numbers and sentences ----- */
    const f2 = function (v) { return v.toFixed(2).replace("-", "\u2212"); };
    function render() {
        const cs = getComputedStyle(document.documentElement);
        ["--ink", "--muted", "--line", "--path", "--path-soft", "--path-mid", "--mean", "--loss", "--loss-soft", "--surface", "--bg"].forEach(function (k) { colors[k] = cs.getPropertyValue(k).trim(); });

        const p = modelParams(), cf = cashflows(values.coupon), preview = formulaPreview(p, cf);
        const exactZero = FACE * vasicekZero(p.r0, T_BOND, p);
        const n = run ? run.done : 0, target = run ? run.nTarget : parseInt($("paths").value, 10);
        $("oExact").textContent = f2(exactZero);
        $("oN").textContent = n.toLocaleString("en-US");
        $("meterText").textContent = n.toLocaleString("en-US") + " of " + target.toLocaleString("en-US") + " paths";
        $("meterBar").style.width = (n / target * 100) + "%";
        const verdict = $("verdict"), status = $("status"), varEl = $("oVar"), esEl = $("oEs");

        if (n < 2) {
            $("valueCap").textContent = "Value of the bond today by exact formula, LKR per 100 face value";
            $("bigPrice").textContent = preview.price.toFixed(2); $("bigCi").textContent = "";
            $("oCoupon").textContent = "–"; $("oZero").textContent = "–"; $("oEs").textContent = "–"; $("oPloss").textContent = "–";
            $("oYtm").textContent = pct(yieldToMaturity(preview.price, cf));
            varEl.innerHTML = f2(preview.valueAtRisk) + " <small>formula</small>";
            varEl.className = preview.valueAtRisk < 0 ? "gain" : "lossy"; esEl.className = "";
            verdict.className = "verdict idle";
            verdict.textContent = "Press Start to test these formula values with a simulation and to see the full range of outcomes.";
            $("pathsCap").textContent = "Before simulation: the expected path of the rate and the range that should hold 90% of futures. Dashed line: long-run rate θ.";
            status.className = "status";
            status.textContent = "Agreement with the notebook is checked after a run with the Calm period values and at least " + REF.minPaths.toLocaleString("en-US") + " paths.";
        } else {
            const z = estimate(run.sumZ, run.sumZ2, n), c = estimate(run.sumC, run.sumC2, n);
            $("valueCap").textContent = "Value of the bond today by simulation, LKR per 100 face value";
            $("bigPrice").textContent = c.mean.toFixed(2); $("bigCi").textContent = "± " + c.half.toFixed(2);
            $("oCoupon").textContent = f2(c.mean) + " ± " + c.half.toFixed(2);
            $("oZero").innerHTML = f2(z.mean) + " <small>± " + z.half.toFixed(2) + "</small>";
            $("oYtm").textContent = pct(yieldToMaturity(c.mean, run.cf));
            $("pathsCap").textContent = "Each line is one simulated future for the rate. Shaded bands hold 50% and 90% of paths; the bold line is the median.";
            if (run.risk) {
                const k = run.risk, gain = k.valueAtRisk < 0;
                varEl.textContent = f2(k.valueAtRisk); varEl.className = gain ? "gain" : "lossy";
                esEl.textContent = f2(k.shortfall); esEl.className = k.shortfall < 0 ? "gain" : "lossy";
                $("oPloss").textContent = (k.pLoss * 100).toFixed(1) + "%";
                verdict.className = "verdict" + (gain ? "" : " loss");
                verdict.textContent = gain
                    ? "Even in the worst 5% of simulated years, the holder ends at least " + Math.abs(k.valueAtRisk).toFixed(2) + " ahead per 100 of face value: the coupons outweigh the fall in price."
                    : "In the worst 5% of simulated years, the holder loses at least " + k.valueAtRisk.toFixed(2) + " per 100 of face value, and " + k.shortfall.toFixed(2) + " on average.";
            }
            if (activePreset === "calm" && n >= REF.minPaths) {
                const diff = z.mean - REF.zeroPrice, ok = Math.abs(diff) <= REF.tolerance;
                status.className = "status " + (ok ? "ok" : "bad");
                status.textContent = (ok ? "Matches the notebook. " : "Does not match the notebook. ") +
                    "Zero-coupon price here " + z.mean.toFixed(3) + ", notebook " + REF.zeroPrice.toFixed(2) +
                    ", difference " + f2(diff) + " (tolerance " + REF.tolerance + " LKR).";
            } else {
                status.className = "status";
                status.textContent = activePreset === "calm"
                    ? "Agreement with the notebook is checked once " + REF.minPaths.toLocaleString("en-US") + " paths have run."
                    : "To check agreement with the notebook, choose the Calm period values and run at least " + REF.minPaths.toLocaleString("en-US") + " paths.";
            }
        }
        drawPaths(p); drawConvergence(exactZero); drawLoss(); drawFlow(p, cf, preview);
    }

    /* ----- Canvas helpers ----- */
    function setup(canvas, m) {
        const ratio = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;
        if (canvas.width !== Math.round(w * ratio) || canvas.height !== Math.round(h * ratio)) { canvas.width = Math.round(w * ratio); canvas.height = Math.round(h * ratio); }
        const g = canvas.getContext("2d");
        g.setTransform(ratio, 0, 0, ratio, 0, 0);
        g.clearRect(0, 0, w, h);
        g.font = "12px system-ui, sans-serif"; g.globalAlpha = 1;
        return Object.assign({ g: g, w: w, h: h }, m);
    }
    function niceTicks(lo, hi, count) {
        const raw = (hi - lo) / count, mag = Math.pow(10, Math.floor(Math.log10(raw)));
        const step = [1, 2, 2.5, 5, 10].map(function (k) { return k * mag; }).find(function (s) { return s >= raw; });
        const ticks = [];
        for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) ticks.push(Math.abs(v) < 1e-12 ? 0 : v);
        return ticks;
    }
    function axes(c, x0, x1, y0, y1, xLabel, yFmt, xFmt, xTicks) {
        const g = c.g;
        const X = function (x) { return c.l + (x - x0) / (x1 - x0) * (c.w - c.l - c.r); };
        const Y = function (y) { return c.h - c.b - (y - y0) / (y1 - y0) * (c.h - c.t - c.b); };
        g.strokeStyle = colors["--line"]; g.fillStyle = colors["--muted"]; g.lineWidth = 1;
        g.textAlign = "right"; g.textBaseline = "middle";
        niceTicks(y0, y1, 5).forEach(function (v) {
            g.beginPath(); g.moveTo(c.l, Math.round(Y(v)) + 0.5); g.lineTo(c.w - c.r, Math.round(Y(v)) + 0.5); g.stroke();
            g.fillText(yFmt(v), c.l - 7, Y(v));
        });
        g.textAlign = "center"; g.textBaseline = "top";
        (xTicks || niceTicks(x0, x1, 5)).forEach(function (v) { g.fillText(xFmt(v), X(v), c.h - c.b + 6); });
        g.fillText(xLabel, c.l + (c.w - c.l - c.r) / 2, c.h - 15);
        return {
            X: X, Y: Y, x0: x0, x1: x1, y0: y0, y1: y1,
            invX: function (px) { return x0 + (px - c.l) / (c.w - c.l - c.r) * (x1 - x0); },
            invY: function (py) { return y0 + (c.h - c.b - py) / (c.h - c.t - c.b) * (y1 - y0); }
        };
    }
    function line(g, pts, color, width, dash) {
        g.save(); g.strokeStyle = color; g.lineWidth = width; if (dash) g.setLineDash(dash);
        g.beginPath(); pts.forEach(function (q, i) { if (i === 0) g.moveTo(q[0], q[1]); else g.lineTo(q[0], q[1]); }); g.stroke(); g.restore();
    }
    function band(g, upper, lower, color) {
        g.fillStyle = color; g.beginPath();
        upper.forEach(function (q, i) { if (i === 0) g.moveTo(q[0], q[1]); else g.lineTo(q[0], q[1]); });
        for (let i = lower.length - 1; i >= 0; i--) g.lineTo(lower[i][0], lower[i][1]);
        g.closePath(); g.fill();
    }
    function emptyNote(c, text) { c.g.fillStyle = colors["--muted"]; c.g.textAlign = "center"; c.g.textBaseline = "middle"; c.g.fillText(text, c.l + (c.w - c.l - c.r) / 2, c.t + (c.h - c.t - c.b) / 2); }
    function tag(g, text, x, y, color, align) {
        g.font = "600 12px system-ui, sans-serif";
        const w = g.measureText(text).width + 12, left = align === "right" ? x - w : x;
        g.fillStyle = color; g.beginPath();
        if (g.roundRect) g.roundRect(left, y - 10, w, 20, 5); else g.rect(left, y - 10, w, 20);
        g.fill();
        g.fillStyle = colors["--surface"]; g.textAlign = "left"; g.textBaseline = "middle"; g.fillText(text, left + 6, y + 0.5);
        g.font = "12px system-ui, sans-serif";
    }

    /* ----- Chart 1: rate paths ----- */
    function rateRange(p) {
        const sd = p.sigma / Math.sqrt(2 * p.kappa);                       // long-run spread of the rate
        let lo = Math.min(p.r0, p.theta) - 2.7 * sd, hi = Math.max(p.r0, p.theta) + 2.7 * sd;
        const pad = Math.max(0.006, (hi - lo) * 0.1);
        return [lo - pad, hi + pad];
    }
    function drawPaths(p) {
        const c = setup($("cPaths"), { l: 50, r: 16, t: 16, b: 44 }), g = c.g;
        const range = frozenRange || (run ? run.range : rateRange(p));
        const a = axes(c, 0, T_BOND, range[0] * 100, range[1] * 100, "years from today", function (v) { return v.toFixed(1) + "%"; }, String, [0, 1, 2, 3, 4, 5]);
        hit.paths = { a: a, c: c };
        const X = a.X, Y = function (r) { return a.Y(r * 100); };

        // one-year marker: where the risk is measured
        line(g, [[X(HORIZON), c.t], [X(HORIZON), c.h - c.b]], colors["--muted"], 1, [2, 4]);
        g.fillStyle = colors["--muted"]; g.textAlign = "left"; g.textBaseline = "top";
        g.fillText("1 year: risk is measured here", X(HORIZON) + 6, c.t);

        g.save(); g.beginPath(); g.rect(c.l, c.t, c.w - c.l - c.r, c.h - c.t - c.b); g.clip();
        if (!run || run.stored.length === 0) {
            const up = [], dn = [], mid = [];
            for (let i = 0; i <= 100; i++) {
                const t = T_BOND * i / 100, m = rateMean(t, p), s = Z_ALPHA * rateSd(t, p);
                up.push([X(t), Y(m + s)]); dn.push([X(t), Y(m - s)]); mid.push([X(t), Y(m)]);
            }
            band(g, up, dn, colors["--path-soft"]);
            line(g, mid, colors["--path"], 2.5);
        } else {
            const tEdge = run.reveal * T_BOND;
            g.save(); g.beginPath(); g.rect(c.l, c.t, X(tEdge) - c.l, c.h); g.clip();
            if (run.fan) {
                const f = run.fan, pt = function (arr) { return arr.map(function (v, i) { return [X(f.step[i] / STEPS_PER_YEAR), Y(v)]; }); };
                band(g, pt(f.q95), pt(f.q05), colors["--path-soft"]);
                band(g, pt(f.q75), pt(f.q25), colors["--path-mid"]);
            }
            g.strokeStyle = colors["--path"]; g.lineWidth = 1; g.globalAlpha = 0.3;
            for (let i = 0; i < Math.min(DRAWN_PATHS, run.stored.length); i++) {
                const path = run.stored[i];
                g.beginPath();
                for (let n = 0; n < path.length; n++) { const x = X(n / STEPS_PER_YEAR), y = Y(path[n]); if (n === 0) g.moveTo(x, y); else g.lineTo(x, y); }
                g.stroke();
            }
            g.globalAlpha = 1;
            if (run.fan) line(g, run.fan.q50.map(function (v, i) { return [X(run.fan.step[i] / STEPS_PER_YEAR), Y(v)]; }), colors["--path"], 2.5);
            g.restore();
            if (run.reveal < 1) {                                             // glowing tips at the leading edge
                const n = Math.round(tEdge * STEPS_PER_YEAR);
                g.fillStyle = colors["--path"];
                for (let i = 0; i < Math.min(DRAWN_PATHS, run.stored.length); i++) { g.beginPath(); g.arc(X(tEdge), Y(run.stored[i][n]), 2.2, 0, 2 * Math.PI); g.fill(); }
            }
        }
        // long-run rate: draggable dashed line
        line(g, [[c.l, Y(p.theta)], [c.w - c.r, Y(p.theta)]], colors["--mean"], 2, [7, 5]);
        // hover crosshair
        if (hoverT !== null && !drag) line(g, [[X(hoverT), c.t], [X(hoverT), c.h - c.b]], colors["--ink"], 1);
        g.restore();

        tag(g, "θ " + pct(p.theta), c.w - c.r - 4, Math.min(c.h - c.b - 12, Math.max(c.t + 12, Y(p.theta) - 14)), colors["--mean"], "right");
        // starting rate: draggable dot
        const y0 = Math.min(c.h - c.b, Math.max(c.t, Y(p.r0)));
        g.fillStyle = colors["--surface"]; g.strokeStyle = colors["--ink"]; g.lineWidth = 3;
        g.beginPath(); g.arc(X(0), y0, 7, 0, 2 * Math.PI); g.fill(); g.stroke();
        tag(g, "r₀ " + pct(p.r0), X(0) + 14, y0 + (p.r0 >= p.theta ? -16 : 16), colors["--ink"], "left");
    }

    /* ----- Chart 2: running price estimate ----- */
    function drawConvergence(exact) {
        const c = setup($("cConv"), { l: 50, r: 10, t: 10, b: 44 }), g = c.g;
        hit.conv = null;
        const pts = run ? run.conv.filter(function (q) { return q[0] >= 500; }) : [];
        if (pts.length < 2) {
            const a0 = axes(c, 0, parseInt($("paths").value, 10), exact - 0.5, exact + 0.5, "paths simulated", function (v) { return v.toFixed(2); }, kfmt);
            line(g, [[a0.X(a0.x0), a0.Y(exact)], [a0.X(a0.x1), a0.Y(exact)]], colors["--mean"], 1.8, [7, 5]);
            return;
        }
        let lo = exact, hi = exact;
        pts.forEach(function (q) { lo = Math.min(lo, q[1] - q[2]); hi = Math.max(hi, q[1] + q[2]); });
        const pad = (hi - lo) * 0.1 || 0.01;
        const a = axes(c, 0, run.nTarget, lo - pad, hi + pad, "paths simulated", function (v) { return v.toFixed(2); }, kfmt);
        hit.conv = { a: a, pts: pts, exact: exact };
        band(g, pts.map(function (q) { return [a.X(q[0]), a.Y(q[1] + q[2])]; }), pts.map(function (q) { return [a.X(q[0]), a.Y(q[1] - q[2])]; }), colors["--path-mid"]);
        line(g, [[a.X(0), a.Y(exact)], [a.X(run.nTarget), a.Y(exact)]], colors["--mean"], 1.8, [7, 5]);
        line(g, pts.map(function (q) { return [a.X(q[0]), a.Y(q[1])]; }), colors["--path"], 2);
        const last = pts[pts.length - 1];
        g.fillStyle = colors["--path"]; g.beginPath(); g.arc(a.X(last[0]), a.Y(last[1]), 4, 0, 2 * Math.PI); g.fill();
    }
    function kfmt(v) { return v >= 1000 ? (v / 1000) + "k" : String(v); }

    /* ----- Chart 3: one-year gain or loss ----- */
    function drawLoss() {
        const c = setup($("cLoss"), { l: 44, r: 10, t: 24, b: 44 }), g = c.g;
        hit.loss = null;
        if (!run || !run.risk) { emptyNote(c, "Run the simulation to see the outcomes"); return; }
        const s = run.risk.sorted, n = s.length, v = run.risk.valueAtRisk;
        let lo = s[Math.floor(0.001 * (n - 1))], hi = s[Math.floor(0.999 * (n - 1))];
        if (hi - lo < 1e-6) { lo -= 0.5; hi += 0.5; }
        const bins = 40, counts = new Array(bins).fill(0), width = (hi - lo) / bins;
        for (let i = 0; i < n; i++) { const b = Math.floor((s[i] - lo) / width); if (b >= 0 && b < bins) counts[b]++; }
        const top = Math.max.apply(null, counts);
        const a = axes(c, lo, hi, 0, top * 1.05, "loss after one year (negative = gain)", function (y) { return (y / n * 100).toFixed(0) + "%"; }, function (x) { return f2(x).replace(/\.00$/, ""); });
        hit.loss = { a: a, lo: lo, width: width, counts: counts, n: n, bins: bins };
        if (v > lo && v < hi) { g.fillStyle = colors["--loss-soft"]; g.fillRect(a.X(v), c.t, a.X(hi) - a.X(v), c.h - c.t - c.b); }
        for (let b = 0; b < bins; b++) {
            const x0 = a.X(lo + b * width), x1 = a.X(lo + (b + 1) * width);
            g.fillStyle = (lo + (b + 0.5) * width) > 0 ? colors["--loss"] : colors["--path"];
            g.globalAlpha = 0.85;
            g.fillRect(x0, a.Y(counts[b]), Math.max(1, x1 - x0 - 1), a.Y(0) - a.Y(counts[b]));
        }
        g.globalAlpha = 1;
        if (lo < 0 && hi > 0) line(g, [[a.X(0), c.t], [a.X(0), c.h - c.b]], colors["--muted"], 1, [2, 3]);
        if (v > lo && v < hi) {
            line(g, [[a.X(v), c.t - 6], [a.X(v), c.h - c.b]], colors["--ink"], 2);
            tag(g, "95% VaR " + f2(v), a.X(v) - 4, c.t - 12 + 10, colors["--ink"], "right");
        }
    }

    /* ----- Chart 4: discount factors at each payment date ----- */
    function drawFlow(p, cf, preview) {
        const c = setup($("cFlow"), { l: 40, r: 10, t: 10, b: 44 }), g = c.g;
        const simulated = run && run.done >= 2;
        const df = simulated ? Array.from(run.sumD, function (s) { return s / run.done; }) : preview.df;
        const a = axes(c, 0, T_BOND + 0.25, 0, 1, "years until the payment", function (v) { return v.toFixed(1); }, String, [1, 2, 3, 4, 5]);
        hit.flow = { a: a, df: df, cf: cf };
        const w = (a.X(0.5) - a.X(0)) * 0.62;
        cf.times.forEach(function (t, j) {
            g.fillStyle = colors["--path"]; g.globalAlpha = simulated ? 0.85 : 0.35;
            g.fillRect(a.X(t) - w / 2, a.Y(df[j]), w, a.Y(0) - a.Y(df[j]));
            g.globalAlpha = 1;
            if (simulated) line(g, [[a.X(t) - w / 2 - 3, a.Y(preview.df[j])], [a.X(t) + w / 2 + 3, a.Y(preview.df[j])]], colors["--mean"], 2.5);
        });
    }

    /* ----- Hover read-outs ----- */
    function showTip(id, canvas, e, html) {
        const tip = $(id);
        if (!html) { tip.style.display = "none"; return; }
        tip.innerHTML = html; tip.style.display = "block";
        const box = canvas.getBoundingClientRect(), x = e.clientX - box.left, y = e.clientY - box.top;
        tip.style.left = Math.max(0, Math.min(x + 14, box.width - tip.offsetWidth)) + "px";
        tip.style.top = (canvas.offsetTop + Math.max(0, y - tip.offsetHeight - 10)) + "px";
    }
    function local(canvas, e) { const b = canvas.getBoundingClientRect(); return { x: e.clientX - b.left, y: e.clientY - b.top }; }

    const cPaths = $("cPaths");
    function nearHandle(pos) {
        const h = hit.paths; if (!h) return null;
        const p = modelParams(), yR0 = h.a.Y(p.r0 * 100), yTheta = h.a.Y(p.theta * 100);
        if (Math.hypot(pos.x - h.a.X(0), pos.y - yR0) < 16) return "r0";
        if (pos.x > h.c.l && Math.abs(pos.y - yTheta) < 9) return "theta";
        return null;
    }
    cPaths.addEventListener("pointerdown", function (e) {
        if (e.pointerType === "touch") return;                 // on touch screens the sliders are used instead
        const k = nearHandle(local(cPaths, e));
        if (!k) return;
        drag = k; frozenRange = run ? run.range : rateRange(modelParams());
        cPaths.setPointerCapture(e.pointerId); showTip("tipPaths", cPaths, e, null);
        e.preventDefault();
    });
    cPaths.addEventListener("pointermove", function (e) {
        const pos = local(cPaths, e), h = hit.paths; if (!h) return;
        if (drag) {
            const s = sliders[drag], raw = h.a.invY(pos.y);
            const v = Math.max(parseFloat(s.el.min), Math.min(parseFloat(s.el.max), Math.round(raw * 100) / 100));
            setValue(drag, v / 100);
            return;
        }
        cPaths.style.cursor = nearHandle(pos) ? "ns-resize" : "crosshair";
        const t = h.a.invX(pos.x);
        if (t < 0 || t > T_BOND) { hoverT = null; showTip("tipPaths", cPaths, e, null); drawPaths(modelParams()); return; }
        hoverT = t;
        const p = modelParams(); let html;
        if (run && run.fan && t <= run.reveal * T_BOND) {
            const f = run.fan, i = Math.min(f.step.length - 1, Math.round(t * STEPS_PER_YEAR / 4));
            html = "<b>After " + t.toFixed(1) + " years</b><br>Median rate " + pct(f.q50[i]) + "<br>Half of paths: " + pct(f.q25[i]) + " to " + pct(f.q75[i]) + "<br>90% of paths: " + pct(f.q05[i]) + " to " + pct(f.q95[i]);
        } else {
            const m = rateMean(t, p), sd = Z_ALPHA * rateSd(t, p);
            html = "<b>After " + t.toFixed(1) + " years</b><br>Expected rate " + pct(m) + "<br>90% range: " + pct(m - sd) + " to " + pct(m + sd) + "<br><i>from the exact formulas</i>";
        }
        drawPaths(p); showTip("tipPaths", cPaths, e, html);
    });
    function endDrag() { if (drag) { drag = null; frozenRange = null; render(); } }
    cPaths.addEventListener("pointerup", endDrag);
    cPaths.addEventListener("pointercancel", endDrag);
    cPaths.addEventListener("pointerleave", function (e) { hoverT = null; showTip("tipPaths", cPaths, e, null); if (!drag) drawPaths(modelParams()); });

    const cConv = $("cConv");
    cConv.addEventListener("pointermove", function (e) {
        const h = hit.conv; if (!h) return showTip("tipConv", cConv, e, "Dashed line: the exact price. The simulated estimate appears here after Start.");
        const nAt = h.a.invX(local(cConv, e).x);
        let best = h.pts[0];
        h.pts.forEach(function (q) { if (Math.abs(q[0] - nAt) < Math.abs(best[0] - nAt)) best = q; });
        showTip("tipConv", cConv, e, "<b>After " + best[0].toLocaleString("en-US") + " paths</b><br>Estimate " + best[1].toFixed(3) + " ± " + best[2].toFixed(3) + "<br>Exact " + h.exact.toFixed(3));
    });
    cConv.addEventListener("pointerleave", function (e) { showTip("tipConv", cConv, e, null); });

    const cLoss = $("cLoss");
    cLoss.addEventListener("pointermove", function (e) {
        const h = hit.loss; if (!h) return showTip("tipLoss", cLoss, e, null);
        const b = Math.floor((h.a.invX(local(cLoss, e).x) - h.lo) / h.width);
        if (b < 0 || b >= h.bins) return showTip("tipLoss", cLoss, e, null);
        const from = h.lo + b * h.width, to = from + h.width, share = (h.counts[b] / h.n * 100).toFixed(1);
        const what = to <= 0 ? "a gain of " + Math.abs(to).toFixed(2) + " to " + Math.abs(from).toFixed(2)
            : from >= 0 ? "a loss of " + from.toFixed(2) + " to " + to.toFixed(2)
                : "between a gain of " + Math.abs(from).toFixed(2) + " and a loss of " + to.toFixed(2);
        showTip("tipLoss", cLoss, e, "<b>" + share + "% of scenarios</b><br>end the year with " + what);
    });
    cLoss.addEventListener("pointerleave", function (e) { showTip("tipLoss", cLoss, e, null); });

    const cFlow = $("cFlow");
    cFlow.addEventListener("pointermove", function (e) {
        const h = hit.flow; if (!h) return;
        const j = Math.round(h.a.invX(local(cFlow, e).x) * COUPONS_PER_YEAR) - 1;
        if (j < 0 || j >= h.df.length) return showTip("tipFlow", cFlow, e, null);
        const amt = h.cf.amounts[j];
        showTip("tipFlow", cFlow, e, "<b>Payment in " + h.cf.times[j] + " years</b><br>" + amt.toFixed(2) + " paid then is worth " + (amt * h.df[j]).toFixed(2) + " today<br>One rupee then = " + h.df[j].toFixed(3) + " today");
    });
    cFlow.addEventListener("pointerleave", function (e) { showTip("tipFlow", cFlow, e, null); });

    $("refPrice").textContent = REF.zeroPrice.toFixed(2);
    $("refTol").textContent = String(REF.tolerance);
    window.addEventListener("resize", render);

    /* ----- Light / dark switch (white is the default; the choice is remembered on this device) ----- */
    const themeButton = $("themeToggle");
    function applyTheme(theme) {
        if (theme === "dark") document.documentElement.setAttribute("data-theme", "dark");
        else document.documentElement.removeAttribute("data-theme");
        themeButton.textContent = theme === "dark" ? "Light mode" : "Dark mode";
        themeButton.setAttribute("aria-pressed", String(theme === "dark"));
    }
    themeButton.addEventListener("click", function () {
        const next = document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark";
        applyTheme(next);
        try { localStorage.setItem("theme", next); } catch (err) { /* storage unavailable: the choice lasts until the page closes */ }
        render();
    });
    let savedTheme = "light";
    try { savedTheme = localStorage.getItem("theme") || "light"; } catch (err) { savedTheme = "light"; }
    applyTheme(savedTheme);
    loadPreset("calm");
})();