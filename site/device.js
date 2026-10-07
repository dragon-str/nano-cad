/* Device scale explorer (M29-11).
 *
 * One page, one continuous zoom, from the whole respirocyte (about 1 um)
 * down to the carrier heme site (about 2 nm). The hull and its banks come
 * from /api/system; the atomistic detail comes from
 * /api/part?id=carrier_site. The camera pans from the hull centre to the
 * site as the view shrinks, and each level fades in or out by its size in
 * pixels. Everything is drawn from server data; this file only lays it out.
 * Simulated geometry.
 */
(function () {
  "use strict";

  var BANK_COLORS = { load: "#4f9cf0", water: "#3fc1b0", discharge: "#f0a04f", cascade: "#b48cf0", shell: "#39424f" };
  var ELEMENT_COLORS = {
    H: "#e8eef5", C: "#8b98a5", N: "#5b8def", O: "#f05a5a",
    F: "#7ee787", Si: "#c9a227", P: "#ff9e64", S: "#e0c000", B: "#ffb86c", Fe: "#e0662c",
  };
  var ATOM_RADII = { H: 0.31e-10, C: 0.76e-10, N: 0.71e-10, O: 0.66e-10, Fe: 1.32e-10 };

  // The zoom ladder, in log10 of the view width in metres.
  var LOG_WIDE = -4.7;   // 20 um: the whole hull, small
  var LOG_TIGHT = -9.7;  // 0.2 nm: one atom
  var HULL_FADE = [-8.4, -6.4];   // hull alpha 0 .. 1
  var SITE_FADE = [-6.2, -7.4];   // site alpha 0 .. 1
  var PAN_FROM = -5.2, PAN_TO = -7.0;  // camera pans from hull centre to the site

  var system = null;
  var part = null;
  var cam = { yaw: 0.6, pitch: 0.32, logWidth: -5.15, drag: null };
  var sitePos = null;
  var siteCenter = [0, 0, 0];
  var bonds = [];
  var dirty = true;

  function $(id) { return document.getElementById(id); }
  function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
  function smoothstep(a, b, x) { var t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); }
  function lerp(a, b, t) { return a + (b - a) * t; }

  function fmtLen(m) {
    var a = Math.abs(m);
    if (a >= 1e-6) return (m * 1e6).toPrecision(3) + " um";
    if (a >= 1e-9) return (m * 1e9).toPrecision(3) + " nm";
    return (m * 1e10).toPrecision(3) + " A";
  }

  function elementColor(e) { return ELEMENT_COLORS[e] || "#b9c2cf"; }
  function atomRadius(e) { return ATOM_RADII[e] || 0.7e-10; }
  function hexToRgb(hex) {
    var n = parseInt(hex.slice(1), 16);
    return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
  }

  function load() {
    return Promise.all([
      fetch("/api/system").then(function (r) { return r.json(); }),
      fetch("/api/part?id=carrier_site").then(function (r) { return r.json(); }),
    ]).then(function (res) {
      system = res[0];
      part = res[1].scene;
      sitePos = [0, 0, system.layout.hull_radius_m];
      var atoms = part.atomistic.atoms;
      var lo = [1e9, 1e9, 1e9], hi = [-1e9, -1e9, -1e9];
      atoms.forEach(function (a) {
        for (var k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], a.position_m[k]); hi[k] = Math.max(hi[k], a.position_m[k]); }
      });
      for (var k = 0; k < 3; k++) siteCenter[k] = sitePos[k] + 0.5 * (lo[k] + hi[k]);
      buildBonds(atoms);
      dirty = true;
    });
  }

  function buildBonds(atoms) {
    var n = atoms.length, p = atoms.map(function (a) { return a.position_m; });
    var cut = 1.75e-10;
    for (var i = 0; i < n; i++) {
      for (var j = i + 1; j < n; j++) {
        var dx = p[i][0] - p[j][0], dy = p[i][1] - p[j][1], dz = p[i][2] - p[j][2];
        if (dx * dx + dy * dy + dz * dz < cut * cut) bonds.push([i, j]);
      }
    }
  }

  function viewState() {
    var cv = $("stage");
    var W = cv.clientWidth, H = cv.clientHeight;
    var s = W / Math.pow(10, cam.logWidth);
    var focus = smoothstep(PAN_FROM, PAN_TO, cam.logWidth);
    var c = [0, 0, 0];
    for (var k = 0; k < 3; k++) c[k] = lerp(0, siteCenter[k], focus);
    return { W: W, H: H, s: s, c: c, focus: focus };
  }

  function makeProjector(v) {
    var cy = Math.cos(cam.yaw), sy = Math.sin(cam.yaw);
    var cp = Math.cos(cam.pitch), sp = Math.sin(cam.pitch);
    return function (p) {
      var dx = p[0] - v.c[0], dy = p[1] - v.c[1], dz = p[2] - v.c[2];
      var x = cy * dx + sy * dy;
      var y = -sy * dx + cy * dy;
      var y2 = cp * y - sp * dz;
      var z2 = sp * y + cp * dz;
      return [v.W / 2 + v.s * x, v.H / 2 - v.s * z2, y2];
    };
  }

  function sph(r, th, ph) {
    return [r * Math.sin(th) * Math.cos(ph), r * Math.sin(th) * Math.sin(ph), r * Math.cos(th)];
  }

  function drawHull(ctx, v, pr, alpha) {
    var L = system.layout, R = L.hull_radius_m;
    if (R * v.s < 6) return;  // smaller than a few pixels
    var nph = 96, dth = 0.06, quads = [];
    var bands = L.bands.slice(), maxT = 0;
    L.bands.forEach(function (b) { if (b.theta1_rad > maxT) maxT = b.theta1_rad; });
    if (maxT < Math.PI - 0.02) bands.push({ kind: "shell", theta0_rad: maxT, theta1_rad: Math.PI });
    bands.forEach(function (b) {
      var nth = Math.max(1, Math.ceil((b.theta1_rad - b.theta0_rad) / dth));
      for (var i = 0; i < nth; i++) {
        var t0 = b.theta0_rad + (b.theta1_rad - b.theta0_rad) * i / nth;
        var t1 = b.theta0_rad + (b.theta1_rad - b.theta0_rad) * (i + 1) / nth;
        for (var j = 0; j < nph; j++) {
          var p0 = 2 * Math.PI * j / nph, p1 = 2 * Math.PI * (j + 1) / nph;
          var q = [sph(R, t0, p0), sph(R, t0, p1), sph(R, t1, p1), sph(R, t1, p0)].map(pr);
          var minx = Math.min(q[0][0], q[1][0], q[2][0], q[3][0]);
          var maxx = Math.max(q[0][0], q[1][0], q[2][0], q[3][0]);
          var miny = Math.min(q[0][1], q[1][1], q[2][1], q[3][1]);
          var maxy = Math.max(q[0][1], q[1][1], q[2][1], q[3][1]);
          if (maxx < -20 || minx > v.W + 20 || maxy < -20 || miny > v.H + 20) continue;
          quads.push({ pts: q, depth: (q[0][2] + q[2][2]) / 2, color: BANK_COLORS[b.kind] || "#8b98a5" });
        }
      }
    });
    quads.sort(function (a, b) { return a.depth - b.depth; });
    ctx.globalAlpha = alpha * 0.16;
    quads.forEach(function (q) {
      if (q.depth > 0) return;
      ctx.fillStyle = q.color;
      fillPoly(ctx, q.pts);
    });
    // The tank and the cascade chambers, behind the front surface.
    ctx.globalAlpha = alpha * 0.35;
    var t = pr(L.tank.center_m);
    ctx.fillStyle = "#6a7d99";
    ctx.beginPath(); ctx.arc(t[0], t[1], L.tank.radius_m * v.s, 0, 2 * Math.PI); ctx.fill();
    L.chambers.forEach(function (c) {
      var q = pr(c.center_m);
      ctx.fillStyle = BANK_COLORS.cascade;
      ctx.beginPath(); ctx.arc(q[0], q[1], Math.max(1.5, c.radius_m * v.s), 0, 2 * Math.PI); ctx.fill();
    });
    ctx.globalAlpha = alpha;
    ctx.strokeStyle = "rgba(120,132,150,0.55)";
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(v.W / 2 - v.s * (v.c[0] * Math.cos(cam.yaw) + v.c[1] * Math.sin(cam.yaw)),
      v.H / 2, R * v.s, 0, 2 * Math.PI); ctx.stroke();
    quads.forEach(function (q) {
      if (q.depth <= 0) return;
      ctx.fillStyle = q.color;
      ctx.globalAlpha = alpha * 0.8;
      fillPoly(ctx, q.pts);
    });
    ctx.globalAlpha = 1;
  }

  function fillPoly(ctx, pts) {
    ctx.beginPath();
    pts.forEach(function (p, k) { if (k) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]); });
    ctx.fill();
  }

  function drawSite(ctx, v, pr, alpha) {
    var atoms = part.atomistic.atoms;
    var spacing = 1.54e-10 * v.s;
    var draws = [];
    for (var i = 0; i < atoms.length; i++) {
      var a = atoms[i];
      var w = [sitePos[0] + a.position_m[0], sitePos[1] + a.position_m[1], sitePos[2] + a.position_m[2]];
      var q = pr(w);
      if (q[0] < -30 || q[0] > v.W + 30 || q[1] < -30 || q[1] > v.H + 30) continue;
      draws.push({ x: q[0], y: q[1], d: q[2], r: Math.max(0.6, 0.42 * atomRadius(a.element) * v.s), c: elementColor(a.element) });
    }
    draws.sort(function (a, b) { return a.d - b.d; });
    if (spacing > 3.5) {
      ctx.globalAlpha = alpha * 0.5;
      ctx.strokeStyle = "#7d8798";
      ctx.lineWidth = Math.max(0.5, 0.14e-10 * v.s);
      ctx.beginPath();
      bonds.forEach(function (b) {
        var p0 = pr([sitePos[0] + atoms[b[0]].position_m[0], sitePos[1] + atoms[b[0]].position_m[1], sitePos[2] + atoms[b[0]].position_m[2]]);
        var p1 = pr([sitePos[0] + atoms[b[1]].position_m[0], sitePos[1] + atoms[b[1]].position_m[1], sitePos[2] + atoms[b[1]].position_m[2]]);
        ctx.moveTo(p0[0], p0[1]); ctx.lineTo(p1[0], p1[1]);
      });
      ctx.stroke();
    }
    draws.forEach(function (a) {
      ctx.globalAlpha = alpha;
      ctx.fillStyle = a.c;
      ctx.beginPath(); ctx.arc(a.x, a.y, a.r, 0, 2 * Math.PI); ctx.fill();
      if (a.r > 3) {
        ctx.globalAlpha = alpha * 0.5;
        ctx.fillStyle = "#0b0d11";
        ctx.beginPath(); ctx.arc(a.x - a.r * 0.3, a.y - a.r * 0.3, a.r * 0.5, 0, 2 * Math.PI); ctx.fill();
      }
    });
    ctx.globalAlpha = 1;
  }

  function draw() {
    dirty = false;
    var cv = $("stage");
    var W = cv.clientWidth, H = cv.clientHeight;
    if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
    var ctx = cv.getContext("2d");
    ctx.fillStyle = "#0b0d11";
    ctx.fillRect(0, 0, W, H);
    if (!system || !part) return;
    var v = viewState(), pr = makeProjector(v);
    var hullA = smoothstep(HULL_FADE[0], HULL_FADE[1], cam.logWidth);
    var siteA = smoothstep(SITE_FADE[0], SITE_FADE[1], cam.logWidth);
    drawHull(ctx, v, pr, hullA);
    drawSite(ctx, v, pr, siteA);
    drawOverlay(ctx, v);
  }

  function drawOverlay(ctx, v) {
    var W = v.W, H = v.H;
    // Scale bar.
    var targetPx = 120, metres = targetPx / v.s;
    var nice = [1, 2, 5, 10, 20, 50, 100, 200, 500];
    var unit = 1e-9, label = "nm";
    if (metres >= 1e-6) { unit = 1e-6; label = "um"; }
    else if (metres < 1e-9) { unit = 1e-10; label = "A"; }
    var val = metres / unit, pick = nice[0];
    nice.forEach(function (n) { if (n <= val) pick = n; });
    var bar = pick * unit * v.s;
    ctx.strokeStyle = "#d8dde6"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(16, H - 58); ctx.lineTo(16 + bar, H - 58); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(16, H - 63); ctx.lineTo(16, H - 53); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(16 + bar, H - 63); ctx.lineTo(16 + bar, H - 53); ctx.stroke();
    ctx.fillStyle = "#d8dde6"; ctx.font = "12px system-ui";
    ctx.fillText(pick + " " + label, 16, H - 64);
    // Level label.
    var level = v.focus < 0.5 ? "Respirocyte" : "Carrier site";
    ctx.font = "bold 15px system-ui";
    ctx.fillText(level, 16, 26);
    ctx.font = "12px system-ui";
    ctx.fillStyle = "#9aa4b4";
    var sub = v.focus < 0.5
      ? "hull " + fmtLen(2 * system.layout.hull_radius_m) + " across; view " + fmtLen(Math.pow(10, cam.logWidth))
      : "Fe-porphine carrier, " + part.atomistic.atom_count + " atoms; view " + fmtLen(Math.pow(10, cam.logWidth));
    ctx.fillText(sub, 16, 44);
    // Legend.
    if (v.focus < 0.5) {
      var y = 70;
      system.layout.bands.forEach(function (b) {
        ctx.fillStyle = BANK_COLORS[b.kind] || "#8b98a5";
        ctx.fillRect(16, y - 9, 10, 10);
        ctx.fillStyle = "#d8dde6";
        ctx.fillText(b.kind + ": " + b.count.toLocaleString() + " sections", 32, y);
        y += 17;
      });
      ctx.fillStyle = "#9aa4b4";
      ctx.fillText("the carrier site is one molecule; the hull holds about 1e20", 16, y + 4);
    } else {
      var y2 = 70, seen = {};
      part.atomistic.atoms.forEach(function (a) { seen[a.element] = (seen[a.element] || 0) + 1; });
      Object.keys(seen).sort().forEach(function (e) {
        ctx.fillStyle = elementColor(e);
        ctx.beginPath(); ctx.arc(21, y2 - 4, 6, 0, 2 * Math.PI); ctx.fill();
        ctx.fillStyle = "#d8dde6";
        ctx.fillText(e + " x" + seen[e], 34, y2);
        y2 += 17;
      });
    }
  }

  function animate() {
    if (dirty) draw();
    requestAnimationFrame(animate);
  }

  function init() {
    var cv = $("stage");
    cv.addEventListener("pointerdown", function (e) { cam.drag = [e.clientX, e.clientY, cam.yaw, cam.pitch]; cv.setPointerCapture(e.pointerId); });
    cv.addEventListener("pointermove", function (e) {
      if (!cam.drag) return;
      cam.yaw = cam.drag[2] + (e.clientX - cam.drag[0]) * 0.01;
      cam.pitch = clamp(cam.drag[3] + (e.clientY - cam.drag[1]) * 0.01, -1.5, 1.5);
      dirty = true;
    });
    window.addEventListener("pointerup", function () { cam.drag = null; });
    cv.addEventListener("wheel", function (e) {
      e.preventDefault();
      cam.logWidth = clamp(cam.logWidth + (e.deltaY < 0 ? -0.12 : 0.12), LOG_TIGHT, LOG_WIDE);
      dirty = true;
    }, { passive: false });
    window.addEventListener("keydown", function (e) {
      if (e.key === "=" || e.key === "+") cam.logWidth = clamp(cam.logWidth - 0.2, LOG_TIGHT, LOG_WIDE);
      else if (e.key === "-" || e.key === "_") cam.logWidth = clamp(cam.logWidth + 0.2, LOG_TIGHT, LOG_WIDE);
      else return;
      dirty = true;
    });
    $("zoom-in").addEventListener("click", function () { cam.logWidth = clamp(cam.logWidth - 0.6, LOG_TIGHT, LOG_WIDE); dirty = true; });
    $("zoom-out").addEventListener("click", function () { cam.logWidth = clamp(cam.logWidth + 0.6, LOG_TIGHT, LOG_WIDE); dirty = true; });
    $("go-device").addEventListener("click", function () { cam.logWidth = -5.15; dirty = true; });
    $("go-site").addEventListener("click", function () { cam.logWidth = -9.0; dirty = true; });
    window.addEventListener("resize", function () { dirty = true; });
    animate();
    load();
  }

  window.NanoCadDevice = {
    ready: function () { return !!(system && part); },
    view: function () { return cam; },
    setWidth: function (lw) { cam.logWidth = clamp(lw, LOG_TIGHT, LOG_WIDE); dirty = true; },
    level: function () { return viewState().focus < 0.5 ? "respirocyte" : "carrier_site"; },
    siteCenter: function () { return siteCenter.slice(); },
    redraw: function () { dirty = true; },
    scale: function () { return viewState().s; },
  };
  document.addEventListener("DOMContentLoaded", init);
})();
