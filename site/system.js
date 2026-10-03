/* Respirocyte system model page (M29-01).
 *
 * Talks to /api/system (GET for the defaults, POST with edits). Draws the
 * real-size hull with its rotor banks as bands of their measured area, the
 * tank and the cascade chambers as spheres of their volumes, the fill and
 * discharge time course, and the indicators of the central, pessimistic and
 * optimistic cases. Every number comes from the server; this file only lays
 * it out. All results are simulated estimates.
 */
(function () {
  "use strict";

  var COLORS = { load: "#4f9cf0", water: "#3fc1b0", discharge: "#f0a04f", cascade: "#b48cf0" };
  var GROUP_FIELDS = {
    hull: ["diameter_m", "hull_fill_max"],
    tank: ["target_atm", "volume_m3", "cushion_atm", "chamber_volume_m3"],
    power: ["peak_w"],
  };
  var BANK_KEY_PARAMS = ["s_back_nm3", "s_fwd_nm3", "binding_rate_per_s", "energy_per_pass_kj_mol",
    "water_per_pass", "max_rev_per_s"];

  var edits = {};
  var report = null;
  var view = { yaw: 0.6, pitch: 0.35 };

  function $(id) { return document.getElementById(id); }

  function fmt(x, digits) {
    if (x === null || x === undefined) return "n/a";
    if (typeof x !== "number") return String(x);
    if (x === 0) return "0";
    var a = Math.abs(x);
    if (a >= 1e4 || a < 1e-2) return x.toExponential(digits === undefined ? 2 : digits);
    return x.toFixed(digits === undefined ? (a >= 100 ? 0 : a >= 10 ? 1 : 2) : digits);
  }

  function tag(p) {
    return '<span class="tag tag-' + p.tag + '" title="' + (p.source || "no source") +
      (p.note ? " - " + p.note : "") + '">' + p.tag + "</span>";
  }

  function compute() {
    $("indicators").textContent = "Computing ...";
    return fetch("/api/system", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(edits),
    }).then(function (r) { return r.json(); }).then(function (data) {
      if (data.error) throw new Error(data.error);
      report = data;
      render();
    }).catch(function (e) {
      $("indicators").innerHTML = '<span class="bad">' + e.message + "</span>";
    });
  }

  function setGroup(group, name, value) {
    edits.groups = edits.groups || {};
    edits.groups[group] = edits.groups[group] || {};
    edits.groups[group][name] = value;
  }

  function setBank(kind, stage, field, value) {
    edits.banks = edits.banks || [];
    var e = edits.banks.find(function (b) { return b.kind === kind && b.stage === stage; });
    if (!e) { e = { kind: kind, stage: stage }; edits.banks.push(e); }
    e[field] = value;
  }

  function setBankParam(kind, name, value) {
    edits.bank_params = edits.bank_params || {};
    edits.bank_params[kind] = edits.bank_params[kind] || {};
    edits.bank_params[kind][name] = value;
  }

  function input(value, onchange) {
    var el = document.createElement("input");
    el.type = "number";
    el.value = value;
    el.step = "any";
    el.addEventListener("change", function () { onchange(parseFloat(el.value)); });
    return el;
  }

  function row(label, el, extra) {
    var r = document.createElement("div");
    r.className = "row";
    var l = document.createElement("label");
    l.innerHTML = label + (extra || "");
    r.appendChild(l);
    r.appendChild(el);
    return r;
  }

  function renderDesign() {
    var dev = report.device;
    $("stages").value = dev.n_stages;
    var g = $("groups");
    g.innerHTML = "";
    Object.keys(GROUP_FIELDS).forEach(function (group) {
      GROUP_FIELDS[group].forEach(function (name) {
        var p = dev.groups[group][name];
        g.appendChild(row(group + "." + name + " <span class='hint'>" + p.unit + "</span>",
          input(p.value, function (v) { setGroup(group, name, v); }), tag(p)));
      });
    });
    var bs = $("banks");
    bs.innerHTML = "";
    dev.banks.forEach(function (b) {
      var box = document.createElement("div");
      box.className = "bank";
      var head = document.createElement("div");
      head.innerHTML = '<span class="swatch" style="background:' + (COLORS[b.kind] || "#999") + '"></span><b>' +
        b.kind + "</b>" + (b.kind === "load" || b.kind === "cascade" ? " stage " + b.stage : "") +
        ' <span class="hint">(' + b.placement + ")</span>";
      box.appendChild(head);
      box.appendChild(row("sections", input(b.count, function (v) { setBank(b.kind, b.stage, "count", v); })));
      box.appendChild(row("design rev/s", input(b.rev_per_s, function (v) { setBank(b.kind, b.stage, "rev_per_s", v); })));
      var det = document.createElement("details");
      det.innerHTML = "<summary class='hint'>parameters (all " + b.kind + " banks)</summary>";
      BANK_KEY_PARAMS.forEach(function (name) {
        var p = b.params[name];
        if (!p) return;
        det.appendChild(row(name + " <span class='hint'>" + p.unit + "</span>",
          input(p.value, function (v) { setBankParam(b.kind, name, v); }), tag(p)));
      });
      box.appendChild(det);
      bs.appendChild(box);
    });
  }

  function renderIndicators() {
    var c = report.cases;
    var rows = [
      ["Hull used", function (x) { return fmt(100 * x.budget.hull_fraction_used, 1) + " %"; }],
      ["Sections that fit", function (x) { return fmt(x.budget.sections_that_fit); }],
      ["Tank capacity (O2)", function (x) { return fmt(x.budget.tank_capacity_o2_at_target); }],
      ["Cascade ceiling (atm)", function (x) { return fmt(x.budget.cascade_ceiling_atm); }],
      ["Top-stage limit, net O2/pass", function (x) { return fmt(x.budget.infinite_cascade_net_at_target); }],
      ["Fill result", function (x) { return x.events.fill_end; }],
      ["Fill time (s)", function (x) { return fmt(x.events.fill_time_s); }],
      ["Tank after fill (atm)", function (x) { return fmt(x.events.tank_atm_after_fill); }],
      ["Time to 500 atm (s)", function (x) { return fmt(x.events.t_tank_500_atm_s); }],
      ["Discharge time (s)", function (x) { return fmt(x.events.discharge_time_s); }],
      ["O2 delivered", function (x) { return fmt(x.events.o2_delivered); }],
      ["Power demand at design (W)", function (x) { return fmt(x.budget.power_demand_at_design_w); }],
      ["Power peak (W)", function (x) { return fmt(x.budget.power_peak_w); }],
      ["Lowest duty cycle", function (x) { return fmt(x.events.min_duty); }],
      ["Mean power (W)", function (x) { return fmt(x.events.mean_power_w); }],
      ["Energy per delivered O2 (J)", function (x) { return fmt(x.events.energy_per_delivered_o2_j); }],
      ["Water peak in sump", function (x) { return fmt(x.events.water_peak); }],
      ["Sump fill at peak", function (x) { return fmt(100 * x.events.water_peak_fraction_of_chamber, 1) + " %"; }],
    ];
    var html = "<table class='ind'><tr><th></th><th>central</th><th>pessim.</th><th>optim.</th></tr>";
    rows.forEach(function (r) {
      html += "<tr><td>" + r[0] + "</td>" + ["central", "pessimistic", "optimistic"].map(function (k) {
        return "<td>" + r[1](c[k]) + "</td>";
      }).join("") + "</tr>";
    });
    html += "</table><div class='hint'>Pessimistic and optimistic put every uncertain value at the end of " +
      "its range that hurts or helps, power supply included. Fill from the 10 atm cushion; discharge to it.</div>";
    $("indicators").innerHTML = html;
    var w = report.warnings || [];
    var b = c.central.budget;
    var lines = w.map(function (s) { return '<div class="warn">' + s + "</div>"; });
    if (!b.hull_fits) lines.push('<div class="bad">The hull banks take more than the allowed area.</div>');
    $("warnings").innerHTML = lines.join("");
  }

  function renderProvenance() {
    $("provenance").innerHTML = report.provenance.map(function (p) {
      return "<div class='prov-row'><b>" + p.component + "." + p.name + "</b> = " + fmt(p.value) + " " + p.unit +
        " [" + fmt(p.low) + ", " + fmt(p.high) + "] " + tag(p) +
        "<div class='src'>" + (p.source || "no source") + (p.note ? " - " + p.note : "") + "</div></div>";
    }).join("");
  }

  // --- Hull drawing: orthographic, drag to turn. ---
  function project(p, cx, cy, s) {
    var cyaw = Math.cos(view.yaw), syaw = Math.sin(view.yaw);
    var cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
    var x = cyaw * p[0] + syaw * p[1];
    var y = -syaw * p[0] + cyaw * p[1];
    var z = p[2];
    var y2 = cp * y - sp * z;
    var z2 = sp * y + cp * z;
    return [cx + s * x, cy - s * z2, y2];
  }

  function sph(r, th, ph) {
    return [r * Math.sin(th) * Math.cos(ph), r * Math.sin(th) * Math.sin(ph), r * Math.cos(th)];
  }

  function drawHull() {
    var cv = $("hull");
    var w = cv.clientWidth, h = cv.clientHeight;
    cv.width = w; cv.height = h;
    var ctx = cv.getContext("2d");
    ctx.fillStyle = "#0b0d11";
    ctx.fillRect(0, 0, w, h);
    var L = report.layout;
    var R = L.hull_radius_m;
    var s = 0.42 * Math.min(w, h) / R;
    var cx = w / 2, cy = h / 2;
    var quads = [];
    var nph = 48;
    L.bands.forEach(function (b) {
      var nth = Math.max(1, Math.ceil((b.theta1_rad - b.theta0_rad) / 0.12));
      for (var i = 0; i < nth; i++) {
        var t0 = b.theta0_rad + (b.theta1_rad - b.theta0_rad) * i / nth;
        var t1 = b.theta0_rad + (b.theta1_rad - b.theta0_rad) * (i + 1) / nth;
        for (var j = 0; j < nph; j++) {
          var p0 = 2 * Math.PI * j / nph, p1 = 2 * Math.PI * (j + 1) / nph;
          var pts = [sph(R, t0, p0), sph(R, t0, p1), sph(R, t1, p1), sph(R, t1, p0)].map(function (q) {
            return project(q, cx, cy, s);
          });
          var depth = (pts[0][2] + pts[2][2]) / 2;
          quads.push({ pts: pts, depth: depth, color: COLORS[b.kind] || "#999" });
        }
      }
    });
    quads.sort(function (a, b) { return a.depth - b.depth; });
    // Back half first, dim.
    quads.forEach(function (q) {
      if (q.depth > 0) return;
      ctx.globalAlpha = 0.18;
      ctx.fillStyle = q.color;
      ctx.beginPath();
      q.pts.forEach(function (p, k) { if (k) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]); });
      ctx.fill();
    });
    // Interior: tank and chambers.
    ctx.globalAlpha = 0.35;
    var t = project(L.tank.center_m, cx, cy, s);
    ctx.fillStyle = "#6a7d99";
    ctx.beginPath(); ctx.arc(t[0], t[1], L.tank.radius_m * s, 0, 2 * Math.PI); ctx.fill();
    ctx.globalAlpha = 0.9;
    ctx.fillStyle = "#cfd7e3";
    ctx.font = "11px system-ui";
    ctx.fillText("O2 tank " + fmt(report.device.groups.tank.volume_m3.value) + " m^3", t[0] - 50, t[1]);
    L.chambers.forEach(function (c) {
      var q = project(c.center_m, cx, cy, s);
      ctx.globalAlpha = 0.6;
      ctx.fillStyle = COLORS.cascade;
      ctx.beginPath(); ctx.arc(q[0], q[1], Math.max(2, c.radius_m * s), 0, 2 * Math.PI); ctx.fill();
      ctx.globalAlpha = 0.9;
      ctx.fillStyle = "#cfd7e3";
      ctx.fillText("chamber " + c.stage, q[0] + 6, q[1] - 4);
    });
    // Hull outline and the front half of the bands.
    ctx.globalAlpha = 1.0;
    ctx.strokeStyle = "#5c6678";
    ctx.beginPath(); ctx.arc(cx, cy, R * s, 0, 2 * Math.PI); ctx.stroke();
    quads.forEach(function (q) {
      if (q.depth <= 0) return;
      ctx.globalAlpha = 0.75;
      ctx.fillStyle = q.color;
      ctx.beginPath();
      q.pts.forEach(function (p, k) { if (k) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]); });
      ctx.fill();
    });
    ctx.globalAlpha = 1.0;
    // Legend and scale bar.
    var y = 16;
    ctx.font = "12px system-ui";
    L.bands.forEach(function (b) {
      ctx.fillStyle = COLORS[b.kind]; ctx.fillRect(10, y - 9, 10, 10);
      ctx.fillStyle = "#d8dde6";
      ctx.fillText(b.kind + ": " + fmt(b.count) + " sections, " + fmt(b.area_m2) + " m^2", 26, y);
      y += 16;
    });
    L.internal_banks.forEach(function (b) {
      ctx.fillStyle = COLORS.cascade; ctx.fillRect(10, y - 9, 10, 10);
      ctx.fillStyle = "#d8dde6";
      ctx.fillText("stage " + b.stage + " (internal bulkhead): " + fmt(b.count) + " sections", 26, y);
      y += 16;
    });
    var bar = 100e-9 * s;
    ctx.strokeStyle = "#d8dde6";
    ctx.beginPath(); ctx.moveTo(10, h - 14); ctx.lineTo(10 + bar, h - 14); ctx.stroke();
    ctx.fillStyle = "#d8dde6";
    ctx.fillText("100 nm (hull diameter " + fmt(2 * R * 1e9) + " nm; one rotor section is 13.7 x 6.9 nm)",
      10, h - 20);
  }

  function drawPlot() {
    var cv = $("plot");
    var w = cv.clientWidth, h = cv.clientHeight;
    cv.width = w; cv.height = h;
    var ctx = cv.getContext("2d");
    ctx.fillStyle = "#0b0d11";
    ctx.fillRect(0, 0, w, h);
    var log = report.cases.central.log;
    if (!log || !log.t_s.length) return;
    var pad = { l: 48, r: 44, t: 14, b: 22 };
    var tmax = log.t_s[log.t_s.length - 1];
    var pmax = 0;
    log.tank_atm.forEach(function (p, i) {
      pmax = Math.max(pmax, p);
      log.chamber_atm[i].forEach(function (q) { pmax = Math.max(pmax, q); });
    });
    pmax = Math.max(pmax * 1.05, 10);
    function X(t) { return pad.l + (w - pad.l - pad.r) * t / tmax; }
    function Y(p) { return h - pad.b - (h - pad.t - pad.b) * p / pmax; }
    function Yd(d) { return h - pad.b - (h - pad.t - pad.b) * d; }
    // Discharge shading.
    var i0 = log.phase.indexOf("discharge");
    if (i0 >= 0) {
      ctx.fillStyle = "#1c1a12";
      ctx.fillRect(X(log.t_s[i0]), pad.t, X(tmax) - X(log.t_s[i0]), h - pad.t - pad.b);
    }
    ctx.strokeStyle = "#3a4250";
    ctx.strokeRect(pad.l, pad.t, w - pad.l - pad.r, h - pad.t - pad.b);
    ctx.fillStyle = "#8b95a7";
    ctx.font = "11px system-ui";
    [0, 0.25, 0.5, 0.75, 1].forEach(function (f) {
      ctx.fillText(fmt(f * pmax), 4, Y(f * pmax) + 4);
      ctx.fillText(fmt(f * tmax) + " s", X(f * tmax) - 10, h - 6);
    });
    ctx.fillText("atm", pad.l + 4, pad.t + 26);
    ctx.fillText("duty", w - pad.r + 6, pad.t + 8);
    function line(vals, color, width, ys) {
      ctx.strokeStyle = color; ctx.lineWidth = width; ctx.beginPath();
      vals.forEach(function (v, k) {
        var x = X(log.t_s[k]), y = (ys || Y)(v);
        if (k) ctx.lineTo(x, y); else ctx.moveTo(x, y);
      });
      ctx.stroke();
    }
    var nch = log.chamber_atm[0].length;
    for (var c = 0; c < nch; c++) {
      line(log.chamber_atm.map(function (a) { return a[c]; }), COLORS.cascade, 1);
    }
    line(log.duty, "#5a6474", 1, Yd);
    line(log.tank_atm, "#f2f4f8", 2);
    ctx.lineWidth = 1;
    ctx.fillStyle = "#f2f4f8";
    ctx.fillText("tank pressure (white), cascade chambers (violet), duty cycle (grey); shaded: discharge in the tissues",
      pad.l + 6, pad.t + 12);
  }

  function render() {
    renderDesign();
    renderIndicators();
    renderProvenance();
    drawHull();
    drawPlot();
  }

  function init() {
    var cv = $("hull");
    var drag = null;
    cv.addEventListener("mousedown", function (e) { drag = [e.clientX, e.clientY, view.yaw, view.pitch]; });
    window.addEventListener("mouseup", function () { drag = null; });
    window.addEventListener("mousemove", function (e) {
      if (!drag || !report) return;
      view.yaw = drag[2] + (e.clientX - drag[0]) * 0.01;
      view.pitch = Math.max(-1.5, Math.min(1.5, drag[3] + (e.clientY - drag[1]) * 0.01));
      drawHull();
    });
    window.addEventListener("resize", function () { if (report) { drawHull(); drawPlot(); } });
    $("stages").addEventListener("change", function () { edits.n_stages = parseInt($("stages").value, 10); });
    $("run").addEventListener("click", compute);
    $("reset").addEventListener("click", function () { edits = {}; compute(); });
    compute();
  }

  window.NanoCadSystem = {
    report: function () { return report; },
    edits: function () { return edits; },
    compute: function (e) { edits = e || {}; return compute(); },
  };
  document.addEventListener("DOMContentLoaded", init);
})();
