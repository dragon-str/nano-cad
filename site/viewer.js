/*
 * Nano-CAD static three-scale viewer.
 *
 * Draws the three layers of a `nanocad.scene` document:
 *   - atomistic (L1): one dot per atom,
 *   - device (L2): one marker per rigid body plus the joints,
 *   - coarse (handoff): gear pitch circles and bounding cylinders.
 *
 * The viewer is offline. It uses the Canvas 2D API only. It has no external
 * dependency. The drawing is a schematic projection, not a physical render.
 */
(function () {
  "use strict";

  var canvas = document.getElementById("canvas");
  var ctx = canvas.getContext("2d");
  var slider = document.getElementById("layer-slider");
  var sliderLabel = document.getElementById("slider-label");
  var readout = document.getElementById("readout");
  var message = document.getElementById("message");
  var fallback = document.getElementById("load-fallback");
  var fileInput = document.getElementById("file-input");

  var toggles = {
    atomistic: document.getElementById("toggle-atomistic"),
    device: document.getElementById("toggle-device"),
    coarse: document.getElementById("toggle-coarse"),
  };

  var ROLE_COLORS = {
    ground: "#6e7681",
    sun: "#f2c14e",
    planet: "#4ec9b0",
    ring: "#b083f0",
    carrier: "#f0883e",
  };

  var ELEMENT_COLORS = {
    H: "#e8eef5",
    C: "#8b98a5",
    N: "#5b8def",
    O: "#f05a5a",
    F: "#7ee787",
    Si: "#c9a227",
    P: "#ff9e64",
    S: "#e0c000",
    B: "#ffb86c",
  };

  var view = { yaw: -0.65, tilt: 1.02, zoom: 1 };
  var fit = { c: [0, 0, 0], r: 1 };
  var scene = null;
  var drag = null;

  /* ---------- Small vector helpers ---------- */

  function sub(a, b) {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  }

  function cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
  }

  function norm(a) {
    var len = Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) || 1;
    return [a[0] / len, a[1] / len, a[2] / len];
  }

  function elementColor(element) {
    if (ELEMENT_COLORS[element]) {
      return ELEMENT_COLORS[element];
    }
    var z = parseInt(element.replace(/[^0-9]/g, ""), 10) || 20;
    return "hsl(" + ((z * 47) % 360) + ", 55%, 60%)";
  }

  /* ---------- Fit and projection ---------- */

  function addExtent(point) {
    for (var i = 0; i < 3; i += 1) {
      if (point[i] < fitMin[i]) {
        fitMin[i] = point[i];
      }
      if (point[i] > fitMax[i]) {
        fitMax[i] = point[i];
      }
    }
  }

  var fitMin = [0, 0, 0];
  var fitMax = [0, 0, 0];

  function computeFit() {
    fitMin = [Infinity, Infinity, Infinity];
    fitMax = [-Infinity, -Infinity, -Infinity];

    var atoms = scene.atomistic.atoms;
    for (var i = 0; i < atoms.length; i += 1) {
      addExtent(atoms[i].position_m);
    }
    for (var b = 0; b < scene.device.bodies.length; b += 1) {
      addExtent(scene.device.bodies[b].position_m);
    }
    for (var k = 0; k < scene.coarse.bodies.length; k += 1) {
      var coarse = scene.coarse.bodies[k];
      if (coarse.pitch_circle) {
        var pc = coarse.pitch_circle;
        addExtent([
          pc.center_m[0] - pc.radius_m,
          pc.center_m[1] - pc.radius_m,
          pc.center_m[2],
        ]);
        addExtent([
          pc.center_m[0] + pc.radius_m,
          pc.center_m[1] + pc.radius_m,
          pc.center_m[2],
        ]);
      }
      if (coarse.bounding_cylinder) {
        var cy = coarse.bounding_cylinder;
        addExtent([
          cy.center_m[0] - cy.radius_m,
          cy.center_m[1] - cy.radius_m,
          cy.center_m[2] - cy.half_length_m,
        ]);
        addExtent([
          cy.center_m[0] + cy.radius_m,
          cy.center_m[1] + cy.radius_m,
          cy.center_m[2] + cy.half_length_m,
        ]);
      }
    }

    fit.c = [
      0.5 * (fitMin[0] + fitMax[0]),
      0.5 * (fitMin[1] + fitMax[1]),
      0.5 * (fitMin[2] + fitMax[2]),
    ];

    var radius = 0;
    for (var q = 0; q < scene.atomistic.atoms.length; q += 1) {
      var d = sub(scene.atomistic.atoms[q].position_m, fit.c);
      radius = Math.max(radius, Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
    }
    for (var s = 0; s < scene.coarse.bodies.length; s += 1) {
      var body = scene.coarse.bodies[s];
      if (body.pitch_circle) {
        var cd = sub(body.pitch_circle.center_m, fit.c);
        radius = Math.max(
          radius,
          Math.sqrt(cd[0] * cd[0] + cd[1] * cd[1] + cd[2] * cd[2]) +
            body.pitch_circle.radius_m
        );
      }
    }
    fit.r = radius > 0 ? radius : 1;
  }

  function pixelSize() {
    var w = canvas.clientWidth || 800;
    var h = canvas.clientHeight || 600;
    return Math.min(w, h) * 0.42 * view.zoom;
  }

  function camera(point) {
    var cy = Math.cos(view.yaw);
    var sy = Math.sin(view.yaw);
    var ct = Math.cos(view.tilt);
    var st = Math.sin(view.tilt);

    var nx = (point[0] - fit.c[0]) / fit.r;
    var ny = (point[1] - fit.c[1]) / fit.r;
    var nz = (point[2] - fit.c[2]) / fit.r;

    var x1 = nx * cy - ny * sy;
    var y1 = nx * sy + ny * cy;
    var y2 = y1 * ct - nz * st;
    var depth = y1 * st + nz * ct;

    var camDistance = 3.2;
    var persp = camDistance / (camDistance - depth * 0.6);
    var scale = pixelSize();

    return [w / 2 + x1 * scale * persp, h / 2 + y2 * scale * persp, depth, persp];
  }

  function pixRadius(radius_m, persp) {
    return (radius_m / fit.r) * pixelSize() * persp;
  }

  var w = 0;
  var h = 0;

  function resize() {
    var dpr = window.devicePixelRatio || 1;
    w = canvas.clientWidth || 800;
    h = canvas.clientHeight || 600;
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
  }

  /* ---------- Geometry helpers ---------- */

  function circlePoints(center, axis, radius, segments) {
    var a = norm(axis);
    var ref = Math.abs(a[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    var u = norm(cross(a, ref));
    var v = cross(a, u);
    var points = [];
    for (var i = 0; i <= segments; i += 1) {
      var t = (2 * Math.PI * i) / segments;
      var ct = Math.cos(t) * radius;
      var st = Math.sin(t) * radius;
      points.push([
        center[0] + u[0] * ct + v[0] * st,
        center[1] + u[1] * ct + v[1] * st,
        center[2] + u[2] * ct + v[2] * st,
      ]);
    }
    return points;
  }

  function strokePath(points) {
    ctx.beginPath();
    for (var i = 0; i < points.length; i += 1) {
      var p = camera(points[i]);
      if (i === 0) {
        ctx.moveTo(p[0], p[1]);
      } else {
        ctx.lineTo(p[0], p[1]);
      }
    }
    ctx.stroke();
  }

  /* ---------- Involute gear profile (port of scripts/gear_profile.py) ---------- */

  var PRESSURE_ANGLE_RAD = (20.0 * Math.PI) / 180.0;
  var ADDENDUM_COEFF = 1.0;
  var DEDENDUM_COEFF = 1.25;
  var C_C_BOND_M = 1.544e-10;

  function involuteFn(a) {
    return Math.tan(a) - a;
  }

  function pitchRadiusJS(moduleM, teeth) {
    return (moduleM * teeth) / 2;
  }

  function baseRadiusJS(moduleM, teeth) {
    return pitchRadiusJS(moduleM, teeth) * Math.cos(PRESSURE_ANGLE_RAD);
  }

  function outerRadiusJS(moduleM, teeth) {
    return pitchRadiusJS(moduleM, teeth) + ADDENDUM_COEFF * moduleM;
  }

  function rootRadiusJS(moduleM, teeth) {
    return pitchRadiusJS(moduleM, teeth) - DEDENDUM_COEFF * moduleM;
  }

  function flankAngleRad(radiusM, sign, moduleM, teeth) {
    var base = baseRadiusJS(moduleM, teeth);
    if (radiusM < base) {
      throw new Error("radius inside the base circle");
    }
    var cosine = Math.max(-1, Math.min(1, base / radiusM));
    var alphaR = Math.acos(cosine);
    var halfTooth = Math.PI / (2 * teeth);
    return (
      sign * (halfTooth + involuteFn(PRESSURE_ANGLE_RAD) - involuteFn(alphaR))
    );
  }

  function outlinePointsJS(moduleM, teeth, flankSamples, arcSamples) {
    var flanks = Math.max(1, flankSamples);
    var arcs = Math.max(1, arcSamples);
    var toothPitch = (2 * Math.PI) / teeth;
    var halfPitch = Math.PI / teeth;
    var outerM = outerRadiusJS(moduleM, teeth);
    var rootM = rootRadiusJS(moduleM, teeth);
    var baseM = baseRadiusJS(moduleM, teeth);
    var tipHalf = halfPitch + involuteFn(PRESSURE_ANGLE_RAD)
      - involuteFn(Math.acos(Math.max(-1, Math.min(1, baseM / outerM))));
    var flankStart = Math.max(rootM, baseM);
    var flankStartAngle = flankAngleRad(flankStart, -1, moduleM, teeth);
    var local = [];
    var i;
    var t;
    if (rootM < baseM) {
      for (i = 0; i <= arcs; i += 1) {
        t = i / arcs;
        var a0 = -halfPitch + (halfPitch + flankStartAngle) * t;
        local.push([rootM * Math.cos(a0), rootM * Math.sin(a0)]);
      }
    }
    local.push([
      flankStart * Math.cos(flankStartAngle),
      flankStart * Math.sin(flankStartAngle),
    ]);
    for (i = 1; i <= flanks; i += 1) {
      t = i / flanks;
      var r1 = flankStart + (outerM - flankStart) * t;
      var a1 = flankAngleRad(r1, -1, moduleM, teeth);
      local.push([r1 * Math.cos(a1), r1 * Math.sin(a1)]);
    }
    for (i = 0; i <= arcs; i += 1) {
      t = i / arcs;
      var a2 = -tipHalf + 2 * tipHalf * t;
      local.push([outerM * Math.cos(a2), outerM * Math.sin(a2)]);
    }
    for (i = flanks - 1; i >= 0; i -= 1) {
      t = i / flanks;
      var r2 = flankStart + (outerM - flankStart) * t;
      var a3 = flankAngleRad(r2, 1, moduleM, teeth);
      local.push([r2 * Math.cos(a3), r2 * Math.sin(a3)]);
    }
    if (rootM < baseM) {
      for (i = 0; i <= arcs; i += 1) {
        t = i / arcs;
        var a4 = -flankStartAngle + (halfPitch + flankStartAngle) * t;
        local.push([rootM * Math.cos(a4), rootM * Math.sin(a4)]);
      }
    }
    var points = [];
    for (var tooth = 0; tooth < teeth; tooth += 1) {
      var rotation = tooth * toothPitch;
      var cosR = Math.cos(rotation);
      var sinR = Math.sin(rotation);
      for (var k = 0; k < local.length; k += 1) {
        points.push([
          local[k][0] * cosR - local[k][1] * sinR,
          local[k][0] * sinR + local[k][1] * cosR,
        ]);
      }
    }
    return points;
  }

  function reflectToInternalJS(points, pitchM) {
    return points.map(function (p) {
      var radius = Math.sqrt(p[0] * p[0] + p[1] * p[1]);
      if (radius <= 0) {
        return [p[0], p[1]];
      }
      var scale = (2 * pitchM - radius) / radius;
      return [p[0] * scale, p[1] * scale];
    });
  }

  function gearOutlineWorldPoints(role, index, design) {
    var moduleM = design.module_m;
    var teeth;
    var internal = false;
    var rotation = 0;
    if (role === "sun") {
      teeth = design.sun_teeth;
    } else if (role === "planet") {
      teeth = design.planet_teeth;
      rotation =
        (index * 2 * Math.PI) / design.planet_count + Math.PI / design.planet_teeth;
    } else if (role === "ring") {
      teeth = design.ring_teeth;
      internal = true;
    } else {
      return null;
    }
    var points = outlinePointsJS(moduleM, teeth, 4, 2);
    if (internal) {
      points = reflectToInternalJS(points, pitchRadiusJS(moduleM, teeth));
    }
    if (rotation === 0) {
      return points;
    }
    var cosR = Math.cos(rotation);
    var sinR = Math.sin(rotation);
    return points.map(function (p) {
      return [p[0] * cosR - p[1] * sinR, p[0] * sinR + p[1] * cosR];
    });
  }

  /* ---------- Layer drawing ---------- */

  /* ---------- Level of detail ---------- */

  function atomSpacingPx() {
    return (C_C_BOND_M / Math.max(fit.r, 1e-30)) * pixelSize();
  }

  function shadeColor(hex, factor) {
    var n = parseInt(hex.slice(1), 16);
    var r = Math.min(255, Math.round(((n >> 16) & 255) * factor));
    var g = Math.min(255, Math.round(((n >> 8) & 255) * factor));
    var b = Math.min(255, Math.round((n & 255) * factor));
    return "rgb(" + r + "," + g + "," + b + ")";
  }

  function drawGearSchematic(alpha) {
    if (!scene.design || !scene.design.module_m) {
      return;
    }
    var planetIndex = 0;
    ctx.globalAlpha = alpha;
    ctx.lineWidth = 1.4;
    scene.device.bodies.forEach(function (body) {
      var role = body.role;
      var index = 0;
      if (role === "planet") {
        index = planetIndex;
        planetIndex += 1;
      }
      var local = gearOutlineWorldPoints(role, index, scene.design);
      if (!local) {
        return;
      }
      ctx.strokeStyle = ROLE_COLORS[role] || "#8b98a5";
      ctx.beginPath();
      for (var i = 0; i < local.length; i += 1) {
        var p = camera([
          body.position_m[0] + local[i][0],
          body.position_m[1] + local[i][1],
          body.position_m[2],
        ]);
        if (i === 0) {
          ctx.moveTo(p[0], p[1]);
        } else {
          ctx.lineTo(p[0], p[1]);
        }
      }
      ctx.closePath();
      ctx.stroke();
    });
    ctx.globalAlpha = 1;
  }

  function drawAtoms(alpha) {
    var atoms = scene.atomistic.atoms;
    var spacing = atomSpacingPx();
    var radius = Math.max(0.6, spacing * 0.30);
    var groups = new Map();
    for (var i = 0; i < atoms.length; i += 1) {
      var p = camera(atoms[i].position_m);
      var bin = Math.max(0, Math.min(3, Math.floor((p[2] + 1) * 2)));
      var key = atoms[i].element + ":" + bin;
      if (!groups.has(key)) {
        groups.set(key, { element: atoms[i].element, bin: bin, points: [] });
      }
      groups.get(key).points.push(p);
    }
    ctx.globalAlpha = alpha;
    groups.forEach(function (group) {
      var base = elementColor(group.element);
      var factor = 0.45 + 0.18 * group.bin;
      ctx.fillStyle = base.charAt(0) === "#" ? shadeColor(base, factor) : base;
      ctx.beginPath();
      for (var j = 0; j < group.points.length; j += 1) {
        var q = group.points[j];
        ctx.moveTo(q[0] + radius, q[1]);
        ctx.arc(q[0], q[1], radius, 0, 2 * Math.PI);
      }
      ctx.fill();
    });
    ctx.globalAlpha = 1;
  }

  function drawAtomistic(alpha) {
    if (atomSpacingPx() < 6) {
      drawGearSchematic(alpha);
      return;
    }
    drawAtoms(alpha);
  }

  function drawCoarse(alpha) {
    ctx.globalAlpha = alpha;
    ctx.lineWidth = 1.6;
    for (var i = 0; i < scene.coarse.bodies.length; i += 1) {
      var body = scene.coarse.bodies[i];
      var color = ROLE_COLORS[body.role] || "#8b949e";
      ctx.strokeStyle = color;

      if (body.pitch_circle) {
        ctx.setLineDash([6, 4]);
        strokePath(
          circlePoints(
            body.pitch_circle.center_m,
            body.pitch_circle.axis,
            body.pitch_circle.radius_m,
            72
          )
        );
      }
      if (body.bounding_cylinder) {
        ctx.setLineDash([]);
        var c = body.bounding_cylinder;
        var half = c.half_length_m;
        var low = [c.center_m[0], c.center_m[1], c.center_m[2] - half];
        var high = [c.center_m[0], c.center_m[1], c.center_m[2] + half];
        strokePath(circlePoints(low, c.axis, c.radius_m, 48));
        strokePath(circlePoints(high, c.axis, c.radius_m, 48));
        var a = camera(low);
        var b = camera(high);
        ctx.beginPath();
        ctx.moveTo(a[0], a[1]);
        ctx.lineTo(b[0], b[1]);
        ctx.stroke();
      }
    }
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
  }

  function drawDevice(alpha) {
    var bodies = scene.device.bodies;
    var coarseByIndex = new Map();
    for (var k = 0; k < scene.coarse.bodies.length; k += 1) {
      coarseByIndex.set(scene.coarse.bodies[k].index, scene.coarse.bodies[k]);
    }

    ctx.globalAlpha = alpha;

    ctx.strokeStyle = "rgba(230, 237, 243, 0.35)";
    ctx.lineWidth = 1.2;
    for (var j = 0; j < scene.device.joints.length; j += 1) {
      var joint = scene.device.joints[j];
      var a = camera(bodies[joint.body_a].position_m);
      var b = camera(bodies[joint.body_b].position_m);
      ctx.beginPath();
      ctx.moveTo(a[0], a[1]);
      ctx.lineTo(b[0], b[1]);
      ctx.stroke();

      if (joint.anchor_world_m) {
        var anchor = camera(joint.anchor_world_m);
        ctx.fillStyle = "#e6edf3";
        ctx.beginPath();
        ctx.arc(anchor[0], anchor[1], 2.5, 0, 2 * Math.PI);
        ctx.fill();
      }
    }

    for (var i = 0; i < bodies.length; i += 1) {
      var body = bodies[i];
      var color = ROLE_COLORS[body.role] || "#8b949e";
      var p = camera(body.position_m);
      var coarse = coarseByIndex.get(body.index);
      var radius = 4;
      if (coarse && coarse.bounding_cylinder) {
        radius = Math.max(3, pixRadius(coarse.bounding_cylinder.radius_m, p[3]));
      }
      ctx.fillStyle = color;
      ctx.globalAlpha = alpha * 0.35;
      ctx.beginPath();
      ctx.arc(p[0], p[1], radius, 0, 2 * Math.PI);
      ctx.fill();
      ctx.globalAlpha = alpha;
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(p[0], p[1], radius, 0, 2 * Math.PI);
      ctx.stroke();

      if (body.fixed) {
        ctx.fillStyle = "#e6edf3";
        ctx.beginPath();
        ctx.arc(p[0], p[1], 2, 0, 2 * Math.PI);
        ctx.fill();
      }
      ctx.fillStyle = "rgba(230, 237, 243, 0.85)";
      ctx.font = "11px ui-monospace, Menlo, monospace";
      ctx.fillText(body.name, p[0] + radius + 3, p[1] - radius - 2);
    }
    ctx.globalAlpha = 1;
  }

  /* ---------- Layer weights and frame ---------- */

  function layerWeights() {
    var v = parseFloat(slider.value);
    var raw = [
      Math.max(0, 1 - Math.abs(v - 0)),
      Math.max(0, 1 - Math.abs(v - 1)),
      Math.max(0, 1 - Math.abs(v - 2)),
    ];
    return [
      toggles.atomistic.checked ? raw[0] : 0,
      toggles.device.checked ? raw[1] : 0,
      toggles.coarse.checked ? raw[2] : 0,
    ];
  }

  var LAYER_NAMES = ["Atomistic (L1)", "Device (L2)", "Coarse (handoff)"];

  function updateSliderLabel() {
    var v = parseFloat(slider.value);
    var nearest = Math.round(v);
    sliderLabel.textContent =
      LAYER_NAMES[nearest] + (Math.abs(v - nearest) > 0.01 ? " (blend)" : "");
  }

  function draw() {
    ctx.clearRect(0, 0, w, h);
    if (!scene) {
      return;
    }
    var weights = layerWeights();
    if (weights[2] > 0) {
      drawCoarse(weights[2]);
    }
    if (weights[1] > 0) {
      drawDevice(weights[1]);
    }
    if (weights[0] > 0) {
      drawAtomistic(weights[0]);
    }
    drawCaption();
  }

  function drawCaption() {
    ctx.globalAlpha = 1;
    ctx.fillStyle = "rgba(240, 180, 41, 0.9)";
    ctx.font = "bold 12px ui-monospace, Menlo, monospace";
    ctx.fillText("simulated / schematic", 12, h - 14);
    ctx.fillStyle = "rgba(139, 148, 158, 0.9)";
    ctx.font = "11px ui-monospace, Menlo, monospace";
    ctx.fillText(
      "yaw " + view.yaw.toFixed(2) + " rad  tilt " + view.tilt.toFixed(2) + " rad",
      12,
      h - 30
    );
  }

  /* ---------- Readout ---------- */

  function setReadout() {
    var d = scene.design;
    var lines = [
      "schema: " + scene.schema + " v" + scene.version,
      "sun " + d.sun_teeth + "t  planet " + d.planet_teeth + "t  ring " + d.ring_teeth + "t",
      "planets: " + d.planet_count,
      "gear ratio (ring fixed): " + d.gear_ratio,
      "bodies: " + scene.device.body_count,
      "joints: " + scene.device.joints.length,
      "atoms: " + scene.atomistic.atom_count,
      "module: " + d.module_m.toExponential(3) + " m",
    ];
    readout.textContent = lines.join("\n");
  }

  /* ---------- Loading ---------- */

  function setScene(data) {
    if (!data || data.schema !== "nanocad.scene") {
      throw new Error("not a nanocad.scene document");
    }
    scene = data;
    if (scene.design && scene.design.layers !== undefined) {
      currentParams = paramsFromDesign(scene.design);
      syncParamInputs();
    }
    message.textContent = "";
    fallback.hidden = true;
    computeFit();
    setReadout();
    resize();
  }

  function showFallback(reason) {
    message.textContent = "scene.json not loaded: " + reason;
    fallback.hidden = false;
  }

  /*
   * The companion `scene.data.js` sets `window.NANOCAD_SCENE`. A `<script>`
   * tag loads a local file from `file://`, but `fetch` does not. Try the
   * companion first, then fall back to fetch and to the file picker.
   */
  function loadFromCompanion() {
    return new Promise(function (resolve) {
      if (window.NANOCAD_SCENE) {
        resolve(true);
        return;
      }
      var script = document.createElement("script");
      script.src = "scene.data.js";
      script.onload = function () {
        resolve(Boolean(window.NANOCAD_SCENE));
      };
      script.onerror = function () {
        resolve(false);
      };
      document.head.appendChild(script);
    });
  }

  function load() {
    loadFromCompanion()
      .then(function (fromCompanion) {
        if (fromCompanion) {
          setScene(window.NANOCAD_SCENE);
          return null;
        }
        if (!window.fetch) {
          throw new Error("fetch is unavailable");
        }
        return fetch("scene.json", { cache: "no-store" }).then(function (response) {
          if (!response.ok) {
            throw new Error("HTTP " + response.status);
          }
          return response.json().then(setScene);
        });
      })
      .catch(function (error) {
        showFallback(error.message);
      });
  }

  fileInput.addEventListener("change", function () {
    var file = fileInput.files && fileInput.files[0];
    if (!file) {
      return;
    }
    var reader = new FileReader();
    reader.onload = function () {
      try {
        setScene(JSON.parse(reader.result));
      } catch (error) {
        showFallback(error.message);
      }
    };
    reader.readAsText(file);
  });

  window.addEventListener("dragover", function (event) {
    event.preventDefault();
  });
  window.addEventListener("drop", function (event) {
    event.preventDefault();
    var file = event.dataTransfer && event.dataTransfer.files[0];
    if (!file) {
      return;
    }
    var reader = new FileReader();
    reader.onload = function () {
      try {
        setScene(JSON.parse(reader.result));
      } catch (error) {
        showFallback(error.message);
      }
    };
    reader.readAsText(file);
  });

  /* ---------- Interaction ---------- */

  canvas.addEventListener("pointerdown", function (event) {
    drag = { x: event.clientX, y: event.clientY };
    canvas.setPointerCapture(event.pointerId);
  });

  canvas.addEventListener("pointermove", function (event) {
    if (!drag) {
      return;
    }
    view.yaw += (event.clientX - drag.x) * 0.01;
    view.tilt += (event.clientY - drag.y) * 0.01;
    view.tilt = Math.max(0.05, Math.min(Math.PI - 0.05, view.tilt));
    drag = { x: event.clientX, y: event.clientY };
    draw();
  });

  canvas.addEventListener("pointerup", function (event) {
    drag = null;
    if (canvas.hasPointerCapture(event.pointerId)) {
      canvas.releasePointerCapture(event.pointerId);
    }
  });

  canvas.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      view.zoom *= event.deltaY < 0 ? 1.1 : 0.9;
      view.zoom = Math.max(0.2, Math.min(8, view.zoom));
      draw();
    },
    { passive: false }
  );

  canvas.addEventListener("dblclick", function () {
    view.yaw = -0.65;
    view.tilt = 1.02;
    view.zoom = 1;
    draw();
  });

  slider.addEventListener("input", function () {
    updateSliderLabel();
    draw();
  });

  Object.keys(toggles).forEach(function (name) {
    toggles[name].addEventListener("change", draw);
  });

  window.addEventListener("resize", resize);

  /* ---------- Live parameters and chat (served by app/server.py) ---------- */

  var paramsEl = document.getElementById("params");
  var paramReset = document.getElementById("param-reset");
  var chatLog = document.getElementById("chat-log");
  var chatForm = document.getElementById("chat-form");
  var chatInput = document.getElementById("chat-input");
  var chatHint = document.getElementById("chat-hint");
  var meta = null;
  var live = false;
  var currentParams = null;
  var busy = false;

  function paramsFromDesign(design) {
    return {
      module_m: design.module_m,
      sun_teeth: design.sun_teeth,
      planet_teeth: design.planet_teeth,
      planet_count: design.planet_count,
      layers: design.layers,
    };
  }

  function setMessage(text) {
    message.textContent = text;
  }

  function applyScene(payload) {
    if (payload && payload.scene) {
      setScene(payload.scene);
      currentParams = paramsFromDesign(payload.scene.design);
      syncParamInputs();
      draw();
      setMessage("");
    }
  }

  function syncParamInputs() {
    if (!meta || !currentParams) {
      return;
    }
    meta.parameters.forEach(function (entry) {
      var input = document.getElementById("param-" + entry.key);
      if (!input) {
        return;
      }
      var shown = currentParams[entry.key] / entry.scale;
      input.value = String(Number(shown.toPrecision(6)));
    });
  }

  function buildParamInputs() {
    paramsEl.innerHTML = "";
    meta.parameters.forEach(function (entry) {
      var row = document.createElement("label");
      row.className = "param-row";
      var span = document.createElement("span");
      span.textContent = entry.label + (entry.unit ? " (" + entry.unit + ")" : "");
      var input = document.createElement("input");
      input.type = "number";
      input.id = "param-" + entry.key;
      var limits = meta.limits[entry.key];
      input.min = String(limits[0] / entry.scale);
      input.max = String(limits[1] / entry.scale);
      input.step = entry.scale === 1 ? "1" : "any";
      input.addEventListener("change", buildFromInputs);
      row.appendChild(span);
      row.appendChild(input);
      paramsEl.appendChild(row);
    });
    syncParamInputs();
  }

  function buildFromInputs() {
    if (!meta || !live) {
      return;
    }
    var params = {};
    meta.parameters.forEach(function (entry) {
      var input = document.getElementById("param-" + entry.key);
      if (!input) {
        return;
      }
      var value = parseFloat(input.value);
      if (isFinite(value)) {
        params[entry.key] = value * entry.scale;
      }
    });
    requestBuild(params);
  }

  function requestBuild(params) {
    if (busy || !live) {
      return;
    }
    busy = true;
    setMessage("building ...");
    var query = Object.keys(params)
      .map(function (key) {
        return encodeURIComponent(key) + "=" + encodeURIComponent(params[key]);
      })
      .join("&");
    fetch("api/build?" + query)
      .then(function (response) {
        return response.json();
      })
      .then(function (payload) {
        if (payload.error) {
          setMessage("rejected by the engine: " + payload.error);
          return;
        }
        applyScene(payload);
      })
      .catch(function (error) {
        setMessage("build failed: " + error.message);
      })
      .then(function () {
        busy = false;
      });
  }

  function addChat(who, text) {
    var line = document.createElement("div");
    line.className = "chat-" + who;
    line.textContent = text;
    chatLog.appendChild(line);
    chatLog.scrollTop = chatLog.scrollHeight;
  }

  function sendChat(event) {
    event.preventDefault();
    var text = chatInput.value.trim();
    if (!text) {
      return;
    }
    if (!live) {
      addChat("bot", "Start python3 app/server.py to edit the design.");
      return;
    }
    addChat("you", text);
    chatInput.value = "";
    var body = {
      message: text,
      params: currentParams || paramsFromDesign(scene.design),
    };
    fetch("api/chat", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
    })
      .then(function (response) {
        return response.json();
      })
      .then(function (payload) {
        addChat("bot", payload.reply);
        currentParams = payload.params;
        if (payload.scene) {
          applyScene(payload);
        } else {
          syncParamInputs();
        }
      })
      .catch(function (error) {
        addChat("bot", "error: " + error.message);
      });
  }

  function initApp() {
    if (!window.fetch || window.location.protocol === "file:") {
      paramsEl.textContent = "Start python3 app/server.py, then reload this page.";
      chatHint.textContent =
        "Live edit needs the local server: python3 app/server.py";
      return;
    }
    chatForm.addEventListener("submit", sendChat);
    paramReset.addEventListener("click", function () {
      if (meta) {
        requestBuild(meta.defaults);
      }
    });
    fetch("api/meta")
      .then(function (response) {
        return response.json();
      })
      .then(function (data) {
        meta = data;
        live = true;
        chatHint.textContent =
          'Try "one atomic layer thicker", "6 atoms thick", ' +
          '"add more teeth", "4 planets", or "reset".';
        buildParamInputs();
        if (!currentParams && scene) {
          currentParams = paramsFromDesign(scene.design);
        }
        syncParamInputs();
      })
      .catch(function () {
        chatHint.textContent =
          "Live edit needs the local server: python3 app/server.py";
      });
  }

  updateSliderLabel();
  load();
  initApp();
})();
