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
  var readout = document.getElementById("readout");
  var message = document.getElementById("message");
  var fallback = document.getElementById("load-fallback");
  var fileInput = document.getElementById("file-input");

  var pickEl = document.getElementById("pick");
  var measureEl = document.getElementById("measure");
  var legendEl = document.getElementById("legend");
  var aoToggle = document.getElementById("toggle-ao");
  var atomSizeInput = document.getElementById("atom-size");
  var clipInput = document.getElementById("clip");
  var speedInput = document.getElementById("speed");
  var scaleBarLine = document.getElementById("scale-bar-line");
  var scaleBarLabel = document.getElementById("scale-bar-label");
  var triad = {
    x: document.getElementById("triad-x"),
    y: document.getElementById("triad-y"),
    z: document.getElementById("triad-z"),
  };
  var buttons = {
    reset: document.getElementById("view-reset"),
    iso: document.getElementById("view-iso"),
    top: document.getElementById("view-top"),
    front: document.getElementById("view-front"),
    play: document.getElementById("view-play"),
    turntable: document.getElementById("view-turntable"),
    full: document.getElementById("view-full"),
    save: document.getElementById("view-save"),
  };

  var toggles = {
    atomistic: document.getElementById("toggle-atomistic"),
    device: document.getElementById("toggle-device"),
    coarse: document.getElementById("toggle-coarse"),
  };

  var panelEl = document.getElementById("panel");
  var panelToggleEl = document.getElementById("panel-toggle");
  if (panelEl && panelToggleEl) {
    panelToggleEl.addEventListener("click", function () {
      var collapsed = panelEl.classList.toggle("collapsed");
      panelToggleEl.textContent = collapsed ? "Show" : "Hide";
    });
  }

  var glCanvas = document.getElementById("gl-canvas");
  var renderer = null;
  try {
    if (window.NanoCadAtoms && glCanvas) {
      renderer = window.NanoCadAtoms.createAtomRenderer(glCanvas);
    }
  } catch (error) {
    renderer = null;
  }

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

  var view = { yaw: -0.65, tilt: 1.02, zoom: 1, panX: 0, panY: 0 };
  var fit = { c: [0, 0, 0], r: 1 };
  var scene = null;

  var display = {
    ao: true,
    atomSize: 1,
    clip: 1,
    hiddenElements: {},
  };
  var atomCache = null;
  var motion = { playing: false, turntable: false, time: 0, speed: 1, last: 0 };
  var hoverAtom = -1;
  var measure = [];
  var livePositions = null;

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

    return [w / 2 + x1 * scale * persp + view.panX, h / 2 + y2 * scale * persp + view.panY, depth, persp];
  }

  function pixRadius(radius_m, persp) {
    return (radius_m / fit.r) * pixelSize() * persp;
  }

  /* ---------- WebGL atom pass ---------- */

  var ATOM_RADII = {
    H: 0.31e-10,
    C: 0.76e-10,
    N: 0.71e-10,
    O: 0.66e-10,
    F: 0.57e-10,
    Si: 1.11e-10,
    P: 1.07e-10,
    S: 1.05e-10,
    B: 0.84e-10,
  };

  function atomRadius(element) {
    return ATOM_RADII[element] || 0.7e-10;
  }

  function hexToRgb(hex) {
    var n = parseInt(hex.slice(1), 16);
    return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
  }

  function glState() {
    var zmin = fitMin[2];
    var zmax = fitMax[2];
    var clipZ = zmin + display.clip * (zmax - zmin);
    return {
      center: fit.c,
      invR: 1 / Math.max(fit.r, 1e-30),
      yaw: view.yaw,
      tilt: view.tilt,
      camDistance: 3.2,
      perspK: 0.6,
      scale: pixelSize(),
      viewport: [w, h],
      pan: [view.panX, view.panY],
      near: 2.4,
      far: 4.0,
      radiusScale: display.atomSize,
      ballStick: display.atomSize < 0.72,
      keyDir: [-0.42, 0.55, 0.72],
      clipZ: clipZ,
      clipOn: display.clip < 0.999,
      showBonds: false,
      ao: display.ao,
      aoRadius: 0.02,
      aoStrength: 1.0,
    };
  }

  function buildAtomBuffers() {
    var atoms = scene.atomistic.atoms;
    var position = new Float32Array(atoms.length * 3);
    var radii = new Float32Array(atoms.length);
    var colors = new Float32Array(atoms.length * 3);
    for (var i = 0; i < atoms.length; i += 1) {
      position[i * 3] = atoms[i].position_m[0];
      position[i * 3 + 1] = atoms[i].position_m[1];
      position[i * 3 + 2] = atoms[i].position_m[2];
      radii[i] = atomRadius(atoms[i].element);
      var rgb = hexToRgb(elementColor(atoms[i].element));
      colors[i * 3] = rgb[0];
      colors[i * 3 + 1] = rgb[1];
      colors[i * 3 + 2] = rgb[2];
    }
    var design = scene.design || {};
    var ns = design.sun_teeth || 24;
    var nr = design.ring_teeth || 60;
    motion.w_sun = -0.34;
    motion.w_carrier = (motion.w_sun * ns) / (ns + nr);
    motion.w_planet = (-2 / 3) * motion.w_sun;
    atomCache = {
      atoms: atoms,
      base: position,
      position: position.slice(),
      radii: radii,
      colors: colors,
    };
    uploadAtoms();
  }

  function hasHiddenElements() {
    for (var key in display.hiddenElements) {
      if (display.hiddenElements[key]) {
        return true;
      }
    }
    return false;
  }

  function uploadAtoms() {
    if (!renderer || !atomCache) {
      return;
    }
    var cache = atomCache;
    if (!hasHiddenElements()) {
      renderer.setAtoms(cache.position, cache.radii, cache.colors);
      livePositions = cache.position;
      return;
    }
    var kept = [];
    for (var i = 0; i < cache.atoms.length; i += 1) {
      if (!display.hiddenElements[cache.atoms[i].element]) {
        kept.push(i);
      }
    }
    var position = new Float32Array(kept.length * 3);
    var radii = new Float32Array(kept.length);
    var colors = new Float32Array(kept.length * 3);
    for (var k = 0; k < kept.length; k += 1) {
      var j = kept[k];
      position[k * 3] = cache.position[j * 3];
      position[k * 3 + 1] = cache.position[j * 3 + 1];
      position[k * 3 + 2] = cache.position[j * 3 + 2];
      radii[k] = cache.radii[j];
      colors[k * 3] = cache.colors[j * 3];
      colors[k * 3 + 1] = cache.colors[j * 3 + 1];
      colors[k * 3 + 2] = cache.colors[j * 3 + 2];
    }
    renderer.setAtoms(position, radii, colors);
    livePositions = position;
  }

  /* Advance the planetary kinematics. The ring is fixed. */
  function updateAnimation(dt) {
    if (!atomCache) {
      return;
    }
    motion.time += dt;
    var t = motion.time;
    var ca = Math.cos(motion.w_carrier * t);
    var sa = Math.sin(motion.w_carrier * t);
    var cb = Math.cos(motion.w_planet * t);
    var sb = Math.sin(motion.w_planet * t);
    var cc = Math.cos(motion.w_sun * t);
    var sc = Math.sin(motion.w_sun * t);
    var base = atomCache.base;
    var live = atomCache.position;
    var atoms = atomCache.atoms;
    var bodies = scene.device.bodies;
    for (var i = 0; i < atoms.length; i += 1) {
      var body = atoms[i].body;
      var x = base[i * 3];
      var y = base[i * 3 + 1];
      var ox;
      var oy;
      if (body === 1) {
        ox = x * cc - y * sc;
        oy = x * sc + y * cc;
      } else if (body >= 2 && body <= 4) {
        var center = bodies[body].position_m;
        var dx = x - center[0];
        var dy = y - center[1];
        ox = center[0] * ca - center[1] * sa + (dx * cb - dy * sb);
        oy = center[0] * sa + center[1] * ca + (dx * sb + dy * cb);
      } else if (body === 6) {
        ox = x * ca - y * sa;
        oy = x * sa + y * ca;
      } else {
        ox = x;
        oy = y;
      }
      live[i * 3] = ox;
      live[i * 3 + 1] = oy;
    }
  }

  /* ---------- Legend ---------- */

  function buildLegend() {
    legendEl.textContent = "";
    var seen = [];
    for (var i = 0; i < scene.atomistic.atoms.length; i += 1) {
      var el = scene.atomistic.atoms[i].element;
      if (seen.indexOf(el) < 0) {
        seen.push(el);
      }
    }
    seen.forEach(function (element) {
      var row = document.createElement("label");
      row.className = "legend-row";
      var box = document.createElement("input");
      box.type = "checkbox";
      box.setAttribute("data-element", element);
      box.checked = !display.hiddenElements[element];
      box.addEventListener("change", function () {
        display.hiddenElements[element] = !box.checked;
        uploadAtoms();
        draw();
      });
      var swatch = document.createElement("span");
      swatch.className = "swatch";
      swatch.style.background = elementColor(element);
      var text = document.createElement("span");
      text.textContent = element;
      row.appendChild(box);
      row.appendChild(swatch);
      row.appendChild(text);
      legendEl.appendChild(row);
    });
  }

  /* ---------- Overlays ---------- */

  function updateTriad() {
    var cy = Math.cos(view.yaw);
    var sy = Math.sin(view.yaw);
    var ct = Math.cos(view.tilt);
    var st = Math.sin(view.tilt);
    var axes = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] };
    Object.keys(axes).forEach(function (key) {
      var n = axes[key];
      var x1 = n[0] * cy - n[1] * sy;
      var y1 = n[0] * sy + n[1] * cy;
      var y2 = y1 * ct - n[2] * st;
      var line = triad[key];
      line.setAttribute("x2", (32 + x1 * 18).toFixed(1));
      line.setAttribute("y2", (32 - y2 * 18).toFixed(1));
    });
  }

  function updateScaleBar() {
    var scale = pixelSize();
    if (!(scale > 0) || !(fit.r > 0)) {
      return;
    }
    var metersPerPixel = fit.r / scale;
    var nm = (metersPerPixel * 110) / 1e-9;
    var pow = Math.pow(10, Math.floor(Math.log10(nm)));
    var mult = nm / pow;
    var nice = mult < 1.5 ? 1 : mult < 3.5 ? 2 : mult < 7.5 ? 5 : 10;
    var nmNice = nice * pow;
    var pixels = (nmNice * 1e-9) / metersPerPixel;
    scaleBarLine.style.width = Math.max(20, Math.round(pixels)) + "px";
    scaleBarLabel.textContent =
      (nmNice >= 1 ? nmNice.toFixed(0) : nmNice.toFixed(2)) + " nm";
  }

  function atomPoint(index) {
    var p = livePositions || (atomCache && atomCache.position);
    if (!p) {
      return [0, 0, 0];
    }
    return [p[index * 3], p[index * 3 + 1], p[index * 3 + 2]];
  }

  function pickAtom(mx, my) {
    if (!atomCache) {
      return -1;
    }
    var cy = Math.cos(view.yaw);
    var sy = Math.sin(view.yaw);
    var ct = Math.cos(view.tilt);
    var st = Math.sin(view.tilt);
    var scale = pixelSize();
    var halfW = w / 2;
    var halfH = h / 2;
    var invR = 1 / Math.max(fit.r, 1e-30);
    var camDistance = 3.2;
    var c0 = fit.c[0];
    var c1 = fit.c[1];
    var c2 = fit.c[2];
    var p = livePositions || atomCache.position;
    var best = -1;
    var bestDistance = 14 * 14;
    for (var i = 0; i < atomCache.atoms.length; i += 1) {
      var nx = (p[i * 3] - c0) * invR;
      var ny = (p[i * 3 + 1] - c1) * invR;
      var nz = (p[i * 3 + 2] - c2) * invR;
      var x1 = nx * cy - ny * sy;
      var y1 = nx * sy + ny * cy;
      var y2 = y1 * ct - nz * st;
      var depth = y1 * st + nz * ct;
      var persp = camDistance / (camDistance - depth * 0.6);
      var sx = halfW + x1 * scale * persp + view.panX;
      var sy2 = halfH + y2 * scale * persp + view.panY;
      var dx = sx - mx;
      var dy = sy2 - my;
      var d = dx * dx + dy * dy;
      if (d < bestDistance) {
        bestDistance = d;
        best = i;
      }
    }
    return best;
  }

  function describeAtom(index) {
    var atom = atomCache.atoms[index];
    var p = atomPoint(index);
    var role = (scene.device.bodies[atom.body] || {}).name || "?";
    return (
      atom.element +
      " #" +
      index +
      "  " +
      role +
      "  (" +
      (p[0] / 1e-9).toFixed(3) +
      ", " +
      (p[1] / 1e-9).toFixed(3) +
      ", " +
      (p[2] / 1e-9).toFixed(3) +
      ") nm"
    );
  }

  function drawOverlay() {
    updateTriad();
    updateScaleBar();
    if (!atomCache) {
      return;
    }
    if (hoverAtom >= 0) {
      var hp = camera(atomPoint(hoverAtom));
      ctx.globalAlpha = 1;
      ctx.strokeStyle = "#f2c14e";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(hp[0], hp[1], 9, 0, 2 * Math.PI);
      ctx.stroke();
    }
    if (measure.length === 2) {
      var a = camera(atomPoint(measure[0]));
      var b = camera(atomPoint(measure[1]));
      ctx.strokeStyle = "#f2c14e";
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      ctx.moveTo(a[0], a[1]);
      ctx.lineTo(b[0], b[1]);
      ctx.stroke();
      ctx.setLineDash([]);
    }
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
    if (renderer) {
      glCanvas.width = Math.round(w * dpr);
      glCanvas.height = Math.round(h * dpr);
      renderer.resize(w * dpr, h * dpr);
    }
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
    var t = motion.time;
    var ca = Math.cos(motion.w_carrier * t);
    var sa = Math.sin(motion.w_carrier * t);
    var cb = Math.cos(motion.w_planet * t);
    var sb = Math.sin(motion.w_planet * t);
    var cc = Math.cos(motion.w_sun * t);
    var sc = Math.sin(motion.w_sun * t);
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
      var cx = body.position_m[0];
      var cy = body.position_m[1];
      var rcx = cx * ca - cy * sa;
      var rcy = cx * sa + cy * ca;
      ctx.strokeStyle = ROLE_COLORS[role] || "#8b98a5";
      ctx.beginPath();
      for (var i = 0; i < local.length; i += 1) {
        var lx = local[i][0];
        var ly = local[i][1];
        var x;
        var y;
        if (role === "sun") {
          x = cx + lx;
          y = cy + ly;
          var sx = x * cc - y * sc;
          var sy = x * sc + y * cc;
          x = sx;
          y = sy;
        } else if (role === "planet") {
          x = rcx + (lx * cb - ly * sb);
          y = rcy + (lx * sb + ly * cb);
        } else {
          x = cx + lx;
          y = cy + ly;
        }
        var p = camera([x, y, body.position_m[2]]);
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
    var radius = Math.max(0.6, spacing * 0.30 * display.atomSize);
    var clipZ = null;
    if (display.clip < 0.999) {
      clipZ = fitMin[2] + display.clip * (fitMax[2] - fitMin[2]);
    }
    var groups = new Map();
    for (var i = 0; i < atoms.length; i += 1) {
      if (display.hiddenElements[atoms[i].element]) {
        continue;
      }
      var world = atomPoint(i);
      if (clipZ !== null && world[2] > clipZ) {
        continue;
      }
      var p = camera(world);
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
      var fill = base.charAt(0) === "#" ? shadeColor(base, factor) : base;
      ctx.beginPath();
      for (var j = 0; j < group.points.length; j += 1) {
        var q = group.points[j];
        ctx.moveTo(q[0] + radius, q[1]);
        ctx.arc(q[0], q[1], radius, 0, 2 * Math.PI);
      }
      ctx.fillStyle = fill;
      ctx.fill();
      ctx.strokeStyle = "rgba(8, 12, 17, 0.45)";
      ctx.lineWidth = Math.max(0.4, radius * 0.12);
      ctx.stroke();
    });
    ctx.globalAlpha = 1;
  }

  function drawAtomistic(alpha) {
    if (renderer && atomSpacingPx() >= 6) {
      return;
    }
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
      ctx.fillStyle = "rgba(230, 237, 243, 0.95)";
      ctx.font = "11px ui-monospace, Menlo, monospace";
      ctx.fillText(body.name, p[0] + radius + 3, p[1] - radius - 2);
    }
    ctx.globalAlpha = 1;
  }

  /* ---------- Layer weights and frame ---------- */

  function layerWeights() {
    return [
      toggles.atomistic.checked ? 1 : 0,
      toggles.device.checked ? 1 : 0,
      toggles.coarse.checked ? 1 : 0,
    ];
  }

  function draw() {
    ctx.clearRect(0, 0, w, h);
    if (!scene) {
      return;
    }
    var weights = layerWeights();
    var useGl = Boolean(renderer) && weights[0] > 0 && atomSpacingPx() >= 6;
    if (renderer) {
      if (useGl) {
        renderer.render(glState());
      } else {
        renderer.clear();
      }
    }
    if (weights[2] > 0) {
      drawCoarse(weights[2]);
    }
    if (weights[1] > 0) {
      drawDevice(weights[1]);
    }
    if (weights[0] > 0) {
      drawAtomistic(weights[0]);
    }
    drawOverlay();
    drawCaption();
  }

  function drawCaption() {
    ctx.globalAlpha = 1;
    ctx.fillStyle = "rgba(240, 180, 41, 0.9)";
    ctx.font = "bold 12px ui-monospace, Menlo, monospace";
    ctx.fillText("simulated / schematic", 12, 20);
    ctx.fillStyle = "rgba(139, 148, 158, 0.9)";
    ctx.font = "11px ui-monospace, Menlo, monospace";
    ctx.fillText(
      "yaw " + view.yaw.toFixed(2) + " rad  tilt " + view.tilt.toFixed(2) + " rad",
      12,
      36
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
    buildAtomBuffers();
    buildLegend();
    hoverAtom = -1;
    measure = [];
    pickEl.textContent = "";
    measureEl.textContent = "";
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

  var pointerState = { active: false, mode: "rotate", x: 0, y: 0, moved: 0, lastPick: 0 };

  canvas.addEventListener("pointerdown", function (event) {
    pointerState.active = true;
    pointerState.mode = event.shiftKey || event.button === 1 ? "pan" : "rotate";
    pointerState.x = event.clientX;
    pointerState.y = event.clientY;
    pointerState.moved = 0;
    canvas.setPointerCapture(event.pointerId);
  });

  canvas.addEventListener("pointermove", function (event) {
    if (!pointerState.active) {
      var now = performance.now();
      if (pointerState.lastPick && now - pointerState.lastPick < 40) {
        return;
      }
      pointerState.lastPick = now;
      var rect = canvas.getBoundingClientRect();
      var found = pickAtom(event.clientX - rect.left, event.clientY - rect.top);
      if (found !== hoverAtom) {
        hoverAtom = found;
        pickEl.textContent = found >= 0 ? describeAtom(found) : "";
        draw();
      }
      return;
    }
    var dx = event.clientX - pointerState.x;
    var dy = event.clientY - pointerState.y;
    pointerState.moved += Math.abs(dx) + Math.abs(dy);
    if (pointerState.mode === "pan") {
      view.panX += dx;
      view.panY += dy;
    } else {
      view.yaw += dx * 0.01;
      view.tilt += dy * 0.01;
      view.tilt = Math.max(0.05, Math.min(Math.PI - 0.05, view.tilt));
    }
    pointerState.x = event.clientX;
    pointerState.y = event.clientY;
    draw();
  });

  canvas.addEventListener("pointerup", function (event) {
    var wasClick = pointerState.active && pointerState.moved < 4;
    pointerState.active = false;
    if (canvas.hasPointerCapture(event.pointerId)) {
      canvas.releasePointerCapture(event.pointerId);
    }
    if (wasClick) {
      var rect = canvas.getBoundingClientRect();
      var index = pickAtom(event.clientX - rect.left, event.clientY - rect.top);
      if (index >= 0) {
        measure.push(index);
        if (measure.length > 2) {
          measure = [index];
        }
      } else {
        measure = [];
      }
      updateMeasureText();
      draw();
    }
  });

  function updateMeasureText() {
    if (measure.length < 2) {
      measureEl.textContent =
        measure.length === 1 ? "measure: pick a second atom" : "";
      return;
    }
    var a = atomPoint(measure[0]);
    var b = atomPoint(measure[1]);
    var d = Math.sqrt(
      (a[0] - b[0]) * (a[0] - b[0]) +
        (a[1] - b[1]) * (a[1] - b[1]) +
        (a[2] - b[2]) * (a[2] - b[2])
    );
    measureEl.textContent =
      "measure: " + (d * 1e10).toFixed(3) + " A (" + (d * 1e9).toFixed(3) + " nm)";
  }

  canvas.addEventListener(
    "wheel",
    function (event) {
      event.preventDefault();
      if (event.shiftKey) {
        view.panX -= event.deltaY * 0.6;
      } else {
        view.zoom *= event.deltaY < 0 ? 1.1 : 0.9;
        view.zoom = Math.max(0.2, Math.min(8, view.zoom));
      }
      draw();
    },
    { passive: false }
  );

  canvas.addEventListener("dblclick", function () {
    setView("reset");
  });

  function setView(name) {
    var presets = {
      reset: [-0.65, 1.02],
      iso: [-0.65, 1.02],
      top: [0, 0.05],
      front: [0, Math.PI / 2],
    };
    var preset = presets[name] || presets.reset;
    view.yaw = preset[0];
    view.tilt = preset[1];
    view.zoom = 1;
    view.panX = 0;
    view.panY = 0;
    draw();
  }

  buttons.reset.addEventListener("click", function () {
    setView("reset");
  });
  buttons.iso.addEventListener("click", function () {
    setView("iso");
  });
  buttons.top.addEventListener("click", function () {
    setView("top");
  });
  buttons.front.addEventListener("click", function () {
    setView("front");
  });

  buttons.play.addEventListener("click", function () {
    motion.playing = !motion.playing;
    buttons.play.textContent = motion.playing ? "Pause" : "Play";
  });

  buttons.turntable.addEventListener("click", function () {
    motion.turntable = !motion.turntable;
    buttons.turntable.textContent = motion.turntable ? "Stop" : "Turntable";
  });

  buttons.full.addEventListener("click", function () {
    var app = document.getElementById("app");
    if (document.fullscreenElement) {
      document.exitFullscreen();
    } else if (app.requestFullscreen) {
      app.requestFullscreen();
    }
  });

  buttons.save.addEventListener("click", function () {
    var dpr = window.devicePixelRatio || 1;
    var out = document.createElement("canvas");
    out.width = Math.round(w * dpr);
    out.height = Math.round(h * dpr);
    var octx = out.getContext("2d");
    octx.fillStyle = "#0d1117";
    octx.fillRect(0, 0, out.width, out.height);
    if (renderer && glCanvas.width) {
      octx.drawImage(glCanvas, 0, 0, out.width, out.height);
    }
    octx.drawImage(canvas, 0, 0, out.width, out.height);
    out.toBlob(function (blob) {
      if (!blob) {
        return;
      }
      var link = document.createElement("a");
      link.href = URL.createObjectURL(blob);
      link.download = "nano-cad-gearbox.png";
      link.click();
      URL.revokeObjectURL(link.href);
    });
  });

  aoToggle.addEventListener("change", function () {
    display.ao = aoToggle.checked;
    draw();
  });
  atomSizeInput.addEventListener("input", function () {
    display.atomSize = parseFloat(atomSizeInput.value);
    draw();
  });
  clipInput.addEventListener("input", function () {
    display.clip = parseFloat(clipInput.value);
    draw();
  });
  speedInput.addEventListener("input", function () {
    motion.speed = parseFloat(speedInput.value);
  });

  window.addEventListener("keydown", function (event) {
    if (event.target && /input|textarea/i.test(event.target.tagName)) {
      return;
    }
    var panStep = 24;
    if (event.key === "ArrowLeft") {
      view.panX -= panStep;
    } else if (event.key === "ArrowRight") {
      view.panX += panStep;
    } else if (event.key === "ArrowUp") {
      view.panY -= panStep;
    } else if (event.key === "ArrowDown") {
      view.panY += panStep;
    } else if (event.key === "+" || event.key === "=") {
      view.zoom = Math.min(8, view.zoom * 1.1);
    } else if (event.key === "-" || event.key === "_") {
      view.zoom = Math.max(0.2, view.zoom * 0.9);
    } else if (event.key === "r" || event.key === "R") {
      setView("reset");
      return;
    } else if (event.key === "f" || event.key === "F") {
      buttons.full.click();
      return;
    } else if (event.key === "t" || event.key === "T") {
      buttons.turntable.click();
      return;
    } else if (event.key === "p" || event.key === "P") {
      buttons.play.click();
      return;
    } else if (event.key === "Escape") {
      measure = [];
      updateMeasureText();
    } else {
      return;
    }
    draw();
  });

  Object.keys(toggles).forEach(function (name) {
    toggles[name].addEventListener("change", draw);
  });

  window.addEventListener("resize", resize);

  function frame(now) {
    var dt = motion.last ? Math.min(0.05, (now - motion.last) / 1000) : 0;
    motion.last = now;
    if (motion.playing) {
      updateAnimation(dt * motion.speed);
      uploadAtoms();
    }
    if (motion.turntable) {
      view.yaw += dt * 0.5 * motion.speed;
    }
    if (motion.playing || motion.turntable) {
      draw();
    }
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);

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

  window.NanoCadDebug = {
    view: view,
    display: display,
    motion: motion,
    hasRenderer: function () {
      return !!renderer;
    },
    spacingPx: function () {
      return atomSpacingPx();
    },
    zoom: function () {
      return view.zoom;
    },
    fitR: function () {
      return fit.r;
    },
    useGl: function () {
      return !!renderer && layerWeights()[0] > 0 && atomSpacingPx() >= 6;
    },
    positions: function () {
      if (!atomCache) {
        return null;
      }
      return Array.prototype.slice.call(atomCache.position, 0, 6);
    },
    redraw: function () {
      draw();
    },
  };

  load();
  initApp();
})();
