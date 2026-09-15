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

  /* ---------- Layer drawing ---------- */

  function drawAtomistic(alpha) {
    var atoms = scene.atomistic.atoms;
    var groups = new Map();
    for (var i = 0; i < atoms.length; i += 1) {
      var key = atoms[i].element;
      if (!groups.has(key)) {
        groups.set(key, []);
      }
      groups.get(key).push(atoms[i]);
    }
    ctx.globalAlpha = alpha;
    var radius = Math.max(1.0, 0.011 * pixelSize());
    groups.forEach(function (list, element) {
      ctx.fillStyle = elementColor(element);
      ctx.beginPath();
      for (var j = 0; j < list.length; j += 1) {
        var p = camera(list[j].position_m);
        ctx.moveTo(p[0] + radius, p[1]);
        ctx.arc(p[0], p[1], radius, 0, 2 * Math.PI);
      }
      ctx.fill();
    });
    ctx.globalAlpha = 1;
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

  updateSliderLabel();
  load();
})();
