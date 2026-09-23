/*
 * Unit tests for site/solid_mesh.js.
 *
 * Run with: node site/solid_mesh.test.js
 */
"use strict";

var assert = require("assert");
var solid = require("./solid_mesh.js");

var passed = 0;
var failed = 0;

function test(name, fn) {
  try {
    fn();
    passed += 1;
    console.log("ok   " + name);
  } catch (error) {
    failed += 1;
    console.log("FAIL " + name + "\n     " + error.message);
  }
}

var A = 1.0e-10;

// Every edge of a closed, consistently wound mesh appears once in each direction.
function edgeReport(mesh) {
  var directed = new Map();
  var idx = mesh.indices;
  for (var t = 0; t < idx.length; t += 3) {
    for (var e = 0; e < 3; e += 1) {
      var a = idx[t + e];
      var b = idx[t + ((e + 1) % 3)];
      var key = a + ":" + b;
      directed.set(key, (directed.get(key) || 0) + 1);
    }
  }
  var open = 0;
  var doubled = 0;
  directed.forEach(function (n, key) {
    var parts = key.split(":");
    if (n > 1) doubled += 1;
    if (!directed.has(parts[1] + ":" + parts[0])) open += 1;
  });
  return { open: open, doubled: doubled };
}

function centroid(points) {
  var c = [0, 0, 0];
  var n = points.length / 3;
  for (var i = 0; i < n; i += 1) {
    c[0] += points[i * 3] / n;
    c[1] += points[i * 3 + 1] / n;
    c[2] += points[i * 3 + 2] / n;
  }
  return c;
}

// Signed volume by the divergence theorem. Positive when faces wind outward.
function signedVolume(mesh) {
  var p = mesh.positions;
  var idx = mesh.indices;
  var v = 0;
  for (var t = 0; t < idx.length; t += 3) {
    var a = idx[t] * 3;
    var b = idx[t + 1] * 3;
    var c = idx[t + 2] * 3;
    v +=
      (p[a] * (p[b + 1] * p[c + 2] - p[b + 2] * p[c + 1]) -
        p[a + 1] * (p[b] * p[c + 2] - p[b + 2] * p[c]) +
        p[a + 2] * (p[b] * p[c + 1] - p[b + 1] * p[c])) /
      6;
  }
  return v;
}

test("a lone atom gives a closed sphere of the predicted radius", function () {
  var mesh = solid.buildSurface([1.0e-9, -2.0e-9, 0.5e-9], null, { spacingA: 0.25 });
  assert.ok(mesh.triangleCount > 100, "triangles " + mesh.triangleCount);
  var expected = solid.isolatedRadiusA(1) * A;
  assert.ok(Math.abs(expected - 1.73e-10) < 0.01e-10, "radius " + expected);
  for (var i = 0; i < mesh.vertexCount; i += 1) {
    var dx = mesh.positions[i * 3] - 1.0e-9;
    var dy = mesh.positions[i * 3 + 1] + 2.0e-9;
    var dz = mesh.positions[i * 3 + 2] - 0.5e-9;
    var r = Math.hypot(dx, dy, dz);
    assert.ok(Math.abs(r - expected) < 0.08 * A, "vertex radius " + r);
  }
  var edges = edgeReport(mesh);
  assert.strictEqual(edges.open, 0, "open edges");
  assert.strictEqual(edges.doubled, 0, "doubled edges");
});

test("faces and normals point outward", function () {
  var mesh = solid.buildSurface([0, 0, 0], null, { spacingA: 0.3 });
  var r = solid.isolatedRadiusA(1) * A;
  var volume = signedVolume(mesh);
  var sphere = (4 / 3) * Math.PI * r * r * r;
  assert.ok(volume > 0.9 * sphere && volume < 1.1 * sphere, "volume " + volume / sphere);
  for (var i = 0; i < mesh.vertexCount; i += 1) {
    var dot =
      mesh.positions[i * 3] * mesh.normals[i * 3] +
      mesh.positions[i * 3 + 1] * mesh.normals[i * 3 + 1] +
      mesh.positions[i * 3 + 2] * mesh.normals[i * 3 + 2];
    assert.ok(dot > 0, "normal " + i + " points inward");
  }
});

test("a bonded chain merges into one envelope", function () {
  var points = [];
  for (var i = 0; i < 6; i += 1) points.push(i * 1.54e-10, 0, 0);
  var mesh = solid.buildSurface(points, null, { spacingA: 0.3 });
  var edges = edgeReport(mesh);
  assert.strictEqual(edges.open, 0);
  var c = centroid(mesh.positions);
  assert.ok(Math.abs(c[0] - 2.5 * 1.54e-10) < 0.05e-10, "centroid x " + c[0]);
  assert.ok(Math.abs(c[1]) < 0.05e-10 && Math.abs(c[2]) < 0.05e-10);
});

test("a weak hydrogen gives a smaller bump than a carbon", function () {
  var h = solid.isolatedRadiusA(0.6);
  var c = solid.isolatedRadiusA(1);
  assert.ok(h > 0 && h < c, "h " + h + " c " + c);
  assert.strictEqual(solid.isolatedRadiusA(0.4), 0, "below iso gives nothing");
});

test("meshes group by body", function () {
  var atoms = [
    { element: "C", position_m: [0, 0, 0], body: 0 },
    { element: "C", position_m: [1.0e-9, 0, 0], body: 3 },
    { element: "H", position_m: [1.1e-9, 0, 0], body: 3 },
  ];
  var meshes = solid.buildBodyMeshes(atoms, { spacingA: 0.4 });
  assert.deepStrictEqual(Object.keys(meshes).sort(), ["0", "3"]);
  for (var i = 0; i < meshes["0"].vertexCount; i += 1) {
    assert.ok(meshes["0"].positions[i * 3] < 0.5e-9, "body 0 stays near its atom");
  }
  assert.ok(meshes["3"].triangleCount > 0);
});

test("an empty body gives an empty mesh", function () {
  var mesh = solid.buildSurface([], null);
  assert.strictEqual(mesh.vertexCount, 0);
  assert.strictEqual(mesh.triangleCount, 0);
});

test("the grid spacing stays in its clamp", function () {
  var points = [];
  for (var i = 0; i < 50; i += 1) points.push(i * 5.0e-10, 0, 0);
  var mesh = solid.buildSurface(points, null);
  assert.ok(mesh.spacingM >= 0.6e-10 - 1e-20 && mesh.spacingM <= 1.6e-10 + 1e-20);
});

console.log("\n" + passed + " passed, " + failed + " failed");
if (failed > 0) {
  process.exit(1);
}
