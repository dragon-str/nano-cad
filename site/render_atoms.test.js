/*
 * Unit tests for the pure math in site/render_atoms.js.
 *
 * Run with: node site/render_atoms.test.js
 *
 * WebGL cannot run here, so this covers the projection and the AO kernel only.
 * It shares the exact functions the shaders use.
 */
"use strict";

var assert = require("assert");
var atoms = require("./render_atoms.js");

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

function state(overrides) {
  var base = {
    center: [0, 0, 0],
    invR: 1,
    yaw: -0.65,
    tilt: 1.02,
    camDistance: 3.2,
    perspK: 0.6,
    scale: 240,
    viewport: [800, 600],
    pan: [0, 0],
    near: 2.0,
    far: 4.5,
  };
  Object.keys(overrides || {}).forEach(function (key) {
    base[key] = overrides[key];
  });
  return base;
}

test("the hemisphere kernel has the requested size", function () {
  assert.strictEqual(atoms.hemisphereKernel(16).length, 16);
  assert.strictEqual(atoms.hemisphereKernel(1).length, 1);
});

test("every kernel sample is in the upper hemisphere", function () {
  atoms.hemisphereKernel(16).forEach(function (v) {
    assert.ok(v[2] >= 0, "z is negative: " + v[2]);
    var len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    assert.ok(len > 0 && len <= 1.0001, "length out of range: " + len);
  });
});

test("the projection puts the fit centre at the viewport centre", function () {
  var p = atoms.projectAtom([0, 0, 0], state());
  assert.ok(Math.abs(p.x - 400) < 1e-9, "x is " + p.x);
  assert.ok(Math.abs(p.y - 300) < 1e-9, "y is " + p.y);
});

test("the pan offset shifts the projection by the same pixels", function () {
  var p = atoms.projectAtom([0, 0, 0], state({ pan: [30, -20] }));
  assert.ok(Math.abs(p.x - 430) < 1e-9, "x is " + p.x);
  assert.ok(Math.abs(p.y - 280) < 1e-9, "y is " + p.y);
});

test("the perspective factor matches the 2D camera", function () {
  var near = atoms.projectAtom([0, 0, 0.5], state());
  var far = atoms.projectAtom([0, 0, -0.5], state());
  var depth = 0.5 * Math.cos(1.02);
  var expectedNear = 3.2 / (3.2 - depth * 0.6);
  assert.ok(Math.abs(near.persp - expectedNear) < 1e-12, "persp is " + near.persp);
  assert.ok(far.persp < near.persp, "the far point must be smaller");
});

test("the radius grows with the perspective factor", function () {
  var p = atoms.projectAtom([0, 0, 0.4], state());
  assert.ok(p.persp > 1, "the front point must be larger than one");
});

test("the yaw is a rotation and keeps the radius", function () {
  var a = atoms.projectAtom([0.3, -0.2, 0.1], state({ tilt: 0, yaw: 0 }));
  var b = atoms.projectAtom([0.3, -0.2, 0.1], state({ tilt: 0, yaw: 0.7 }));
  var ra = Math.sqrt((a.x - 400) * (a.x - 400) + (a.y - 300) * (a.y - 300));
  var rb = Math.sqrt((b.x - 400) * (b.x - 400) + (b.y - 300) * (b.y - 300));
  assert.ok(Math.abs(ra - rb) < 1e-9, "radius changed: " + ra + " vs " + rb);
});

test("the depth reconstruction is centred at the screen centre", function () {
  var v = atoms.viewPosition([0.5, 0.5], 0.5, state());
  assert.ok(Math.abs(v[0]) < 1e-9, "x is " + v[0]);
  assert.ok(Math.abs(v[1]) < 1e-9, "y is " + v[1]);
  assert.ok(v[2] < 0, "the view depth must be behind the camera");
});

test("a nearer surface has a smaller normalised depth", function () {
  var near = atoms.viewPosition([0.4, 0.5], 0.2, state());
  var far = atoms.viewPosition([0.4, 0.5], 0.8, state());
  assert.ok(near[2] > far[2], "the nearer point has a larger z (less negative)");
});

test("normalize3 returns a unit vector", function () {
  var v = atoms.normalize3([3, -4, 12]);
  var len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  assert.ok(Math.abs(len - 1) < 1e-12, "length is " + len);
});

test("cross3 is right-handed", function () {
  var v = atoms.cross3([1, 0, 0], [0, 1, 0]);
  assert.deepStrictEqual(v, [0, 0, 1]);
});

console.log("\n" + passed + " passed, " + failed + " failed");
if (failed > 0) {
  process.exit(1);
}
