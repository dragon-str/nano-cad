/*
 * Solid (CAD-style) surface meshes built from atom positions.
 *
 * Far from the model, single atoms are smaller than a pixel. The viewer then
 * shows each component as a closed surface instead. This file turns one body's
 * atoms into that surface. It runs in the browser and in node, so the tests in
 * site/solid_mesh.test.js exercise the same code.
 *
 * Method: every atom adds a smooth bump k(d) = w * (1 - d^2 / R^2)^2 to a
 * scalar field on a regular grid. The surface is the level set field = ISO.
 * Surface nets place one vertex in each grid cell that the level set crosses,
 * at the mean of the crossings on the cell edges, and join the vertices of the
 * four cells around every crossed grid edge into a quad. The result is a closed,
 * consistently wound mesh that follows the atoms' envelope, not their spheres.
 *
 * Units: positions in and out are metres. The field is built in angstrom
 * because the numbers are easier to reason about.
 */
(function (global) {
  "use strict";

  var ANGSTROM_PER_M = 1.0e10;

  // Bump radius and iso level. A lone carbon gives a sphere of radius
  // R * sqrt(1 - sqrt(ISO)) = 1.73 A; a flat diamond face sits about 1.9 A out.
  var DEFAULTS = {
    radiusA: 3.2,
    iso: 0.5,
    hydrogenWeight: 0.6,
    maxCells: 160,
    minSpacingA: 0.6,
    maxSpacingA: 1.6,
  };

  function option(options, name) {
    return options && options[name] !== undefined ? options[name] : DEFAULTS[name];
  }

  // Radius of the iso surface around one isolated atom of weight w.
  function isolatedRadiusA(weight, options) {
    var radius = option(options, "radiusA");
    var ratio = option(options, "iso") / weight;
    if (ratio >= 1) {
      return 0;
    }
    return radius * Math.sqrt(1 - Math.sqrt(ratio));
  }

  /*
   * Build the mesh of one body.
   *
   * positionsM: flat [x, y, z, ...] in metres.
   * weights: optional per-atom weights (default 1). Use hydrogenWeight for H.
   * Returns { positions, normals, indices, vertexCount, triangleCount,
   *           spacingM, dims } with Float32Array positions and normals in metres.
   */
  function buildSurface(positionsM, weights, options) {
    var count = Math.floor(positionsM.length / 3);
    var radius = option(options, "radiusA");
    var iso = option(options, "iso");
    if (count === 0) {
      return emptyMesh();
    }

    var lo = [Infinity, Infinity, Infinity];
    var hi = [-Infinity, -Infinity, -Infinity];
    for (var a = 0; a < count; a += 1) {
      for (var axis = 0; axis < 3; axis += 1) {
        var value = positionsM[a * 3 + axis] * ANGSTROM_PER_M;
        if (value < lo[axis]) lo[axis] = value;
        if (value > hi[axis]) hi[axis] = value;
      }
    }
    var extent = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) + 2 * radius;
    var h = extent / option(options, "maxCells");
    h = Math.min(Math.max(h, option(options, "minSpacingA")), option(options, "maxSpacingA"));
    if (options && options.spacingA) {
      h = options.spacingA;
    }

    var pad = radius + 2 * h;
    var origin = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
    var nx = Math.ceil((hi[0] - lo[0] + 2 * pad) / h) + 1;
    var ny = Math.ceil((hi[1] - lo[1] + 2 * pad) / h) + 1;
    var nz = Math.ceil((hi[2] - lo[2] + 2 * pad) / h) + 1;
    var field = new Float32Array(nx * ny * nz);
    var strideY = nx;
    var strideZ = nx * ny;
    var radius2 = radius * radius;
    var reach = Math.ceil(radius / h);

    for (var p = 0; p < count; p += 1) {
      var weight = weights ? weights[p] : 1;
      if (!(weight > 0)) {
        continue;
      }
      var px = positionsM[p * 3] * ANGSTROM_PER_M - origin[0];
      var py = positionsM[p * 3 + 1] * ANGSTROM_PER_M - origin[1];
      var pz = positionsM[p * 3 + 2] * ANGSTROM_PER_M - origin[2];
      var ci = Math.round(px / h);
      var cj = Math.round(py / h);
      var ck = Math.round(pz / h);
      var i0 = Math.max(0, ci - reach);
      var i1 = Math.min(nx - 1, ci + reach);
      var j0 = Math.max(0, cj - reach);
      var j1 = Math.min(ny - 1, cj + reach);
      var k0 = Math.max(0, ck - reach);
      var k1 = Math.min(nz - 1, ck + reach);
      for (var k = k0; k <= k1; k += 1) {
        var dz = k * h - pz;
        var dz2 = dz * dz;
        if (dz2 >= radius2) continue;
        for (var j = j0; j <= j1; j += 1) {
          var dy = j * h - py;
          var dyz2 = dy * dy + dz2;
          if (dyz2 >= radius2) continue;
          var row = k * strideZ + j * strideY;
          for (var i = i0; i <= i1; i += 1) {
            var dx = i * h - px;
            var d2 = dx * dx + dyz2;
            if (d2 >= radius2) continue;
            var q = 1 - d2 / radius2;
            field[row + i] += weight * q * q;
          }
        }
      }
    }

    return surfaceNets(field, nx, ny, nz, h, origin, iso);
  }

  function emptyMesh() {
    return {
      positions: new Float32Array(0),
      normals: new Float32Array(0),
      indices: new Uint32Array(0),
      vertexCount: 0,
      triangleCount: 0,
      spacingM: 0,
      dims: [0, 0, 0],
    };
  }

  // Corner offsets of a cell, and the 12 edges as corner pairs.
  var CORNERS = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
    [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
  ];
  var EDGES = [
    [0, 1], [2, 3], [4, 5], [6, 7],
    [0, 2], [1, 3], [4, 6], [5, 7],
    [0, 4], [1, 5], [2, 6], [3, 7],
  ];

  function surfaceNets(field, nx, ny, nz, h, origin, iso) {
    var strideY = nx;
    var strideZ = nx * ny;
    var cx = nx - 1;
    var cy = ny - 1;
    var cz = nz - 1;
    var cellVertex = new Int32Array(cx * cy * cz).fill(-1);
    var verts = [];
    var norms = [];
    var corner = new Float64Array(8);

    function at(i, j, k) {
      return field[k * strideZ + j * strideY + i];
    }

    // Central-difference gradient at a grid node, clamped at the border.
    function gradient(i, j, k, out) {
      var ia = Math.max(i - 1, 0);
      var ib = Math.min(i + 1, nx - 1);
      var ja = Math.max(j - 1, 0);
      var jb = Math.min(j + 1, ny - 1);
      var ka = Math.max(k - 1, 0);
      var kb = Math.min(k + 1, nz - 1);
      out[0] = (at(ib, j, k) - at(ia, j, k)) / ((ib - ia) * h);
      out[1] = (at(i, jb, k) - at(i, ja, k)) / ((jb - ja) * h);
      out[2] = (at(i, j, kb) - at(i, j, ka)) / ((kb - ka) * h);
    }

    var g = [0, 0, 0];
    for (var k = 0; k < cz; k += 1) {
      for (var j = 0; j < cy; j += 1) {
        for (var i = 0; i < cx; i += 1) {
          var mask = 0;
          for (var c = 0; c < 8; c += 1) {
            var o = CORNERS[c];
            corner[c] = at(i + o[0], j + o[1], k + o[2]);
            if (corner[c] > iso) mask |= 1 << c;
          }
          if (mask === 0 || mask === 255) continue;

          var sx = 0;
          var sy = 0;
          var sz = 0;
          var crossings = 0;
          for (var e = 0; e < 12; e += 1) {
            var a = EDGES[e][0];
            var b = EDGES[e][1];
            var insideA = (mask >> a) & 1;
            var insideB = (mask >> b) & 1;
            if (insideA === insideB) continue;
            var t = (iso - corner[a]) / (corner[b] - corner[a]);
            var oa = CORNERS[a];
            var ob = CORNERS[b];
            sx += oa[0] + t * (ob[0] - oa[0]);
            sy += oa[1] + t * (ob[1] - oa[1]);
            sz += oa[2] + t * (ob[2] - oa[2]);
            crossings += 1;
          }
          var fx = sx / crossings;
          var fy = sy / crossings;
          var fz = sz / crossings;

          // Normal: minus the trilinear-interpolated field gradient.
          var nxv = 0;
          var nyv = 0;
          var nzv = 0;
          for (var cc = 0; cc < 8; cc += 1) {
            var oc = CORNERS[cc];
            var wgt =
              (oc[0] ? fx : 1 - fx) * (oc[1] ? fy : 1 - fy) * (oc[2] ? fz : 1 - fz);
            gradient(i + oc[0], j + oc[1], k + oc[2], g);
            nxv -= wgt * g[0];
            nyv -= wgt * g[1];
            nzv -= wgt * g[2];
          }
          var length = Math.hypot(nxv, nyv, nzv) || 1;

          cellVertex[k * cx * cy + j * cx + i] = verts.length / 3;
          verts.push(
            (origin[0] + (i + fx) * h) / ANGSTROM_PER_M,
            (origin[1] + (j + fy) * h) / ANGSTROM_PER_M,
            (origin[2] + (k + fz) * h) / ANGSTROM_PER_M
          );
          norms.push(nxv / length, nyv / length, nzv / length);
        }
      }
    }

    function cell(i, j, k) {
      return cellVertex[k * cx * cy + j * cx + i];
    }

    var indices = [];
    function quad(v0, v1, v2, v3, flip) {
      if (v0 < 0 || v1 < 0 || v2 < 0 || v3 < 0) return;
      if (flip) {
        indices.push(v0, v3, v2, v0, v2, v1);
      } else {
        indices.push(v0, v1, v2, v0, v2, v3);
      }
    }

    // One quad per grid edge that the surface crosses. The field is larger
    // inside, so an edge whose far end is inside needs the reversed winding.
    for (var kk = 0; kk < nz; kk += 1) {
      for (var jj = 0; jj < ny; jj += 1) {
        for (var ii = 0; ii < nx; ii += 1) {
          var inside0 = at(ii, jj, kk) > iso;
          if (ii < cx && jj > 0 && kk > 0 && jj < cy && kk < cz) {
            var inside1 = at(ii + 1, jj, kk) > iso;
            if (inside0 !== inside1) {
              quad(
                cell(ii, jj - 1, kk - 1),
                cell(ii, jj, kk - 1),
                cell(ii, jj, kk),
                cell(ii, jj - 1, kk),
                inside1
              );
            }
          }
          if (jj < cy && ii > 0 && kk > 0 && ii < cx && kk < cz) {
            var insideY = at(ii, jj + 1, kk) > iso;
            if (inside0 !== insideY) {
              quad(
                cell(ii - 1, jj, kk - 1),
                cell(ii - 1, jj, kk),
                cell(ii, jj, kk),
                cell(ii, jj, kk - 1),
                insideY
              );
            }
          }
          if (kk < cz && ii > 0 && jj > 0 && ii < cx && jj < cy) {
            var insideZ = at(ii, jj, kk + 1) > iso;
            if (inside0 !== insideZ) {
              quad(
                cell(ii - 1, jj - 1, kk),
                cell(ii, jj - 1, kk),
                cell(ii, jj, kk),
                cell(ii - 1, jj, kk),
                insideZ
              );
            }
          }
        }
      }
    }

    return {
      positions: new Float32Array(verts),
      normals: new Float32Array(norms),
      indices: new Uint32Array(indices),
      vertexCount: verts.length / 3,
      triangleCount: indices.length / 3,
      spacingM: h / ANGSTROM_PER_M,
      dims: [nx, ny, nz],
    };
  }

  /*
   * Build one mesh per body from a scene's atom list.
   * atoms: [{ element, position_m, body }]. Returns a map body -> mesh.
   */
  function buildBodyMeshes(atoms, options) {
    var groups = {};
    for (var i = 0; i < atoms.length; i += 1) {
      var atom = atoms[i];
      var key = atom.body === undefined || atom.body === null ? 0 : atom.body;
      if (!groups[key]) groups[key] = { points: [], weights: [] };
      var group = groups[key];
      group.points.push(atom.position_m[0], atom.position_m[1], atom.position_m[2]);
      group.weights.push(atom.element === "H" ? option(options, "hydrogenWeight") : 1);
    }
    var meshes = {};
    Object.keys(groups).forEach(function (key) {
      meshes[key] = buildSurface(groups[key].points, groups[key].weights, options);
    });
    return meshes;
  }

  var api = {
    buildSurface: buildSurface,
    buildBodyMeshes: buildBodyMeshes,
    isolatedRadiusA: isolatedRadiusA,
    DEFAULTS: DEFAULTS,
  };

  global.NanoCadSolid = api;
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  }
})(typeof window !== "undefined" ? window : globalThis);
