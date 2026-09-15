/*
 * Nano-CAD atomistic renderer (WebGL2).
 *
 * Draws the L1 atoms as instanced sphere impostors with a key/fill/rim lighting
 * model, screen-space ambient occlusion, depth fog, a section clip, and an
 * optional bond layer. It is a still renderer: it does not simulate forces.
 *
 * The projection repeats `camera()` in viewer.js so the atoms land on the same
 * pixels as the 2D schematic. If WebGL2 is unavailable the factory returns
 * null and the viewer keeps its Canvas 2D atom path. No external dependency.
 */
(function (global) {
  "use strict";

  /* ---------- Pure helpers (also used by the tests) ---------- */

  function normalize3(v) {
    var len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) || 1;
    return [v[0] / len, v[1] / len, v[2] / len];
  }

  function cross3(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
  }

  /*
   * A cosine-weighted hemisphere kernel in tangent space (z up). A cosine
   * kernel gives the soft contact darkening of a real ambient term.
   */
  function hemisphereKernel(count) {
    var kernel = [];
    for (var i = 0; i < count; i += 1) {
      var u = (i + 0.5) / count;
      var r = Math.sqrt(u);
      var phi = i * 2.399963229728653;
      var z = Math.sqrt(Math.max(0, 1 - u));
      var scale = 0.35 + 0.65 * u;
      kernel.push([Math.cos(phi) * r * scale, Math.sin(phi) * r * scale, z * scale]);
    }
    return kernel;
  }

  /*
   * Repeat of the vertex projection. `state` carries the camera terms that
   * viewer.js also uses, so the two paths cannot drift.
   */
  function projectAtom(point, state) {
    var cy = Math.cos(state.yaw);
    var sy = Math.sin(state.yaw);
    var ct = Math.cos(state.tilt);
    var st = Math.sin(state.tilt);
    var nx = (point[0] - state.center[0]) * state.invR;
    var ny = (point[1] - state.center[1]) * state.invR;
    var nz = (point[2] - state.center[2]) * state.invR;
    var x1 = nx * cy - ny * sy;
    var y1 = nx * sy + ny * cy;
    var y2 = y1 * ct - nz * st;
    var depth = y1 * st + nz * ct;
    var viewDepth = state.camDistance - depth * state.perspK;
    var persp = state.camDistance / viewDepth;
    return {
      x: state.viewport[0] * 0.5 + x1 * state.scale * persp + state.pan[0],
      y: state.viewport[1] * 0.5 + y2 * state.scale * persp + state.pan[1],
      viewDepth: viewDepth,
      persp: persp,
    };
  }

  /* Reconstruct a view-space point from a normalised depth and a screen uv. */
  function viewPosition(uv, depthNorm, state) {
    var vd = state.near + depthNorm * (state.far - state.near);
    var k = vd / (state.scale * state.camDistance);
    return [
      (uv[0] - 0.5) * state.viewport[0] * k,
      -(uv[1] - 0.5) * state.viewport[1] * k,
      -vd,
    ];
  }

  /* ---------- GLSL ---------- */

  var QUAD_VS = [
    "#version 300 es",
    "out vec2 vUv;",
    "void main() {",
    "  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));",
    "  vUv = p * 0.5;",
    "  gl_Position = vec4(p - 1.0, 0.0, 1.0);",
    "}",
  ].join("\n");

  var ATOM_VS = [
    "#version 300 es",
    "layout(location = 0) in vec3 aPos;",
    "layout(location = 1) in float aRadius;",
    "layout(location = 2) in vec3 aColor;",
    "uniform vec3 uCenter;",
    "uniform float uInvR;",
    "uniform vec2 uYaw;",
    "uniform vec2 uTilt;",
    "uniform float uCamDist;",
    "uniform float uPerspK;",
    "uniform float uScale;",
    "uniform vec2 uViewport;",
    "uniform vec2 uPan;",
    "uniform float uRadiusScale;",
    "uniform float uNear;",
    "uniform float uFar;",
    "uniform float uClipZ;",
    "uniform float uClipOn;",
    "out vec3 vColor;",
    "out vec2 vUV;",
    "out float vViewDepth;",
    "out float vRadiusNorm;",
    "void main() {",
    "  vec3 n = (aPos - uCenter) * uInvR;",
    "  float x1 = n.x * uYaw.x - n.y * uYaw.y;",
    "  float y1 = n.x * uYaw.y + n.y * uYaw.x;",
    "  float y2 = y1 * uTilt.x - n.z * uTilt.y;",
    "  float depth = y1 * uTilt.y + n.z * uTilt.x;",
    "  float viewDepth = uCamDist - depth * uPerspK;",
    "  float persp = uCamDist / viewDepth;",
    "  vec2 screen = vec2(uViewport.x * 0.5, uViewport.y * 0.5)",
    "    + vec2(x1, y2) * uScale * persp + uPan;",
    "  float radiusPx = aRadius * uRadiusScale * uInvR * uScale * persp;",
    "  vec2 corner;",
    "  corner.x = (gl_VertexID == 1 || gl_VertexID == 3) ? 1.0 : -1.0;",
    "  corner.y = (gl_VertexID >= 2) ? 1.0 : -1.0;",
    "  vec2 pos = screen + corner * radiusPx;",
    "  vec2 ndc = vec2(pos.x / uViewport.x * 2.0 - 1.0,",
    "                  1.0 - pos.y / uViewport.y * 2.0);",
    "  gl_Position = vec4(ndc, 0.0, 1.0);",
    "  if (uClipOn > 0.5 && aPos.z > uClipZ) {",
    "    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);",
    "  }",
    "  vUV = corner;",
    "  vColor = aColor;",
    "  vViewDepth = viewDepth;",
    "  vRadiusNorm = aRadius * uRadiusScale * uInvR;",
    "}",
  ].join("\n");

  var ATOM_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec3 vColor;",
    "in vec2 vUV;",
    "in float vViewDepth;",
    "in float vRadiusNorm;",
    "uniform float uNear;",
    "uniform float uFar;",
    "uniform vec3 uKeyDir;",
    "uniform float uBallStick;",
    "out vec4 fragColor;",
    "void main() {",
    "  float r2 = dot(vUV, vUV);",
    "  if (r2 > 1.0) { discard; }",
    "  vec3 normal = vec3(vUV.x, -vUV.y, sqrt(max(0.0, 1.0 - r2)));",
    "  float surfaceDepth = vViewDepth - normal.z * vRadiusNorm;",
    "  gl_FragDepth = (surfaceDepth - uNear) / (uFar - uNear);",
    "  vec3 nrm = normalize(normal);",
    "  vec3 key = normalize(uKeyDir);",
    "  float diff = max(dot(nrm, key), 0.0);",
    "  float ambient = 0.34 + 0.22 * nrm.z;",
    "  vec3 fillDir = normalize(vec3(0.7, -0.4, 0.55));",
    "  float fill = max(dot(nrm, fillDir), 0.0) * 0.25;",
    "  float rim = pow(1.0 - max(nrm.z, 0.0), 3.0) * 0.22;",
    "  vec3 view = vec3(0.0, 0.0, 1.0);",
    "  vec3 halfV = normalize(key + view);",
    "  float spec = pow(max(dot(nrm, halfV), 0.0), 56.0) * (1.0 - uBallStick) * 0.45;",
    "  vec3 col = vColor * (ambient + 0.95 * diff + fill + rim) + vec3(spec);",
    "  float edge = smoothstep(0.0, 0.30, nrm.z);",
    "  col = mix(col * 0.16, col, edge);",
    "  float fogFrac = clamp((surfaceDepth - 2.6) / 1.2, 0.0, 1.0);",
    "  float alpha = mix(0.98, 0.45, fogFrac);",
    "  fragColor = vec4(col, alpha);",
    "}",
  ].join("\n");

  var LINE_VS = [
    "#version 300 es",
    "layout(location = 0) in vec3 aPos;",
    "uniform vec3 uCenter;",
    "uniform float uInvR;",
    "uniform vec2 uYaw;",
    "uniform vec2 uTilt;",
    "uniform float uCamDist;",
    "uniform float uPerspK;",
    "uniform float uScale;",
    "uniform vec2 uViewport;",
    "uniform vec2 uPan;",
    "uniform float uNear;",
    "uniform float uFar;",
    "out float vViewDepth;",
    "void main() {",
    "  vec3 n = (aPos - uCenter) * uInvR;",
    "  float x1 = n.x * uYaw.x - n.y * uYaw.y;",
    "  float y1 = n.x * uYaw.y + n.y * uYaw.x;",
    "  float y2 = y1 * uTilt.x - n.z * uTilt.y;",
    "  float depth = y1 * uTilt.y + n.z * uTilt.x;",
    "  float viewDepth = uCamDist - depth * uPerspK;",
    "  float persp = uCamDist / viewDepth;",
    "  vec2 screen = vec2(uViewport.x * 0.5, uViewport.y * 0.5)",
    "    + vec2(x1, y2) * uScale * persp + uPan;",
    "  vec2 ndc = vec2(screen.x / uViewport.x * 2.0 - 1.0,",
    "                  1.0 - screen.y / uViewport.y * 2.0);",
    "  gl_Position = vec4(ndc, 0.0, 1.0);",
    "  vViewDepth = viewDepth;",
    "}",
  ].join("\n");

  var LINE_FS = [
    "#version 300 es",
    "precision highp float;",
    "in float vViewDepth;",
    "uniform float uNear;",
    "uniform float uFar;",
    "out vec4 fragColor;",
    "void main() {",
    "  gl_FragDepth = (vViewDepth - uNear) / (uFar - uNear);",
    "  fragColor = vec4(0.32, 0.36, 0.42, 0.5);",
    "}",
  ].join("\n");

  var AO_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec2 vUv;",
    "uniform sampler2D uDepth;",
    "uniform vec2 uViewport;",
    "uniform float uScale;",
    "uniform float uCamDist;",
    "uniform float uNear;",
    "uniform float uFar;",
    "uniform float uRadius;",
    "uniform float uStrength;",
    "uniform float uBias;",
    "uniform vec3 uKernel[16];",
    "out vec4 fragColor;",
    "vec3 viewPos(vec2 uv, float dn) {",
    "  float vd = uNear + dn * (uFar - uNear);",
    "  vec2 rel = (uv - 0.5) * uViewport;",
    "  float k = vd / (uScale * uCamDist);",
    "  return vec3(rel.x * k, -rel.y * k, -vd);",
    "}",
    "void main() {",
    "  float dn = texture(uDepth, vUv).r;",
    "  vec3 P = viewPos(vUv, dn);",
    "  vec3 N = normalize(cross(dFdx(P), dFdy(P)));",
    "  if (N.z < 0.0) { N = -N; }",
    "  float angle = fract(sin(dot(vUv * uViewport, vec2(12.9898, 78.233)))",
    "    * 43758.5453) * 6.2831853;",
    "  vec3 rnd = vec3(cos(angle), sin(angle), 0.0);",
    "  vec3 T = normalize(rnd - N * dot(rnd, N));",
    "  vec3 B = cross(N, T);",
    "  mat3 basis = mat3(T, B, N);",
    "  float occ = 0.0;",
    "  for (int i = 0; i < 16; i++) {",
    "    vec3 sp = P + basis * uKernel[i] * uRadius;",
    "    vec4 clip = vec4(",
    "      sp.x / (uScale * uCamDist) / (-sp.z) * uViewport.x + uViewport.x * 0.5,",
    "      uViewport.y * 0.5 - sp.y / (uScale * uCamDist) / (-sp.z) * uViewport.y,",
    "      0.0, 1.0);",
    "    vec2 suv = clip.xy / uViewport;",
    "    if (suv.x < 0.0 || suv.x > 1.0 || suv.y < 0.0 || suv.y > 1.0) { continue; }",
    "    float sd = uNear + texture(uDepth, suv).r * (uFar - uNear);",
    "    vec3 SP = viewPos(suv, texture(uDepth, suv).r);",
    "    float range = smoothstep(0.0, 1.0, uRadius / abs(P.z - SP.z));",
    "    occ += (SP.z >= sp.z + uBias ? 1.0 : 0.0) * range;",
    "  }",
    "  occ = 1.0 - (occ / 16.0) * uStrength;",
    "  fragColor = vec4(occ, occ, occ, 1.0);",
    "}",
  ].join("\n");

  var BLUR_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec2 vUv;",
    "uniform sampler2D uSource;",
    "uniform vec2 uTexel;",
    "out vec4 fragColor;",
    "void main() {",
    "  float sum = 0.0;",
    "  sum += texture(uSource, vUv + uTexel * vec2(-1.0, -1.0)).r * 0.25;",
    "  sum += texture(uSource, vUv + uTexel * vec2(1.0, -1.0)).r * 0.25;",
    "  sum += texture(uSource, vUv + uTexel * vec2(-1.0, 1.0)).r * 0.25;",
    "  sum += texture(uSource, vUv + uTexel * vec2(1.0, 1.0)).r * 0.25;",
    "  fragColor = vec4(sum, sum, sum, 1.0);",
    "}",
  ].join("\n");

  var COMPOSITE_FS = [
    "#version 300 es",
    "precision highp float;",
    "in vec2 vUv;",
    "uniform sampler2D uColor;",
    "uniform sampler2D uAo;",
    "uniform float uAoOn;",
    "uniform float uVignette;",
    "out vec4 fragColor;",
    "void main() {",
    "  vec4 c = texture(uColor, vUv);",
    "  float ao = mix(1.0, texture(uAo, vUv).r, uAoOn);",
    "  vec3 col = c.rgb * ao;",
    "  vec2 q = (vUv - 0.5) * vec2(1.1, 1.02);",
    "  col *= 1.0 - uVignette * dot(q, q);",
    "  fragColor = vec4(col * c.a, c.a);",
    "}",
  ].join("\n");

  /* ---------- GL plumbing ---------- */

  function compile(gl, type, source) {
    var shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      var log = gl.getShaderInfoLog(shader);
      gl.deleteShader(shader);
      throw new Error("shader compile failed: " + log);
    }
    return shader;
  }

  function program(gl, vsSource, fsSource) {
    var vs = compile(gl, gl.VERTEX_SHADER, vsSource);
    var fs = compile(gl, gl.FRAGMENT_SHADER, fsSource);
    var prog = gl.createProgram();
    gl.attachShader(prog, vs);
    gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    gl.deleteShader(vs);
    gl.deleteShader(fs);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
      var log = gl.getProgramInfoLog(prog);
      gl.deleteProgram(prog);
      throw new Error("program link failed: " + log);
    }
    return prog;
  }

  function makeTexture(gl, filter) {
    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return texture;
  }

  /*
   * Build a WebGL2 renderer for `canvas`. Returns null when WebGL2 is missing
   * or a shader fails, so the caller can keep the Canvas 2D path.
   */
  function createAtomRenderer(canvas) {
    var gl = canvas.getContext("webgl2", {
      alpha: true,
      antialias: true,
      premultipliedAlpha: true,
      depth: true,
      preserveDrawingBuffer: true,
    });
    if (!gl) {
      return null;
    }

    var programs;
    try {
      programs = {
        atoms: program(gl, ATOM_VS, ATOM_FS),
        lines: program(gl, LINE_VS, LINE_FS),
        ao: program(gl, AO_FS, QUAD_VS),
        blur: program(gl, BLUR_FS, QUAD_VS),
        composite: program(gl, COMPOSITE_FS, QUAD_VS),
      };
    } catch (error) {
      gl.getExtension("WEBGL_lose_context").loseContext();
      return null;
    }

    var locations = {};
    Object.keys(programs).forEach(function (name) {
      var prog = programs[name];
      var count = gl.getProgramParameter(prog, gl.ACTIVE_UNIFORMS);
      var map = {};
      for (var i = 0; i < count; i += 1) {
        var info = gl.getActiveUniform(prog, i);
        map[info.name.replace(/\[0\]$/, "")] = gl.getUniformLocation(prog, info.name);
      }
      locations[name] = map;
    });

    var atomVao = gl.createVertexArray();
    var atomBuffers = {
      position: gl.createBuffer(),
      radius: gl.createBuffer(),
      color: gl.createBuffer(),
    };
    var lineVao = gl.createVertexArray();
    var lineBuffer = gl.createBuffer();
    var radiiCache = null;
    var colorCache = null;
    var count = 0;
    var lineCount = 0;

    var targets = { color: null, depth: null, ao: null, blur: null };
    var fbo = { scene: gl.createFramebuffer(), ao: gl.createFramebuffer(), blur: gl.createFramebuffer() };
    var size = { w: 0, h: 0 };
    var kernel = hemisphereKernel(16);
    var kernelData = new Float32Array(48);
    kernel.forEach(function (v, i) {
      kernelData[i * 3] = v[0];
      kernelData[i * 3 + 1] = v[1];
      kernelData[i * 3 + 2] = v[2];
    });
    var emptyVao = gl.createVertexArray();

    function attach(target, texture, width, height) {
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texImage2D(
        gl.TEXTURE_2D,
        0,
        target.internal,
        width,
        height,
        0,
        target.format,
        target.type,
        null
      );
    }

    function resize(width, height) {
      var w = Math.max(1, Math.round(width));
      var h = Math.max(1, Math.round(height));
      if (w === size.w && h === size.h) {
        return;
      }
      size.w = w;
      size.h = h;
      ["color", "depth", "ao", "blur"].forEach(function (name) {
        if (targets[name]) {
          gl.deleteTexture(targets[name]);
        }
      });
      var half = Math.max(1, Math.round(w * 0.5));
      var halfH = Math.max(1, Math.round(h * 0.5));
      targets.color = makeTexture(gl, gl.NEAREST);
      targets.depth = makeTexture(gl, gl.NEAREST);
      targets.ao = makeTexture(gl, gl.LINEAR);
      targets.blur = makeTexture(gl, gl.LINEAR);
      attach(
        { internal: gl.RGBA8, format: gl.RGBA, type: gl.UNSIGNED_BYTE },
        targets.color,
        w,
        h
      );
      attach(
        { internal: gl.DEPTH_COMPONENT24, format: gl.DEPTH_COMPONENT, type: gl.UNSIGNED_INT },
        targets.depth,
        w,
        h
      );
      attach({ internal: gl.R8, format: gl.RED, type: gl.UNSIGNED_BYTE }, targets.ao, half, halfH);
      attach({ internal: gl.R8, format: gl.RED, type: gl.UNSIGNED_BYTE }, targets.blur, half, halfH);

      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.scene);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, targets.color, 0);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.TEXTURE_2D, targets.depth, 0);
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.ao);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, targets.ao, 0);
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.blur);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, targets.blur, 0);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }

    function setAtoms(positions, radii, colors) {
      count = positions.length / 3;
      radiiCache = radii;
      colorCache = colors;
      gl.bindVertexArray(atomVao);
      gl.bindBuffer(gl.ARRAY_BUFFER, atomBuffers.position);
      gl.bufferData(gl.ARRAY_BUFFER, positions, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0);
      gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
      gl.vertexAttribDivisor(0, 1);
      gl.bindBuffer(gl.ARRAY_BUFFER, atomBuffers.radius);
      gl.bufferData(gl.ARRAY_BUFFER, radii, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(1);
      gl.vertexAttribPointer(1, 1, gl.FLOAT, false, 0, 0);
      gl.vertexAttribDivisor(1, 1);
      gl.bindBuffer(gl.ARRAY_BUFFER, atomBuffers.color);
      gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(2);
      gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 0, 0);
      gl.vertexAttribDivisor(2, 1);
      gl.bindVertexArray(null);
    }

    function setBonds(segments) {
      lineCount = segments.length / 3;
      gl.bindVertexArray(lineVao);
      gl.bindBuffer(gl.ARRAY_BUFFER, lineBuffer);
      gl.bufferData(gl.ARRAY_BUFFER, segments, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(0);
      gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
      gl.bindVertexArray(null);
    }

    /* Re-upload positions only, for the kinematic animation. */
    function updatePositions(positions) {
      if (positions.length / 3 !== count) {
        setAtoms(positions, radiiCache, colorCache);
        return;
      }
      gl.bindBuffer(gl.ARRAY_BUFFER, atomBuffers.position);
      gl.bufferSubData(gl.ARRAY_BUFFER, 0, positions);
    }

    function setUniforms(map, state) {
      gl.uniform3fv(map.uCenter, state.center);
      gl.uniform1f(map.uInvR, state.invR);
      gl.uniform2f(map.uYaw, Math.cos(state.yaw), Math.sin(state.yaw));
      gl.uniform2f(map.uTilt, Math.cos(state.tilt), Math.sin(state.tilt));
      gl.uniform1f(map.uCamDist, state.camDistance);
      gl.uniform1f(map.uPerspK, state.perspK);
      gl.uniform1f(map.uScale, state.scale);
      gl.uniform2f(map.uViewport, state.viewport[0], state.viewport[1]);
      gl.uniform2f(map.uPan, state.pan[0], state.pan[1]);
      gl.uniform1f(map.uNear, state.near);
      gl.uniform1f(map.uFar, state.far);
    }

    function drawQuad() {
      gl.bindVertexArray(emptyVao);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function render(state) {
      if (!count) {
        return;
      }
      gl.viewport(0, 0, size.w, size.h);

      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.scene);
      gl.clearColor(0, 0, 0, 0);
      gl.clearDepth(1.0);
      gl.enable(gl.DEPTH_TEST);
      gl.depthFunc(gl.LESS);
      gl.disable(gl.BLEND);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

      gl.enable(gl.CULL_FACE);
      gl.cullFace(gl.BACK);

      var atomLoc = locations.atoms;
      gl.useProgram(programs.atoms);
      setUniforms(atomLoc, state);
      gl.uniform1f(atomLoc.uRadiusScale, state.radiusScale);
      gl.uniform3fv(atomLoc.uKeyDir, normalize3(state.keyDir));
      gl.uniform1f(atomLoc.uBallStick, state.ballStick ? 1 : 0);
      gl.uniform1f(atomLoc.uClipZ, state.clipZ);
      gl.uniform1f(atomLoc.uClipOn, state.clipOn ? 1 : 0);
      gl.bindVertexArray(atomVao);
      gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, count);

      if (state.showBonds && lineCount) {
        gl.disable(gl.CULL_FACE);
        gl.useProgram(programs.lines);
        setUniforms(locations.lines, state);
        gl.bindVertexArray(lineVao);
        gl.drawArrays(gl.LINES, 0, lineCount);
      }

      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.ao);
      gl.disable(gl.DEPTH_TEST);
      gl.viewport(0, 0, Math.round(size.w * 0.5), Math.round(size.h * 0.5));
      gl.useProgram(programs.ao);
      var aoLoc = locations.ao;
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, targets.depth);
      gl.uniform1i(aoLoc.uDepth, 0);
      gl.uniform2f(aoLoc.uViewport, state.viewport[0], state.viewport[1]);
      gl.uniform1f(aoLoc.uScale, state.scale);
      gl.uniform1f(aoLoc.uCamDist, state.camDistance);
      gl.uniform1f(aoLoc.uNear, state.near);
      gl.uniform1f(aoLoc.uFar, state.far);
      gl.uniform1f(aoLoc.uRadius, state.aoRadius);
      gl.uniform1f(aoLoc.uStrength, state.aoStrength);
      gl.uniform1f(aoLoc.uBias, 0.0004);
      gl.uniform3fv(aoLoc.uKernel, kernelData);
      drawQuad();

      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.blur);
      gl.useProgram(programs.blur);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, targets.ao);
      gl.uniform1i(locations.blur.uSource, 0);
      gl.uniform2f(locations.blur.uTexel, 1 / (size.w * 0.5), 1 / (size.h * 0.5));
      drawQuad();

      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, size.w, size.h);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
      gl.useProgram(programs.composite);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, targets.color);
      gl.uniform1i(locations.composite.uColor, 0);
      gl.activeTexture(gl.TEXTURE1);
      gl.bindTexture(gl.TEXTURE_2D, targets.blur);
      gl.uniform1i(locations.composite.uAo, 1);
      gl.uniform1f(locations.composite.uAoOn, state.ao ? 1 : 0);
      gl.uniform1f(locations.composite.uVignette, 0.55);
      drawQuad();
      gl.bindVertexArray(null);
    }

    function clear() {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, size.w, size.h);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }

    return {
      resize: resize,
      setAtoms: setAtoms,
      setBonds: setBonds,
      updatePositions: updatePositions,
      render: render,
      clear: clear,
    };
  }

  var api = {
    createAtomRenderer: createAtomRenderer,
    hemisphereKernel: hemisphereKernel,
    projectAtom: projectAtom,
    viewPosition: viewPosition,
    normalize3: normalize3,
    cross3: cross3,
  };

  global.NanoCadAtoms = api;
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  }
})(typeof window !== "undefined" ? window : globalThis);
