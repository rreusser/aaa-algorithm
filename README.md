# [WIP] aaa-algorithm

> The Adaptive Antoulas-Anderson (AAA) Algorithm for rational approximation

## Introduction

Given function `f(z)`, this module takes as input values `f[i]` tabulated at points `z[i]` and computes a rational approximation using the Adaptive Antoulas-Anderson (AAA) algorithm [\[1\]](#ref1). Along with this, the algorithm returns the poles, zeros, and residues of the approximation.

It's a remarkably effective algorithm [\[2\]](#ref2) which both does its job very well and has also found application in a number of different problems.

## Installation

Not (yet?) published to npm

## TO DO

- [x] compute rational approximation
- [x] evaluate approximation
- [ ] compute poles, zeros, residues
- [ ] clean up Froissart doublets

## Example

```javascript
import AAA from "./aaa.js";
await AAA.ready();

const { VectorComplex, aaa } = AAA;
const { PI, sin, cos, sinh, cosh, exp } = Math;

// Z = exp(linspace(-.5, .5+15i*pi, 1000));
const N = 1000;
const z = [...Array(N).keys()].map((i) => {
  const [tr, ti] = [-0.5 + i / (N - 1), 15.0 * PI * (i / (N - 1))];
  return [exp(tr) * cos(ti), exp(tr) * sin(ti)];
});

// F = tan(pi*z/2);
const f = z.map(([zr, zi]) => [
  sin(PI * zr) / (cosh(PI * zi) + cos(PI * zr)),
  sinh(PI * zi) / (cosh(PI * zi) + cos(PI * zr)),
]);

// Compute the AAA approximation
const approx = aaa(z, f, 1e-13, 100);

// Evaluate the result
const fEval = approx.eval([0.5, 0]);
console.log(`tan(π/4) ~ ${fEval[0]} + ${fEval[1]}i`);

// Clean up
approx.delete();
```

Upon running, this prints the approximated values of `tan(π/2 z)`:

```
f(-0.5 + 0i) = -0.999999999999995 + 1.888073305146161e-15i
f(0 + 0i) = 8.654269559363806e-16 + 2.4482177643400466e-15i
f(0.5 + 0i) = 1.0000000000000009 + 2.4281095898961374e-15i
```

## Development

To build the WASM module via Emscripten, all you (hopefully) need is a working Docker installation.

```
cd aaa/src
make build-docker
make build
```

## References

<a name="ref1"></a>1. Nakatsukasa, Y., Sete, O., & Trefethen, L. N. (2018). The AAA algorithm for rational approximation. SIAM Journal on Scientific Computing, 40(3), A1494–A1522. https://doi.org/10.1137/16M1106122

<a name="ref2"></a>2. Nakatsukasa, Y., Sète, O., & Trefethen, L. N. (2023). The first five years of the AAA algorithm. arXiv. https://arxiv.org/abs/2312.03565

## License

&copy; 2025 Ricky Reusser. MPL-2.0 or MIT License.

This project uses the Eigen library, which is licensed under the Mozilla Public License 2.0 (MPL-2.0). While this project is licensed under the MIT License, the parts that rely on Eigen are governed by the MPL-2.0. No modifications to the Eigen library are made by this project.

For more information about the Eigen library and its licensing terms, visit: https://eigen.tuxfamily.org
