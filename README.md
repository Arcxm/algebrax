# algebrax

`algebrax` is a collection of *hypercomplex algebras* implemented in python.

## Why?

I've been fascinated by number systems like these for a long time and wanted to implement them to experiment and explore. Python was chosen because it allows quick prototyping and is used for numerical computing.

Over time, this collection grew to include more systems and functionality for different algebras. Now i want to make it available to the public.

I had a lot of fun learning, researching and experimenting in the REPL. What about you? :)

## What's inside?

For further documentation on the algebras see the [docs](Documentation).
Also see the Docstrings in the code itself as they contain the formulas used and notes on the domain of a function.

| Dimension | Algebras |
| --- | --- |
| 2D | [Dual numbers](Documentation/Duals.md), [Split-complex (including idempotent basis)](Documentation/SplitComplex+SplitComplexIdempotent.md) |
| 3D | [Vec3 (utility for quaternions)](Documentation/Vec3.md) |
| 4D | [Dual complex](Documentation/Duals.md), [Dual split-complex](Documentation/Duals.md), [HyperDual (2nd order duals)](Documentation/HyperDuals.md), Quaternions, Split-quaternions |
| 8D | [Complex hyperduals](Documentation/HyperDuals.md), Dual quaternions, Dual split-quaternions |

Each class implements ([Common Interface Doc](Documentation/Common.md)):

- **Readability:** `__str__`, `__repr__`
- **Comparison operations:** `==`, `!=`
- **Arithmetic operations:** negation, `+`, `-`, `*`, `/`, `(r)pow`
- **Other operations:** conjugate, norm (`__abs__`), inverse
- **Functions:** `exp`, `log`, `sqrt`, `root`

Some classes also provide additional custom functions (e.g. `from_vec3` for quaternions), these are listed in the docs.

Please note that non-commutative algebras do not implement `(r)truediv` and require multiplication with the inverse for clarity (might change in the future).

## Packages

This collection is split into 4 packages.

### dual

The `dual` package contains algebras for automatic differentiation.

Included are the Dual numbers and HyperDuals (2nd order dual numbers).

Dual numbers are available with real, complex and split-complex coefficients.

HyperDuals are available with real and complex coefficients.

### quat

The `quat` package contains algebras for rotation and translation in 3D.

Included is a helper class for 3D Vectors, the Quaternions and the Dual quaternions.

### split

The `split` package contains split-algebras (also called hyperbolic numbers).

Included are the Split-complex numbers (canonical and idempotent basis).

Both classes support mixed arithmetic. The result will be of the type of the left-hand operand. Conversion functions `to_splitcomplex` and `to_idempotent` are available as well.

### split_quat

The `split_quat` package also contains split-algebras.

Included are the Split-quaternions and the Dual split-quaternions.

Split-quaternions extend quaternions to Minkowski space, allowing representation of rotations and boosts in a relativistic context.

## TODO / Future

- Documentation for `Quat`, `SplitQuat`, `DualQuat` and `DualSplitQuat`
- More algebras (e.g. Grassmann)