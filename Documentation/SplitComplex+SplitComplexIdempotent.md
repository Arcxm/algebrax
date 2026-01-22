# SplitComplex and SplitComplexIdempotent

The `SplitComplex` and `SplitComplexIdempotent` classes represent split-complex numbers (hyperbolic numbers) in their canonical and idempotent basis.

Split-complex numbers are commutative.

Both classes support mixed arithmetic. The result will be of the type of the left-hand operand.

## Bases

### Canonical

In the canonical base, a split-complex number has the form $a+bj$ with $a,b\in\mathbb{R}$ and $j^2=1$ but $j\ne\pm1$.

### Idempotent

In the idempotent base, a split-complex number has the form $ae+be^*$ where $e=\frac{1}{2}(1-j)$ and $e^\*=\frac{1}{2}(1+j)$.

Note that $ee=e$, $e^\*e^\*=e^\*$ and $ee^\*=0$.

## Custom functions

- **sinh:** calculates $sinh(z)$
- **cosh:** calculates $cosh(z)$

### `SplitComplex` only

- **to_idempotent:** converts from the canonical to the idempotent base
- **J:** returns $j$ (`SplitComplex(0, 1)`)
- **expJ:** calculates $e^{xj}$

### `SplitComplexIdempotent` only

- **to_splitcomplex:** converts from the idempotent to the canonical base
- **E:** returns $e$ (`SplitComplexIdempotent(1, 0)`)
- **E2:** returns $e^*$ (`SplitComplexIdempotent(0, 1)`)

## Examples

```python
> z1 = SplitComplex(3, 1)
> z2 = SplitComplexIdempotent(2, 4)
>
> z3 = z1 + z2 # converted automatically
> z3
SplitComplex(6.0, 2.0)
>
> z1i = z1.to_idempotent() # conversion
> z1i
SplitComplexIdempotent(2, 4)
```

## Further information

[Wikipedia: Split-complex number](https://en.wikipedia.org/wiki/Split-complex_number).