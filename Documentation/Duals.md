# Dual, DualComplex and DualSplitComplex

The `Dual`, `DualComplex` and `DualSplitComplex` classes represent dual numbers with real, complex or [split-complex](SplitComplex+SplitComplexIdempotent.md) coefficients.

A dual number has the form $ a + b \varepsilon $ where $\varepsilon^2 = 0$ and $\varepsilon \ne 0$.

Dual numbers are commutative.

## Custom functions

### `Dual` only

- **E:** returns $\varepsilon$ (`Dual(0, 1)`)
- **expE:** calculates $e^{x\varepsilon}$

## Automatic Differentiation

Dual numbers can compute the first derivative automatically.

$ f(a+b\varepsilon) = f(a) + b f'(a) \varepsilon $

## Examples

$ f(x) = 3x^2+5x-2 $ therefore $ f'(x) = 6x+5 $.

$ f(5) = 3*5^2+5*5-2 = 98  $ and $ f'(5) = 6*5+5 = 35 $.

```python
> f = lambda x : 3 * x**2 + 5 * x - 2
> f(Dual(5, 1))
Dual(98, 35)
```

## Further information

[Wikipedia: Dual number](https://en.wikipedia.org/wiki/Dual_number).