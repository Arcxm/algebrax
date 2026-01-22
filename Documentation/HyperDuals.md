# HyperDual and HyperDualComplex

The `HyperDual` and `HyperDualComplex` classes represent hyperdual numbers with real or complex coefficients.

A hyperdual number has the form $ a + b \varepsilon_1 + c \varepsilon_2 + d \varepsilon_1 \varepsilon_2 $ where $\varepsilon_1^2 = \varepsilon_2^2 = (\varepsilon_1 \varepsilon_2)^2 = 0$ and $\varepsilon_1 \varepsilon_2 \ne 0$.

Hyperdual numbers are commutative.

## Custom functions

- **E1:** returns $\varepsilon_1$ (`HyperDual(0, 1, 0, 0)`)
- **E2:** returns $\varepsilon_2$ (`HyperDual(0, 0, 1, 0)`)
- **E1E2:** returns $\varepsilon_1 \varepsilon_2$ (`HyperDual(0, 0, 0, 1)`)

## Automatic Differentiation

Hyperdual numbers can compute the first and second derivative automatically.

$ f(a+b\varepsilon_1+c\varepsilon_2+d\varepsilon_1\varepsilon_2) = f(a) + b f'(a) \varepsilon_1 + c f'(a) \varepsilon_2 + (bc f''(a) + d f'(a)) \varepsilon_1 \varepsilon_2 $

Special case with $ d = 0 $:

$ f(a+b\varepsilon_1+c\varepsilon_2+d\varepsilon_1\varepsilon_2) = f(a) + b f'(a) \varepsilon_1 + c f'(a) \varepsilon_2 + bc f''(a) \varepsilon_1 \varepsilon_2 $

## Examples

$ f(x) = 2x^2-3x+12 $ therefore $ f'(x) = 4x-3 $ and $ f''(x) = 4 $.

$ f(7) = 2*7^2-3*7+12 = 89 $ and $ f'(7) = 4*7-3 = 25 $ and $ f''(7) = 4 $.

```python
> f = lambda x : 2 * x**2 - 3*x + 12
> f(HyperDual(7, 1, 1))
HyperDual(89, 25, 25, 4)
```