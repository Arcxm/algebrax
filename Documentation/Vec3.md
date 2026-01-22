# Vec3

The `Vec3` class represents a vector in $ \mathbb{R}^3 $.

A 3D vector consists of 3 components: $x$, $y$ and $z$.

## Custom functions

- **normalize:** normalizes the vector to a length of 1
- **dot:** calculates the dot product of 2 vectors
- **cross:** calculates the cross product of 2 vectors

## Examples

```python
> v1 = Vec3(1, 2, 3)
> v2 = v1.normalize()
> v2
Vec3(0.2672612419124244, 0.5345224838248488, 0.8017837257372732)
> v2.norm()
1.0
```

## Further information

[Wikipedia: Euclidean vector](https://en.wikipedia.org/wiki/Euclidean_vector).