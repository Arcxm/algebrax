# Common Interface

This document describes the common interface all classes implement.

## TOC

- [Readability](#readability)
- [Comparison](#comparison)
- [Arithmetic](#arithmetic)
- [Others](#others)
- [Functions](#functions)

## Readability

Both `__str__` and `__repr__` return the number in a readable tuple form.
Note that `__repr__` adds the name of the class in front.

## Comparison

`__eq__` and `__ne__` allow comparison between numbers.
An object of class C can be compared with:

- an int or float (real number)
- a class D that is a component of C ("real" component, eg. `Quat` in `DualQuat`)
- another object of class C

Note that `SplitComplex` and `SplitComplexIdempotent` allow mixed comparison (convert under the hood).

## Arithmetic

`__neg__`, `__add__`, `__sub__`, `__mul__`, `__truediv__` and `__pow__` as well as their right counterparts (e.g. `__rsub__`) behave as expected.

Again `SplitComplex` and `SplitComplexIdempotent` allow mixed arithmetic, returning an object of the type of the left-hand operand.

Please note that non-commutative algebras do not implement `__truediv__` as well as `__rtruediv__` so you have to use multiplication with the `inverse` for clarity (might change in the future).

## Others

### conjugate

Returns the conjugate of a number, swapping the sign of the non-real components.

### norm

Returns the norm of a number (`__abs__` is implemented and returns `norm`).

Note that the norm is not necessarily Euclidean. Split-algebras use an indefinite norm which may be zero or negative.

### inverse

Returns the inverse of a number.

## Functions

Note that if an operation is mathematically undefined, e.g. outside of its domain, a `ValueError` is raised containing a short note on the domain.

### exp

The exponential function, raises $e$ to the power of a number.

### log

The logarithm, the inverse of the exponential function.

### sqrt

The square root of a number.

### root

Computes the n-th root of a number.

`n` should be a positive integer greater than $0$.