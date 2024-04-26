# Affine Coordinate Transforms

`TwoDimensional` provides types and methods to deal with 2-dimensional affine
coordinate transforms.

An affine 2D transform `A` is defined by 6 real coefficients, `Axx`, `Axy`,
`Ax`, `Ayx`, `Ayy` and `Ay`. Such a transform maps `(x,y)` as `(xp,yp)` given
by:

```julia
xp = Axx*x + Axy*y + Ax
yp = Ayx*x + Ayy*y + Ay
```

coefficients `Axx`, `Axy`, `Ayx`, and `Ayy` are *factors* while coefficients `Ax`
and `Ay` are *offsets*.

The immutable type `AffineTransform{T}` is used to store an affine 2D transform
with coefficients of floating-point type `T`, it can be created by:

```julia
I = AffineTransform{T}() # yields the identity with type T
A = AffineTransform{T}(Axx, Axy, Ax, Ayx, Ayy, Ay)
```

If the parameter `T` is omitted, it is guessed from the types of the
coefficients.


## Operations with affine 2D transforms

Many operations are available to manage or apply affine transforms.


### Perform a coordinate transformation

There are several possibilities to apply an affine transform `A` to coordinates
`(x,y)`:

```julia
(xp, yp) = A(x,y)       # apply affine transform A to coordinates (x,y)
(xp, yp) = A*(x,y)      # idem
(xp, yp) = A(v)         # idem, with v = (x,y)
(xp, yp) = A*v          # idem
```

where `x`, `y`, `xp` and `yp` are reals while `v = (x,y)` is a 2-tuple of real
coordinates.  Coordinates may also be specified as `Point` instances:

```julia
A(Point(x,y)) -> Point(xp, yp)
A*Point(x,y)  -> Point(xp, yp)
```


### Reciprocal transform

The reciprocal transform of `A` is given by:

```julia
B = inv(A)
```

### Compose coordinate transformations

To compose 2 (or more) transforms `A` and `B`, do one of:

```julia
C = compose(A, B, ...)
C = A∘B
C = A*B
C = A⋅B
```

all these statements yields an object `C` which applies `B` then `A`.  Note
that `∘` and `⋅` can be typed by `\\circ<tab>` and `\\cdot<tab>`.

Left and right *"divisions"* of affine tansforms respectively write:

```julia
A/B -> A∘inv(B)
A\B -> inv(A)∘B
```


### Translation of coordinates

An operator `B` which which applies a translations by `(x,y)` *after* a
coordinate transform `A` (possibly the identity) can be created by one of the
following statements:

```julia
B = translate(x, y, A)
B = translate(pnt, A)
B = pnt + A
```

where `pnt` is a 2-tuple of coordinates, `pnt = (x,y)`, a `Point`, `pnt =
Point(x,y)`, or 2-dimensional Cartesian index, `CartesianIndex{2}(x,y)`.

To perform the translation *before* the coordinate transform `A`, do:

```julia
B = translate(A, x, y)
B = translate(A, pnt)
B = A + pnt
```

### Rotating an affine transform

There are two ways to combine a rotation by angle `θ` (in radians
counterclockwise) with an affine transform `A`. Left-rotating as in:

```julia
B = rotate(θ, A)
```

results in rotating the output of the transform; while right-rotating as in:

```julia
C = rotate(A, θ)
```

results in rotating the input of the transform. The above examples are similar
to:

```julia
B = R∘A
C = A∘R
```

where `R` implements rotation by angle `θ` counterclockwise around the origin
at coordinates `(0,0)`. The rotation angle `θ` is assumed to be in radians if
it has no units.


### Scaling of coordinates

There are two ways to combine a scaling by a factor `ρ` with an affine
transform `A`. Left-scaling as in:

```julia
B = scale(ρ, A)
B = ρ*A
```

which yields `B` such that `B(x,y) = ρ*A(x,y)`; or right-scaling as in:

```julia
C = scale(A, ρ)
C = A*ρ
```

which yields `C` such that `C(x,y) = A(ρ*x,ρ*y)`; or right-scaling as in:


### Miscellaneous

`det(A)` returns the determinant of the linear part of the affine transform
`A`.

`jacobian(A)` returns the Jacobian of the affine transform `A`, that is the
absolute value of the determinant of its linear part.


## Type conversion

As a general rule, the floating-point type `T` of an `AffineTransform{T}` is
imposed for all operations and for the result.  The floating-point type of the
composition of several coordinate transforms is the promoted type of the
transforms which have been composed.

The type of the coefficients of the affine transform `A`  is given by:

```julia
eltype(A)
```

To convert the floating-point type of the coefficients of `A` to be `T`, do one
of:

```julia
B = T.(A)
B = AffineTransform{T}(A)
B = convert(AffineTransform{T}, A)
```
