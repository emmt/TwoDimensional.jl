# Affine Coordinate Transforms

TwoDimensional provides types and methods to deal with 2-dimensional affine
coordinate transforms.

An affine 2D transform `C` is defined by 6 real coefficients, `Cxx`, `Cxy`,
`Cx`, `Cyx`, `Cyy` and `Cy`.  Such a transform maps `(x,y)` as `(xp,yp)` given
by:

```julia
xp = Cxx*x + Cxy*y + Cx
yp = Cyx*x + Cyy*y + Cy
```

The immutable type `AffineTransform{T}` is used to store an affine 2D transform
with coefficients of floating-point type `T`, it can be created by:

```julia
I = AffineTransform{T}() # yields the identity with type T
C = AffineTransform{T}(Cxx, Cxy, Cx, Cyx, Cyy, Cy)
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

Left and right divisions respectively write:

```julia
R = A/B   # right division, same as: R = compose(A, inv(B))
L = A\B   # left division, same as: L = compose(inv(A), B)
```


### Translation of coordinates

An operator `B` which which applies a translations by `(x,y)` *after* a
coordinate transform `A` (possibly the identity) can be created by one of the
following statements:

```julia
B = translate(x, y, A)
B = translate(v, A)
B = v + A
```

where `v` is a 2-tuple of coordinates, *i.e.* `v = (x,y)` or a `Point`, *i.e.*
`v = Point(x,y)`.

To perform the translation *before* the coordinate transform `A`, do:

```julia
B = translate(A, x, y)
B = translate(A, v)
B = A + v
```

### Rotation of coordinates

```julia
B = rotate(θ, A)   # B = apply A then rotate by angle θ
C = rotate(A, θ)   # C = rotate by angle θ then apply A
```

### Scaling of coordinates

```julia
B = scale(ρ, A)    # B = apply A then scale by ρ
B = ρ*A            # idem
C = scale(A, ρ)    # C = scale by ρ then apply A
C = A*ρ            # idem
```

### Scaling of coordinates


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
