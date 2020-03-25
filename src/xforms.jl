#
# xforms.jl --
#
# Implementation of affine transforms which are notably useful for coordinate
# transforms.
#
#------------------------------------------------------------------------------
#
# This file if part of the TwoDimensional Julia package licensed under the MIT
# license (https://github.com/emmt/TwoDimensional.jl).
#
# Copyright (C) 2016-2019, Éric Thiébaut.
#

module AffineTransforms

export
    AffineTransform,
    compose,
    intercept,
    jacobian,
    rotate,
    scale,
    translate

using ..TwoDimensional: AbstractPoint, Point

# Imports for extension.
import Base: +, -, *, ∘, /, \, inv, eltype
import Base: Float16, Float32, Float64
import Base.MPFR: BigFloat
import LinearAlgebra: ⋅, det

"""
# Affine 2D Transforms

An affine 2D transform `A` is defined by 6 real coefficients, `Axx`, `Axy`,
`Ax`, `Ayx`, `Ayy` and `Ay`.  Such a transform maps `(x,y)` as `(xp,yp)` given
by:

```julia
xp = Axx*x + Axy*y + Ax
yp = Ayx*x + Ayy*y + Ay
```

The immutable type `AffineTransform` is used to store an affine 2D transform
`A`, it can be created by:

```julia
I = AffineTransform{T}() # yields the identity with type T
A = AffineTransform{T}(Axx, Axy, Ax, Ayx, Ayy, Ay)
```

The parameter `T` above is used to specify the floating-point type for the
coefficients; if omitted it is guessed from the types of the coefficients.


## Operations with affine 2D transforms

Many operations are available to manage or apply affine transforms:

```julia
(xp, yp) = A(x,y)       # apply affine transform A to coordinates (x,y)
(xp, yp) = A*(x,y)      # idem
(xp, yp) = A(v)         # idem, with v = (x,y)
(xp, yp) = A*v          # idem

A(Point(x,y)) -> Point(xp, yp)
A*Point(x,y)  -> Point(xp, yp)

C = compose(A, B, ...)  # compose 2 (or more) transforms, C = apply B then A
C = A∘B                 # idem
C = A*B                 # idem
C = A⋅B                 # idem

B = translate(x, y, A)  # B = apply A then translate by (x,y)
B = translate(v, A)     # idem with v = (x,y)
B = v + A               # idem

B = translate(A, x, y)  # B = translate by (x,y) then apply A
B = translate(A, v)     # idem with v = (x,y)
B = A + v               # idem

B = rotate(θ, A)   # B = apply A then rotate by angle θ
C = rotate(A, θ)   # C = rotate by angle θ then apply A

B = scale(ρ, A)    # B = apply A then scale by ρ
B = ρ*A            # idem
C = scale(A, ρ)    # C = scale by ρ then apply A
C = A*ρ            # idem

B = inv(A)         # reciprocal coordinate transform
C = A/B            # right division, same as: C = compose(A, inv(B))
C = A\\B            # left division, same as: C = compose(inv(A), B)
```

"`∘`" and "`⋅`" can be typed by `\\circ<tab>` and `\\cdot<tab>`.


## Type conversion

As a general rule, the floating-point type `T` of an `AffineTransform{T}` is
imposed for all operations and for the result.  The floating-point type of the
composition of several coordinate transforms is the promoted type of the
transforms which have been composed.

Calling `eltype(A)` yields floating-point type of the coefficients of the 2D
affine transform `A`.  To convert the floating-point type of the coefficients
of `A` to be `T`, do one of:

```julia
B = T.(A)
B = AffineTransform{T}(A)
B = convert(AffineTransform{T}, A)
```

"""
struct AffineTransform{T<:AbstractFloat} <: Function
    xx::T
    xy::T
    x ::T
    yx::T
    yy::T
    y ::T
    AffineTransform{T}() where T = new{T}(1,0,0, 0,1,0)
    function AffineTransform{T}(Axx::Real, Axy::Real, Ax::Real,
                                Ayx::Real, Ayy::Real, Ay::Real) where T
        new{T}(Axx,Axy,Ax, Ayx,Ayy,Ay)
    end
end

eltype(::AffineTransform{T}) where {T} = T

AffineTransform() = AffineTransform{Float64}()
function AffineTransform(Axx::Txx, Axy::Txy, Ax::Tx,
                         Ayx::Tyx, Ayy::Tyy, Ay::Ty
                         ) where {Txx<:Real,Txy<:Real,Tx<:Real,
                                  Tyx<:Real,Tyy<:Real,Ty<:Real}
    T = float(promote_type(Txx, Txy, Tx, Tyx, Tyy, Ty))
    AffineTransform{T}(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

function AffineTransform(Axx::T, Axy::T, Ax::T,
                         Ayx::T, Ayy::T, Ay::T) where {T<:AbstractFloat}
    AffineTransform{T}(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

AffineTransform(A::AffineTransform) = A
AffineTransform{T}(A::AffineTransform{T}) where {T<:AbstractFloat} = A
AffineTransform{T}(A::AffineTransform) where {T<:AbstractFloat} =
    AffineTransform{T}(A.xx, A.xy, A.x,
                       A.yx, A.yy, A.y)

# Allow for `T.(obj)` to work with `T` a floating-point type.
Broadcast.broadcasted(::Type{T}, A::AffineTransform) where {T<:AbstractFloat} =
    AffineTransform{T}(A)

Base.convert(::Type{T}, A::AffineTransform) where {T<:AffineTransform} = T(A)
Base.promote_type(::Type{AffineTransform{T}}, ::Type{AffineTransform{U}}) where {T,U} =
    AffineTransform{promote_type(T,U)}
Base.promote_type(::Type{AffineTransform{T}}, ::Type{AffineTransform{T}}) where {T} =
    AffineTransform{T}

@deprecate(Float16(A::AffineTransform),  AffineTransform{Float16}(A))
@deprecate(Float32(A::AffineTransform),  AffineTransform{Float32}(A))
@deprecate(Float64(A::AffineTransform),  AffineTransform{Float64}(A))
@deprecate(BigFloat(A::AffineTransform), AffineTransform{BigFloat}(A))

#------------------------------------------------------------------------------
# apply the transform to some coordinates:

(A::AffineTransform{T})(x::T, y::T) where {T<:AbstractFloat} =
    (A.xx*x + A.xy*y + A.x,
     A.yx*x + A.yy*y + A.y)

(A::AffineTransform{T})(x::Real, y::Real) where {T<:AbstractFloat} =
    A(convert(T, x), convert(T, y))

(A::AffineTransform)(v::Tuple{Real,Real}) = A(v[1], v[2])

(A::AffineTransform)(v::Point) = Point(A(v.x, v.y)...)

#------------------------------------------------------------------------------
# Combine a translation with an affine transform.

"""
### Translating an affine transform

Affine transforms can be letf- or right-translated.

```julia
translate(x, y, A)
```
or
```julia
translate((x,y), A)
```

yield an affine transform which translate the output of affine transform `A` by
offsets `x` and `y`.

```julia
translate(A, x, y)
```
or
```julia
translate(A, (x,y))
```

yield an affine transform which translate the input of affine transform `A` by
offsets `x` and `y`.

The same results can be obtained with the `+` operator:

```julia
B = (x,y) + A    # same as: B = translate((x,y), A)
B = A + (x,y)    # same as: B = translate(A, (x,y))
```

See also: [`AffineTransform`](@ref), [`rotate`](@ref), [`scale`](@ref).

""" translate

# Left-translating results in translating the output of the transform.
translate(x::T, y::T, A::AffineTransform{T}) where {T<:AbstractFloat} =
    AffineTransform{T}(A.xx, A.xy, A.x + x,
                       A.yx, A.yy, A.y + y)

translate(x::Real, y::Real, A::AffineTransform{T}) where {T<:AbstractFloat} =
    translate(convert(T, x), convert(T, y), A)

translate(v::Tuple{Real,Real}, A::AffineTransform) =
    translate(v[1], v[2], A)

translate(v::AbstractPoint, A::AffineTransform) =
    translate(v.x, v.y, A)

# Right-translating results in translating the input of the transform.
translate(A::AffineTransform{T}, x::T, y::T) where {T<:AbstractFloat} =
    AffineTransform{T}(A.xx, A.xy, A.xx*x + A.xy*y + A.x,
                       A.yx, A.yy, A.yx*x + A.yy*y + A.y)

translate(A::AffineTransform{T}, x::Real, y::Real) where {T<:AbstractFloat} =
    translate(A, convert(T, x), convert(T, y))

translate(A::AffineTransform, v::Tuple{Real,Real}) =
    translate(A, v[1], v[2])

translate(A::AffineTransform, v::AbstractPoint) =
    translate(A, v.x, v.y)

#------------------------------------------------------------------------------
"""
### Scaling an affine transform

There are two ways to combine a scaling by a factor `ρ` with an affine
transform `A`.  Left-scaling as in:

```julia
B = scale(ρ, A)
```

results in scaling the output of the transform; while right-scaling as in:

```julia
C = scale(A, ρ)
```

results in scaling the input of the transform.  The above examples yield
transforms which behave as:

```julia
B(v) = ρ.*A(v)
C(v) = A(ρ.*v)
```

where `v` is any 2-element tuple.

The same results can be obtained with the `*` operator:

```julia
B = ρ*A    # same as: B = scale(ρ, A)
C = A*ρ    # same as: B = scale(A, ρ)
```

See also: [`AffineTransform`](@ref), [`rotate`](@ref), [`translate`](@ref).

"""
scale(ρ::T, A::AffineTransform{T}) where {T<:AbstractFloat} =
    AffineTransform{T}(ρ*A.xx, ρ*A.xy, ρ*A.x,
                       ρ*A.yx, ρ*A.yy, ρ*A.y)

scale(A::AffineTransform{T}, ρ::T) where {T<:AbstractFloat} =
    AffineTransform{T}(ρ*A.xx, ρ*A.xy, A.x,
                       ρ*A.yx, ρ*A.yy, A.y)

#------------------------------------------------------------------------------
"""
### Rotating an affine transform

There are two ways to combine a rotation by angle `θ` (in radians
counterclockwise) with an affine transform `A`.  Left-rotating as in:

```julia
B = rotate(θ, A)
```

results in rotating the output of the transform; while right-rotating as in:

```julia
C = rotate(A, θ)
```

results in rotating the input of the transform.  The above examples are
similar to:

```julia
B = R∘A
C = A∘R
```

where `R` implements rotation by angle `θ` around `(0,0)`.


See also: [`AffineTransform`](@ref), [`scale`](@ref), [`translate`](@ref).

"""
function rotate(θ::T, A::AffineTransform{T}) where {T<:AbstractFloat}
    cs = cos(θ)
    sn = sin(θ)
    return AffineTransform{T}(cs*A.xx - sn*A.yx,
                              cs*A.xy - sn*A.yy,
                              cs*A.x  - sn*A.y,
                              cs*A.yx + sn*A.xx,
                              cs*A.yy + sn*A.xy,
                              cs*A.y  + sn*A.x)
end

function rotate(A::AffineTransform{T}, θ::T) where {T<:AbstractFloat}
    cs = cos(θ)
    sn = sin(θ)
    return AffineTransform{T}(A.xx*cs + A.xy*sn,
                              A.xy*cs - A.xx*sn,
                              A.x,
                              A.yx*cs + A.yy*sn,
                              A.yy*cs - A.yx*sn,
                              A.y)
end

# Make sure the floating-point type of an affine transform is preserved.
for func in (:scale, :rotate)
    @eval begin
        $func(α::Real, A::AffineTransform{T}) where {T<:AbstractFloat} =
            $func(convert(T, α), A)
        $func(A::AffineTransform{T}, α::Real) where {T<:AbstractFloat} =
            $func(A, convert(T, α))
    end
end

#------------------------------------------------------------------------------

"""
`det(A)` returns the determinant of the linear part of the affine
transform `A`.
"""
det(A::AffineTransform) = A.xx*A.yy - A.xy*A.yx

"""
`jacobian(A)` returns the Jacobian of the affine transform `A`, that is the
absolute value of the determinant of its linear part.
"""
jacobian(A::AffineTransform) = abs(det(A))

"""
`inv(A)` returns the inverse of the affine transform `A`.
"""
function inv(A::AffineTransform{T}) where {T}
    Δ = det(A)
    Δ == zero(T) && error("transformation is not invertible")
    α = one(T)/Δ
    Txx =  α*A.yy
    Txy = α*A.xy
    Tyx = α*A.yx
    Tyy =  α*A.xx
    return AffineTransform{T}(+Txx, -Txy, Txy*A.y - Txx*A.x,
                              -Tyx, +Tyy, Tyx*A.x - Tyy*A.y)
end

"""

`compose(A,B)` yields the affine transform which combines the two affine
transforms `A` and `B`, that is the affine transform which applies `B` and then
`A`.  Composition is accessible via: `A∘B`, `A*B` or `A⋅B` ("`∘`" and "`⋅`" can
be typed by `\\circ<tab>` and `\\cdot<tab>`).

It is possible to compose more than two affine transforms.  For instance,
`compose(A,B,C)` yields the affine transform which applies `C` then `B`, then
`A`.

"""
compose(A::AffineTransform) = A

compose(A::AffineTransform{T}, B::AffineTransform{T}) where {T<:AbstractFloat} =
    AffineTransform{T}(A.xx*B.xx + A.xy*B.yx,
                       A.xx*B.xy + A.xy*B.yy,
                       A.xx*B.x  + A.xy*B.y + A.x,
                       A.yx*B.xx + A.yy*B.yx,
                       A.yx*B.xy + A.yy*B.yy,
                       A.yx*B.x  + A.yy*B.y + A.y)

compose(A::AffineTransform, B::AffineTransform) =
    compose(promote(A, B)...)

compose(A::AffineTransform, B::AffineTransform, C::AffineTransform, args::AffineTransform...) =
    compose(compose(A, B), C, args...)

"""

`rightdivide(A,B)` yields `A/B`, the right division of the affine
transform `A` by the affine transform `B`.

"""
rightdivide(A::AffineTransform{T}, B::AffineTransform{T}) where {T<:AbstractFloat} = begin
    Δ = det(B)
    Δ == zero(T) && error("right operand is not invertible")
    α = one(T)/Δ
    Rxx = α*(A.xx*B.yy - A.xy*B.yx)
    Rxy = α*(A.xy*B.xx - A.xx*B.xy)
    Ryx = α*(A.yx*B.yy - A.yy*B.yx)
    Ryy = α*(A.yy*B.xx - A.yx*B.xy)
    AffineTransform{T}(Rxx, Rxy, A.x - (Rxx*B.x + Rxy*B.y),
                       Ryx, Ryy, A.y - (Ryx*B.y + Ryy*B.y))

end

rightdivide(A::AffineTransform, B::AffineTransform) =
    rightdivide(promote(A, B)...)

"""
`leftdivide(A,B)` yields `A\\B`, the left division of the affine
transform `A` by the affine transform `B`.
"""
leftdivide(A::AffineTransform{T}, B::AffineTransform{T}) where {T<:AbstractFloat} = begin
    Δ = det(B)
    Δ == zero(T) && error("left operand is not invertible")
    α = one(T)/Δ
    Txx = α*A.yy
    Txy = α*A.xy
    Tyx = α*A.yx
    Tyy = α*A.xx
    Tx = B.x - A.x
    Ty = B.y - A.y
    AffineTransform{T}(Txx*B.xx - Txy*B.yx,
                       Txx*B.xy - Txy*B.yy,
                       Txx*Tx   - Txy*Ty,
                       Tyy*B.yx - Tyx*B.xx,
                       Tyy*B.yy - Tyx*B.xy,
                       Tyy*Ty   - Tyx*Tx)
end

leftdivide(A::AffineTransform, B::AffineTransform) =
    leftdivide(promote(A, B)...)

"""

`intercept(A)` returns the tuple `(x,y)` such that `A(x,y) = (0,0)`.

"""
function intercept(A::AffineTransform{T}) where {T}
    Δ = det(A)
    Δ == zero(T) && error("transformation is not invertible")
    α = one(T)/Δ
    return (α*(A.xy*A.y - A.yy*A.x), α*(A.yx*A.x - A.xx*A.y))
end

intercept(T::Type{<:Point}, A::AffineTransform) = T(intercept(A)...)

+(v::Union{AbstractPoint,Tuple{Real,Real}}, A::AffineTransform) =
    translate(v, A)

+(A::AffineTransform, v::Union{AbstractPoint,Tuple{Real,Real}}) =
    translate(A, v)

-(A::AffineTransform, v::Tuple{Real,Real}) = A + (-v[1], -v[2])
-(A::AffineTransform, v::AbstractPoint) = A + (-v.x, -v.y)

for op in (:∘, :*, :⋅)
    @eval begin
        $op(A::AffineTransform, B::AffineTransform) = compose(A, B)
    end
end

*(A::AffineTransform, v::Union{Point,Tuple{Real,Real}}) = A(v)

*(ρ::Real, A::AffineTransform) = scale(ρ, A)

*(A::AffineTransform, ρ::Real) = scale(A, ρ)

\(A::AffineTransform, B::AffineTransform) = leftdivide(A, B)

/(A::AffineTransform, B::AffineTransform) = rightdivide(A, B)

Base.show(io::IO, ::MIME"text/plain", A::AffineTransform) =
    print(io, typeof(A),
          "(",   A.xx, ",", A.xy, ",", A.x,
          ",  ", A.yx, ",", A.yy, ",", A.y, ")")

Base.show(io::IO, A::AffineTransform) = show(io, MIME"text/plain"(), A)

Base.print(io::IOBuffer, A::AffineTransform) = show(io, A)

end # module
