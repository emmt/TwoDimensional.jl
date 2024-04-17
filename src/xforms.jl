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
# Copyright (c) 2016-2024, Éric Thiébaut.
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

using Unitless

using ..TwoDimensional: AbstractPoint, Point, PointLike, get_x, get_y

# Imports for extension.
import Base: +, -, *, ∘, /, \, inv
import Base: Float16, Float32, Float64
import Base.MPFR: BigFloat
import LinearAlgebra: ⋅, det

"""
    AffineTransform(Axx, Axy, Ax, Ayx, Ayy, Ay) -> A

yields a callable object implementing a 2-dimensional affine transform such
that:

    A(x, y) -> (Axx*x + Axy*y + Ax, Ayx*x + Ayy*y + Ay)
    A((x, y)) -> (Axx*x + Axy*y + Ax, Ayx*x + Ayy*y + Ay)
    A(pnt::Point) -> Point(A(pnt.x, pnt.y))
    A(pnt::CartesianIndex{2}) -> Point(A(pnt[1], pnt[2]))

The constructor takes one or three parameters:

    AffineTransform{T}(Axx, Axy, Ax, Ayx, Ayy, Ay)
    AffineTransform{T,R,S}(Axx, Axy, Ax, Ayx, Ayy, Ay)

where `T` is the concrete floating-point type of the coefficients (`Float64` by
default), `R` is the type for storing the factors `Axx`, `Axy`, `Ayx`, and
`Ayy`, and `S` is the type for storing the offsets `Ax` and `Ay`. The bare type
of `R` and `S` must be `T` but they may have units.

Changing the floating-point type of an existing 2-dimensional affine transform
`A` can be done by one of:

    B = AffineTransform{T}(A)
    B = convert(AffineTransform{T}, A)
    B = convert_bare_type(T, A)
    B = convert_real_type(T, A)
    B = convert_floating_point_type(T, A)

the 3 last assume `using Unitless`. Using `Unitless` API, the floating-point
type `T` can be retrieved by either of `bare_type(A)`, real_type(A)`, or
`floating_point_type(A)` with `A` a 2-dimensional affine transform instance or
type.

Applying the 2-dimensional affine transform `A` can be done by:

```julia
(xp, yp) = A(x,y)       # apply affine transform A to coordinates (x,y)
(xp, yp) = A*(x,y)      # idem
(xp, yp) = A((x,y))     # idem, with pnt = (x,y)
(xp, yp) = A*(x,y)      # idem

A(Point(x,y)) -> Point(xp, yp)
A*Point(x,y)  -> Point(xp, yp)
A⋅Point(x,y)  -> Point(xp, yp)

C = compose(A, B, ...)  # compose 2 (or more) transforms, C = apply B then A
C = A∘B                 # idem
C = A*B                 # idem
C = A⋅B                 # idem

B = translate(x, y, A)  # B = apply A then translate by (x,y)
B = translate(pnt, A)   # idem with pnt = (x,y)
B = pnt + A             # idem

B = translate(A, x, y)  # B = translate by (x,y) then apply A
B = translate(A, pnt)   # idem with pnt = (x,y)
B = A + pnt             # idem

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

"""
struct AffineTransform{T<:AbstractFloat,R,S} <: Function
    xx::R
    xy::R
    x ::S
    yx::R
    yy::R
    y ::S
    function AffineTransform{T,R,S}(Axx, Axy, Ax, Ayx, Ayy, Ay) where {T<:AbstractFloat,R,S}
        isconcretetype(T) || throw(ArgumentError(
            "type parameter `T = $T` is not a concrete floating-point type"))
        bare_type(R) === T || throw(ArgumentError(
            "bare type of parameter `R = $R` is not `T = $T`, got `bare_type(R) = $(bare_type(R))`"))
        bare_type(S) === T || throw(ArgumentError(
            "bare type of parameter `S = $S` is not `T = $T`, got `bare_type(S) = $(bare_type(S))`"))
        return new{T,R,S}(Axx, Axy, Ax, Ayx, Ayy, Ay)
    end
end

"""
    AffineTransform() -> Id
    AffineTransform{T}() -> Id
    AffineTransform{T,R,S}() -> Id

yields a 2-dimensional affine transform corresponding to the identity (up to
possible change of type and units). Parameter `T` is the floating-point type of
the coefficients (`Float64` by default). Parameters `R` and `S` are the types
of the factors and of the offsets (by default both are assumed to be `T`).
These are shortcuts to:

     AffineTransform(oneunit(R),zero(R),zero(S), zero(R),oneunit(R),zero(S))

"""
AffineTransform() = AffineTransform{Float64}()
AffineTransform{T}() where {T<:AbstractFloat} = AffineTransform{T,T,T}()
AffineTransform{T,R,S}() where {T<:AbstractFloat,R,S} =
    AffineTransform{T,R,S}(oneunit(R),zero(R),zero(S),
                           zero(R),oneunit(R),zero(S))

function AffineTransform(Axx, Axy, Ax,
                         Ayx, Ayy, Ay)
    Axx, Axy, Ayx, Ayy = promote(Axx, Axy, Ayx, Ayy)
    Ax, Ay = promote(Ax, Ay)
    return AffineTransform(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

function AffineTransform{T}(Axx, Axy, Ax,
                            Ayx, Ayy, Ay) where {T<:AbstractFloat}
    Axx, Axy, Ayx, Ayy = promote(Axx, Axy, Ayx, Ayy)
    Ax, Ay = promote(Ax, Ay)
    return AffineTransform{T}(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

function AffineTransform(Axx::R, Axy::R, Ax::S,
                         Ayx::R, Ayy::R, Ay::S) where {R,S}
    T = float(promote_type(real_type(R), real_type(S)))
    return AffineTransform{T}(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

function AffineTransform{T}(Axx::R, Axy::R, Ax::S,
                            Ayx::R, Ayy::R, Ay::S) where {T<:AbstractFloat,R,S}
    Rc = convert_real_type(T, R)
    Sc = convert_real_type(T, S)
    return AffineTransform{T,Rc,Sc}(Axx,Axy,Ax, Ayx,Ayy,Ay)
end

AffineTransform(A::AffineTransform) = A
AffineTransform{T}(A::AffineTransform{T}) where {T<:AbstractFloat} = A
AffineTransform{T}(A::AffineTransform) where {T<:AbstractFloat} =
    AffineTransform{T}(Tuple(A)...)

# Extend Unitless.bare_type, Unitless.real_type, Unitless.floating_point_type,
# Unitless.convert_bare_type, Unitless.convert_real_type, and
# Unitless.convert_floating_point_type,
for func in (:bare_type, :real_type, :floating_point_type)
    conv = Symbol("convert_",func)
    @eval begin
        Unitless.$func(A::AffineTransform) = $func(typeof(A))
        Unitless.$func(::Type{<:AffineTransform{T}}) where {T} = T
        Unitless.$conv(::Type{T}, A::AffineTransform{T}) where {T} = A
        Unitless.$conv(::Type{T}, A::AffineTransform) where {T<:AbstractFloat} =
            AffineTransform{T}(A)
    end
end

"""
    TwoDimensional.factors_type(A) -> R

yields the type of the factors of the 2-dimensional affine transform `A`. The
factors of `A` are the coefficients `A.xx`, `A,xy`, `A.yx`, and `A.yy`.
Argument may also be the type of an affine transform.

"""
factors_type(A::AffineTransform) = factors_type(typeof(A))
factors_type(::Type{AffineTransform{T,R,S}}) where {T,R,S} = R

"""
    TwoDimensional.offsets_type(A) -> R

yields the type of the offsets of the 2-dimensional affine transform `A`. The
offsets of `A` are the coefficients `A.x` and `A.y`. Argument may also be the
type of an affine transform.

"""
offsets_type(A::AffineTransform) = offsets_type(typeof(A))
offsets_type(::Type{AffineTransform{T,R,S}}) where {T,R,S} = S

Base.convert(::Type{T}, A::AffineTransform) where {T<:AffineTransform} = T(A)
Base.promote_type(::Type{AffineTransform{T}}, ::Type{AffineTransform{U}}) where {T,U} =
    AffineTransform{promote_type(T,U)}
Base.promote_type(::Type{AffineTransform{T}}, ::Type{AffineTransform{T}}) where {T} =
    AffineTransform{T}

# Make affine transform objects indexable and iterable.
Base.Tuple(A::AffineTransform) = (A.xx, A.xy, A.x, A.yx, A.yy, A.y)

Base.getindex(A::AffineTransform, i::Integer) =
    i == 1 ? A.xx :
    i == 2 ? A.xy :
    i == 3 ? A.x  :
    i == 4 ? A.yx :
    i == 5 ? A.yy :
    i == 6 ? A.y  : throw(BoundsError(A, i))

Base.iterate(A::AffineTransform, i::Int=1) =
    i == 1 ? (A.xx, 2) :
    i == 2 ? (A.xy, 3) :
    i == 3 ? (A.x,  4) :
    i == 4 ? (A.yx, 5) :
    i == 5 ? (A.yy, 6) :
    i == 6 ? (A.y,  7) : nothing

Base.IteratorSize(A::AffineTransform) = Base.IteratorSize(typeof(A))
Base.IteratorSize(::Type{<:AffineTransform}) = Base.HasLength()
Base.length(A::AffineTransform) = 6

#------------------------------------------------------------------------------
# Apply the transform to some coordinates (promoted to the same type).
(A::AffineTransform)(pnt::Union{Point,CartesianIndex{2}}) = Point(A(get_x(pnt), get_y(pnt)))
(A::AffineTransform)((x, y)::Tuple{Any,Any}) = A(x, y)
(A::AffineTransform)(x, y) = A(promote(x, y)...)
(A::AffineTransform)(x::T, y::T) where {T} = (A.xx*x + A.xy*y + A.x,
                                              A.yx*x + A.yy*y + A.y)

*(A::AffineTransform, pnt::PointLike) = A(pnt)
⋅(A::AffineTransform, pnt::PointLike) = A(pnt)

#------------------------------------------------------------------------------
# Combine a translation with an affine transform.

"""
    B = translate(x, y, A)
    B = translate((x,y), A)
    B = translate(pnt, A)
    B = (x,y) + A
    B = pnt + A

perform a left-translation of the 2-dimensional affine transform `A`. Applying
`B` yields the same result as if coordinates `(x,y)` are added to the output of
`A`. Here, `pnt = Point(x,y)` or `pnt = CartesianIndex(x,y)`.

    C = translate(A, x, y)
    C = translate(A, (x,y))
    C = translate(A, pnt)
    C = A + (x,y)
    C = A + pnt

perform a right-translation of the 2-dimensional affine transform `A`. Applying
`B` yields the same result as if coordinates `(x,y)` are added to the input of
`A`.

See also: [`AffineTransform`](@ref), [`rotate`](@ref), [`scale`](@ref).

""" translate

# Left-translating results in translating the output of the transform.
+(pnt::PointLike, A::AffineTransform) = translate(pnt, A)
-(pnt::PointLike, A::AffineTransform) = pnt + (-A)
translate(pnt::PointLike, A::AffineTransform) = translate(get_x(pnt), get_y(pnt), A)
translate(x, y, A::AffineTransform) = AffineTransform(A.xx, A.xy, A.x + x,
                                                      A.yx, A.yy, A.y + y)

# Right-translating results in translating the input of the transform.
+(A::AffineTransform, pnt::PointLike) = translate(A, pnt)
-(A::AffineTransform, pnt::PointLike) = A + (-get_x(pnt), -get_y(pnt))
translate(A::AffineTransform, pnt::PointLike) = translate(A, get_x(pnt), get_y(pnt))
translate(A::AffineTransform, x, y) = AffineTransform(A.xx, A.xy, A.xx*x + A.xy*y + A.x,
                                                      A.yx, A.yy, A.yx*x + A.yy*y + A.y)

"""
    B = scale(ρ, A)
    B = ρ*A
    C = scale(A, ρ)
    C = A*ρ

yield 2-dmensional affine transforms `B` and `C` such that:

    B(x,y) -> ρ*A(x,y)
    C(x,y) -> A(ρ*x,ρ*y)

See also: [`AffineTransform`](@ref), [`rotate`](@ref), [`translate`](@ref).

""" scale

*(ρ::Number, A::AffineTransform) = scale(ρ, A)
scale(ρ, A::AffineTransform) = AffineTransform(ρ*A.xx, ρ*A.xy, ρ*A.x,
                                               ρ*A.yx, ρ*A.yy, ρ*A.y)

*(A::AffineTransform, ρ::Number) = scale(A, ρ)
scale(A::AffineTransform, ρ) = AffineTransform(ρ*A.xx, ρ*A.xy, A.x,
                                               ρ*A.yx, ρ*A.yy, A.y)

# Negating (unary minus) of a transform amounts to negating its output or to
# left-multiply the transform by -1. Hence, it is sufficient to negate all its
# coefficients.
-(A::AffineTransform) = AffineTransform(map(-, Tuple(A))...)

"""
    B = rotate(θ, A)
    C = rotate(A, θ)

yield 2-dmensional affine transforms `B` and `C` such that:

    B(x,y) = (R∘A)(x,y) = R(A(x,y))
    C(x,y) = (A∘R)(x,y) = A(R(x,y))

where `R` implements rotation by angle `θ` counterclockwise around the origin
at coordimates `(0,0)`. The rotation angle `θ` is assumed to be in radians if
it has no units.

See also: [`AffineTransform`](@ref), [`scale`](@ref), [`translate`](@ref).

"""
function rotate(θ, A::AffineTransform)
    sinθ, cosθ = sincos(θ)
    return AffineTransform(cosθ*A.xx - sinθ*A.yx,
                           cosθ*A.xy - sinθ*A.yy,
                           cosθ*A.x  - sinθ*A.y,
                           cosθ*A.yx + sinθ*A.xx,
                           cosθ*A.yy + sinθ*A.xy,
                           cosθ*A.y  + sinθ*A.x)
end

function rotate(A::AffineTransform, θ)
    sinθ, cosθ = sincos(θ)
    return AffineTransform(A.xx*cosθ + A.xy*sinθ,
                           A.xy*cosθ - A.xx*sinθ,
                           A.x,
                           A.yx*cosθ + A.yy*sinθ,
                           A.yy*cosθ - A.yx*sinθ,
                           A.y)
end

"""
    det(A::TwoDimensional.AffineTransform)

yields the determinant of the linear part of the affine transform `A`.

"""
det(A::AffineTransform) = A.xx*A.yy - A.xy*A.yx

"""
    jacobian(A::TwoDimensional.AffineTransform)

yields the Jacobian of the affine transform `A`, that is the absolute value of
the determinant of its linear part.

"""
jacobian(A::AffineTransform) = abs(det(A))

"""
    inv(A::TwoDimensional.AffineTransform)

yields the inverse of the affine transform `A`.

"""
function inv(A::AffineTransform)
    Δ = det(A)
    iszero(Δ) && error("transformation is not invertible")
    Rxx = A.yy/Δ
    Rxy = A.xy/Δ
    Ryx = A.yx/Δ
    Ryy = A.xx/Δ
    return AffineTransform(+Rxx, -Rxy, Rxy*A.y - Rxx*A.x,
                           -Ryx, +Ryy, Ryx*A.x - Ryy*A.y)
end

"""
    compose(A::TwoDimensional.AffineTransform, B::TwoDimensional.AffineTransform)

yields the affine transform which combines the two affine transforms `A` and
`B`, that is the affine transform which applies `B` and then `A`. Composition
is accessible via: `A∘B`, `A*B` or `A⋅B` ("`∘`" and "`⋅`" can be typed by
`\\circ<tab>` and `\\cdot<tab>`).

It is possible to compose more than two affine transforms. For instance,
`compose(A,B,C)` yields the affine transform which applies `C` then `B`, then
`A`.

"""
compose(A::AffineTransform) = A

compose(A::AffineTransform, B::AffineTransform) =
    AffineTransform(A.xx*B.xx + A.xy*B.yx,
                    A.xx*B.xy + A.xy*B.yy,
                    A.xx*B.x  + A.xy*B.y + A.x,
                    A.yx*B.xx + A.yy*B.yx,
                    A.yx*B.xy + A.yy*B.yy,
                    A.yx*B.x  + A.yy*B.y + A.y)

@inline compose(A::AffineTransform, B::AffineTransform, C::AffineTransform...) =
    compose(compose(A, B), C...)

for op in (:∘, :*, :⋅)
    @eval begin
        $op(A::AffineTransform, B::AffineTransform) = compose(A, B)
    end
end

"""
    rightdivide(A, B)

yields `A/B`, the right division of the affine transform `A` by the affine
transform `B`.

""" rightdivide

/(A::AffineTransform, B::AffineTransform) = rightdivide(A, B)

rightdivide(A::AffineTransform, B::AffineTransform) = begin
    Δ = det(B)
    iszero(Δ) && error("right operand is not invertible")
    Rxx = (A.xx*B.yy - A.xy*B.yx)/Δ
    Rxy = (A.xy*B.xx - A.xx*B.xy)/Δ
    Ryx = (A.yx*B.yy - A.yy*B.yx)/Δ
    Ryy = (A.yy*B.xx - A.yx*B.xy)/Δ
    return AffineTransform(Rxx, Rxy, A.x - (Rxx*B.x + Rxy*B.y),
                           Ryx, Ryy, A.y - (Ryx*B.y + Ryy*B.y))
end

"""

    leftdivide(A, B)

yields `A\\B`, the left division of the affine transform `A` by the affine
transform `B`.

""" leftdivide

\(A::AffineTransform, B::AffineTransform) = leftdivide(A, B)

leftdivide(A::AffineTransform, B::AffineTransform) = begin
    Δ = det(B)
    iszero(Δ) && error("left operand is not invertible")
    Rxx = A.yy/Δ
    Rxy = A.xy/Δ
    Ryx = A.yx/Δ
    Ryy = A.xx/Δ
    Rx = B.x - A.x
    Ry = B.y - A.y
    return AffineTransform(Rxx*B.xx - Rxy*B.yx,
                           Rxx*B.xy - Rxy*B.yy,
                           Rxx*Rx   - Rxy*Ry,
                           Ryy*B.yx - Ryx*B.xx,
                           Ryy*B.yy - Ryx*B.xy,
                           Ryy*Ry   - Ryx*Rx)
end

"""

`intercept(A)` returns the tuple `(x,y)` such that `A(x,y) = (0,0)`.

"""
function intercept(A::AffineTransform)
    Δ = det(A)
    iszero(Δ) && error("transformation is not invertible")
    return ((A.xy*A.y - A.yy*A.x)/Δ, (A.yx*A.x - A.xx*A.y)/Δ)
end

intercept(::Type{T}, A::AffineTransform) where {T<:Point} = T(intercept(A)...)


Base.show(io::IO, ::MIME"text/plain", A::AffineTransform) =
    print(io, typeof(A),
          "(",   A.xx, ",", A.xy, ",", A.x,
          ",  ", A.yx, ",", A.yy, ",", A.y, ")")

Base.show(io::IO, A::AffineTransform) = show(io, MIME"text/plain"(), A)

Base.print(io::IOBuffer, A::AffineTransform) = show(io, A)

end # module
