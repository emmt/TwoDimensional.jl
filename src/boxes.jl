"""
    BoundingBox(xmin,xmax,ymin,ymax)
    BoundingBox((xmin,ymin),(xmax,ymax))

yield an instance of a 2D rectangular bounding-box whose sides are aligned with
the coordinate axes and containing points of coordinates `(x,y)` such that
`xmin ≤ x ≤ xmax` and `ymin ≤ y ≤ ymax`. The box is *empty* if `xmin > xmax` or
`ymin > ymax`.

A bounding-box can be constructed from the first and last points (i.e. at the
lower-left and upper right opposite corners) of the box:

    BoundingBox(P0::Point, P1::Point)
    BoundingBox(I0::CartesianIndex{2}, I1::CartesianIndex{2})

Coordinates can be specified by keywords:

    BoundingBox(xmin=x0, ymin=y0, xmax=x1, ymax=y1)

There are no default values for keywords `xmin`, `xmax`, `ymin` and `ymax` so
all must be specified.

The coordinates of a `BoundingBox`, say `box`, can be retrieved as follows:

    box.xmin  or  box[1]  ->  xmin
    box.xmax  or  box[2]  ->  xmax
    box.ymin  or  box[3]  ->  ymin
    box.ymax  or  box[4]  ->  ymax

or:

    xmin, xmax, ymin, ymax = box

See also [`Point`](@ref), [`interior`](@ref), [`exterior`](@ref).

"""
BoundingBox{T}(xmin, xmax, ymin, ymax) where {T} = BoundingBox{T}((xmin, xmax, ymin, ymax))
BoundingBox(xmin, xmax, ymin, ymax) = BoundingBox(promote(xmin, xmax, ymin, ymax))
BoundingBox(xmin::T, xmax::T, ymin::T, ymax::T) where {T} =
    BoundingBox{T}(xmin, xmax, ymin, ymax)

# Box limits specified by keywords.
BoundingBox(; xmin, xmax, ymin, ymax) = BoundingBox(xmin, xmax, ymin, ymax)
BoundingBox{T}(; xmin, xmax, ymin, ymax) where {T} = BoundingBox{T}(xmin, xmax, ymin, ymax)

# Box specified by a ... box.
BoundingBox(box::BoundingBox) = box
BoundingBox{T}(box::BoundingBox{T}) where {T} = box
BoundingBox{T}(box::BoundingBox) where {T} = BoundingBox{T}(Tuple(box))

# Box limits specified as unit-ranges.
let type = :(AbstractUnitRange{<:Integer}),
    expr = (:(first(xrng)), :(last(xrng)), :(first(yrng)), :(last(yrng)))
    for args in ((:(xrng::$type), :(yrng::$type)),
                 (:((xrng,yrng)::Tuple{$type,$type}),))
        @eval begin
            BoundingBox($(args...)) = BoundingBox($(expr...))
            BoundingBox{T}($(args...)) where {T} = BoundingBox{T}($(expr...))
        end
    end
end

# Box limits specified by a pair of points, of Cartesian indices, or of 2-tuple.
let expr = (:(get_x(min)), :(get_x(max)), :(get_y(min)), :(get_y(max)))
    for type in (:Point, :(CartesianIndex{2}), :(Tuple{Any,Any}))
        for args in ((:(min::$type), :(max::$type)),
                     (:((min,max)::Tuple{$type,$type}),))
            @eval begin
                BoundingBox($(args...)) = BoundingBox($(expr...))
                BoundingBox{T}($(args...)) where {T} = BoundingBox{T}($(expr...))
            end
        end
    end
end

# Box limits specified by Cartesian indices.
BoundingBox(inds::CartesianIndices{2}) = BoundingBox(first(inds), last(inds))
BoundingBox{T}(inds::CartesianIndices{2}) where {T} = BoundingBox{T}(first(inds), last(inds))

# Properties.
Base.propertynames(::BoundingBox) = (:xmin, :xmax, :ymin, :ymax)
Base.getproperty(box::BoundingBox, key::Symbol) =
    key === :xmin ? box[1] :
    key === :xmax ? box[2] :
    key === :ymin ? box[3] :
    key === :ymax ? box[4] : throw(KeyError(key))

function Base.show(io::IO, box::BoundingBox{T}) where {T}
    print(io, "BoundingBox{")
    show(io, T)
    print(io, "}(xmin = "); show(io, box.xmin)
    print(io, ", xmax = "); show(io, box.xmax)
    print(io, ", ymin = "); show(io, box.ymin)
    print(io, ", ymax = "); show(io, box.ymax)
    print(io, ")")
end

# See
# https://stackoverflow.com/questions/9852159/calculate-bounding-box-of-arbitrary-pixel-based-drawing
# for the basic ideas under the following algorithm.
BoundingBox(A::AbstractMatrix{Bool}) = BoundingBox(identity, A)
function BoundingBox(f::Function, A::AbstractMatrix)
    I, J = axes(A)
    i0, i1 = get_axis_bounds(I)
    j0, j1 = get_axis_bounds(J)
    imin = jmin = typemax(Int)
    imax = jmax = typemin(Int)
    # Assuming column-major order, first start by scanning rows to narrow the
    # subsequent searches along columns.
    #
    # 1. Find bottom bound `jmin` by scanning rows from bottom to top.
    flag = false
    @inbounds for j in j0:j1, i in i0:i1
        if f(A[i,j])
            # This definitively set the value of `jmin` and gives limits for
            # the other bounds.
            imin = imax = i
            jmin = jmax = j
            flag = true
            break
        end
    end
    if flag
        # 2. Find top bound `jmax` by scanning rows from top to bottom.  No
        #    needs to go beyond `jmax+1`.
        @inbounds for j in j1:-1:jmax+1, i in i0:i1
            if f(A[i,j])
                jmax = j
                imin = min(imin, i)
                imax = max(imax, i)
                break
            end
        end
        # 3. Find leftmost bound `imin` by scanning columns from left to right.
        #    No needs to go beyond `imin-1`.
        @inbounds for i in i0:imin-1, j in jmin:jmax
            if f(A[i,j])
                imin = i
                imax = max(imax, i)
                break
            end
        end
        # 4. Find rightmost bound `imax` by scanning columns from right to
        #    left.  No needs to go beyond `imax+1`.
        @inbounds for i in i1:-1:imax+1, j in jmin:jmax
            if f(A[i,j])
                imax = i
                break
            end
        end
    end
    return BoundingBox(imin, imax, jmin, jmax)
end

"""
    TwoDimensional.get_axis_bounds(I) = (i0,i1)

yields the bounds `i0` and `i1` of index range `I` as a 2-tuple of `Int`'s and
such that `i0:i1` represents the same indices as `I` (although not in the same
order if `step(I) < 0`). If `step(I)` is not equal to ±1, an `ArgumentError`
exception is thrown.

"""
@inline function get_axis_bounds(I::AbstractRange{<:Integer})
    i0, i1, s = Int(first(I)), Int(last(I)), step(I)
    s == +one(s) && return (i0, i1)
    s == -one(s) && return (i1, i0)
    throw(ArgumentError("expecting a range with a step equal to ±1, got $s"))
end

# Empty bounding and unlimited boxes.
BoundingBox{T}(::Nothing) where {T<:Real} = typemin(BoundingBox{T})
Base.typemin(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemax(T), typemin(T), typemax(T), typemin(T))
Base.typemax(::Type{BoundingBox{T}}) where {T<:Real} =
    BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T))

# Conversion of bounding-boxes to/from Cartesian indices.
Base.CartesianIndices(box::BoundingBox{<:Integer}) =
    CartesianIndices((Int(box.xmin):Int(box.xmax), Int(box.ymin):Int(box.ymax)))

# Lower left and upper right corners of a bounding-box.
Base.first(box::BoundingBox) = Point(box.xmin, box.ymin)
Base.last(box::BoundingBox) = Point(box.xmax, box.ymax)

Base.isempty(box::BoundingBox) = (box.xmin > box.xmax)|(box.ymin > box.ymax)

Base.size(box::BoundingBox{<:Integer}) = map(length, axes(box))
Base.size(box::BoundingBox{<:Integer}, d::Integer) = length(axes(box, d))

Base.axes(box::BoundingBox{<:Integer}) = (UnitRange{Int}(box.xmin, box.xmax),
                                          UnitRange{Int}(box.ymin, box.ymax))
Base.axes(box::BoundingBox{<:Integer}, d::Integer) =
    d == 1 ? UnitRange{Int}(box.xmin, box.xmax) :
    d == 2 ? UnitRange{Int}(box.ymin, box.ymax) :
    d > 2 ? (1:1) : throw_bad_dimension_index()

@noinline throw_bad_dimension_index() =
    error("invalid dimension index")

# Use bounding-boxes to extract a sub-array or a view.
@propagate_inbounds function Base.getindex(A::AbstractMatrix,
                                           B::BoundingBox{<:Integer})
    A[B.xmin:B.xmax, B.ymin:B.ymax]
end

Base.view(A::AbstractMatrix, B::BoundingBox{<:Integer}) =
    view(A, B.xmin:B.xmax, B.ymin:B.ymax)

"""
    interior([T,] box)

yields the largest bounding-box with integer valued bounds and which is
contained by the bounding-box `box`. Optional argument `T` is to specify the
type of the result or of the coordinates of the result which is the same as
`box` by default.

See also: [`exterior`](@ref), [`round`](@ref).

"""
interior(box::BoundingBox{T}) where {T} = interior(T, box)
interior(::Type{BoundingBox{T}}, box::BoundingBox) where {T} = interior(T, box)
interior(::Type{T}, box::BoundingBox{T}) where {T<:Integer} = box
interior(::Type{T}, box::BoundingBox{T}) where {T<:Real} =
    BoundingBox(ceil(box.xmin), floor(box.xmax),
                ceil(box.ymin), floor(box.ymax))
interior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(box)
interior(::Type{T}, box::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(ceil(T, box.xmin), floor(T, box.xmax),
                ceil(T, box.ymin), floor(T, box.ymax))
interior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(box)
interior(::Type{T}, box::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(interior(U, box))

"""
    exterior([T,] box)

yields the smallest bounding-box with integer valued bounds and which contains
the bounding-box `box`.  Optional argument `T` is to specify the type of the
result or of the coordinates of the result which is the same as `box` by default.

See also: [`interior`](@ref), [`round`](@ref).

"""
exterior(box::BoundingBox{T}) where {T} = exterior(T, box)
exterior(::Type{BoundingBox{T}}, box::BoundingBox) where {T} = exterior(T, box)
exterior(::Type{T}, box::BoundingBox{T}) where {T<:Integer} = box
exterior(::Type{T}, box::BoundingBox{T}) where {T<:Real} =
    BoundingBox(floor(box.xmin), ceil(box.xmax),
                floor(box.ymin), ceil(box.ymax))
exterior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Integer} =
    BoundingBox{T}(box)
exterior(::Type{T}, box::BoundingBox{<:Real}) where {T<:Integer} =
    BoundingBox(floor(T, box.xmin), ceil(T, box.xmax),
                floor(T, box.ymin), ceil(T, box.ymax))
exterior(::Type{T}, box::BoundingBox{<:Integer}) where {T<:Real} =
    BoundingBox{T}(box)
exterior(::Type{T}, box::BoundingBox{U}) where {T<:Real,U<:Real} =
    BoundingBox{T}(exterior(U, box))

"""
    center(box::BoundingBox) -> c::Point

yields the central point of the bounding-box `box`.

"""
center(box::BoundingBox) =
    isempty(box) ? throw(ArgumentError("cannot get center of empty box")) :
    Point((box.xmin + box.xmax)/2, (box.ymin + box.ymax)/2)

"""
    area(box)

yields the area of the bounding-box `box`.

"""
area(box::BoundingBox{T}) where {T<:Real} =
    max(box.xmax - box.xmin, zero(T))*max(box.ymax - box.ymin, zero(T))
