"""
    msk = MaskElement{T}(shape::ShapeElement; opaque=true/false)

builds an elementary mask whose shape is given by the elementary geometric
object `shape` and with keyword `opaque` set to `true` for an opaque mask and
to `false` for a transparent mask. `T` is the coordinate type which may be
omitted.

To change the opacity (by default the same opacity is kept) and/or the
coordinate type:

    other_msk = MaskElement{T′}(msk; opaque=...)

All arithmetic operations on elementary geometric objects also apply for a
mask: the operation is applied to the shape of the mask leaving the opacity
unchanged. To toggle the opacity:

    other_msk = inv(msk)

To check the opacity, call
[`TwoDimensional.is_opaque(msk)`](@ref TwoDimensional.is_opaque) or
[`TwoDimensional.is_transparent(msk)`](@ref TwoDimensional.is_transparent).

There are several helper functions to make the code clearer when building
a mask (*obscurations* are opaque mask while *aperture* are transparent masks):
[`TwoDimensional.rectangular_aperture`](@ref TwoDimensional.rectangular_aperture),
[`TwoDimensional.rectangular_obscuration`](@ref TwoDimensional.rectangular_obscuration),
[`TwoDimensional.circular_aperture`](@ref TwoDimensional.circular_aperture),
[`TwoDimensional.circular_obscuration`](@ref TwoDimensional.circular_obscuration),
[`TwoDimensional.polygonal_aperture`](@ref TwoDimensional.polygonal_aperture), or
[`TwoDimensional.polygonal_obscuration`](@ref TwoDimensional.polygonal_obscuration).

A mask is a combination of one or several elementary mask objects (apertures
and/or obscurations). To forge a mask, call
[`TwoDimensional.forge_mask`](@ref TwoDimensional.forge_mask) or
[`TwoDimensional.forge_mask!`](@ref TwoDimensional.forge_mask!).
A forged mask can serve multiple times. To apply a mask once, call
[`TwoDimensional.apply_mask`](@ref TwoDimensional.apply_mask) or
[`TwoDimensional.apply_mask!`](@ref TwoDimensional.apply_mask!).

"""
MaskElement(msk::MaskElement; opaque::Bool = is_opaque(msk)) =
    MaskElement(shape(obj); opaque = opaque)
MaskElement{T}(msk::MaskElement; opaque::Bool = is_opaque(msk)) where {T} =
    MaskElement(convert_coord_type(T, shape(obj)); opaque = opaque)
MaskElement{T}(obj::ShapeElement; opaque::Bool) where {T} =
    MaskElement(convert_coord_type(T, obj), opaque)

# Extend convert to call constructor.
Base.convert(::Type{T}, obj::T) where {T<:GeometricObject} = obj
Base.convert(::Type{T}, obj) where {T<:GeometricObject} = T(obj)

# Impose that collected mask elements have the same coordinate type.
Base.collect(obj::MaskElement, objs::MaskElement...) = collect_mask_elements(obj, objs...)
collect_mask_elements() = error("expecting at least one mask element")
collect_mask_elements(objs::MaskElement{T}...) where {T} = MaskElement{T}[objs...]
function collect_mask_elements(objs::MaskElement...)
    T = promote_type(map(coord_type, objs))
    return MaskElement{T}[obj, objs...]
end

Base.first(obj::RectangularMask) = first(shape(obj))
Base.last(obj::RectangularMask) = last(shape(obj))

"""
    TwoDimensional.is_opaque(msk)

yields whether mask element `msk` is opaque.

"""
is_opaque(obj::MaskElement) = obj.opaque

"""
    TwoDimensional.is_transparent(msk)

yields whether mask element `msk` is transparent.

"""
is_transparent(obj::MaskElement) = !is_opaque(obj)

"""
    TwoDimensional.rectangular_aperture(args...; kwds...)

yields an elementary mask object representing a rectangular aperture defined by
given arguments `args...` and keywords `kwds...` and whose edges are aligned
with the Cartesian axes. See [`TwoDimensional.Rectangle`](@ref) constructor for
possible arguments and keywords. A rectangular aperture is a transparent
rectangular mask.

""" rectangular_aperture

"""
    TwoDimensional.rectangular_obscuration(args...; kwds...)

yields an elementary mask object representing a rectangular obscuration defined
by given arguments `args...` and keywords `kwds...` and whose edges are aligned
with the Cartesian axes. See [`TwoDimensional.Rectangle`](@ref) constructor for
possible arguments and keywords. A rectangular obscuration is an opaque
rectangular mask.

""" rectangular_obscuration

"""
    TwoDimensional.circular_aperture(args...; kwds...)

yields an elementary mask object representing a circular aperture defined by
given arguments `args...` and keywords `kwds...`. See
[`TwoDimensional.Circle`](@ref) constructor for possible arguments and keywords.
A circular aperture is a transparent circular mask.

""" circular_aperture

"""
    TwoDimensional.circular_obscuration(args...; kwds...)

yields an elementary mask object representing a circular obscuration defined by
given arguments `args...` and keywords `kwds...`. See
[`TwoDimensional.Circle`](@ref) constructor for possible arguments and keywords.
A circular obscuration is an opaque circular mask.

""" circular_obscuration

"""
    TwoDimensional.polygonal_aperture(args...; kwds...)

yields an elementary mask object representing a polygonal aperture defined by
given arguments `args...` and keywords `kwds...`. See
[`TwoDimensional.Polygon`](@ref) constructor for possible arguments and
keywords. A polygonal aperture is a transparent polygonal mask.

""" polygonal_aperture

"""
    TwoDimensional.polygonal_obscuration(args...; kwds...)

yields an elementary mask object representing a polygonal obscuration defined
by given arguments `args...` and keywords `kwds...`. See
[`TwoDimensional.Polygon`](@ref) constructor for possible arguments and
keywords. A polygonal obscuration is an opaque polygonal mask.

""" polygonal_obscuration

for (type, adj) in ((:Rectangle, :rectangular),
                    (:Circle,    :circular),
                    (:Polygon,   :polygonal),)
    @eval begin
        $(Symbol(adj,"_aperture"))(args...; kwds...) =
            MaskElement($type(args...; kwds...); opaque = false)
        $(Symbol(adj,"_obscuration"))(args...; kwds...) =
            MaskElement($type(args...; kwds...); opaque = true)
    end
end

"""
    TwoDimensional.apply_mask(Aᵢₙ, args...; kwds...) -> Aₒᵤₜ

multiplies the values of the input 2-dimensional array `Aᵢₙ` by a mask defined
by arguments `args...` and keywords `kwds...` and returns the resulting output
array `Aₒᵤₜ`. The input array `Aᵢₙ` is left unmodified, method
[`TwoDimensional.apply_mask!`](@ref) may be used for in-place operation. See
[`TwoDimensional.forge_mask`](@ref) for how to define a mask.

"""
apply_mask(A::AbstractMatrix, args...; kwds...) = apply_mask!(copy(A), args...; kwds...)

"""
    TwoDimensional.apply_mask!(A, args...) -> A

multiplies in-place the 2-dimensional array `A` by a mask defined by arguments
`args...` and keywords `kwds...` and returns `A`. See
[`TwoDimensional.apply_mask`](@ref) for an out-of-place version and for
details.

"""
apply_mask!(A::AbstractMatrix, args...; kwds...) = multiply!(A, forge_mask(A, args...; kwds...))

# Element-wise multiplication.
function multiply!(A::AbstractArray{<:Any,N}, B::AbstractArray{<:Any,N}) where {N}
    axes(A) == axes(B) || throw(DimensionMismatch("arguments have different axes"))
    @. A *= B
    return A
end

# Yields a range with `n` points centered around zero and step `s/n`. This
# method is used to index a sub-cell.
function subrange(n::Int, s::Number)
    r = (1:n) .* (s/n)
    return r .- (first(r) + last(r))/2
end

"""
    TwoDimensional.forge_mask(A::AbstractMatrix, objs...; kwds...) -> msk

yields a 2-dimensional transmission mask for the 2-dimensional array `A` and
combining elementary mask objects `objs...`.

"""
forge_mask(A::AbstractMatrix, args...; kwds...) =
    forge_mask!(Matrix{floating_point_type(A)}(undef, size(A)), A, args...; kwds...)

function forge_mask!(dst::AbstractMatrix, A::AbstractMatrix, args...; kwds...)
    X = grid_step(A)*RolledCoordinates(grid_size(A))
    Y = X
    return forge_mask!(dst, X, Y, args...; kwds...)
end

"""
    TwoDimensional.forge_mask(X, Y, objs...; kwds...) -> msk

yields a 2-dimensional transmission mask with coordinates given by `X` and `Y`
along the 1st and 2nd dimensions and combining aperture/obscuration objects
`objs...`. The following *painting* algorithm is used:

- The mask is initially filled with the transparent or opaque value depending
  on whether the first component is opaque or transparent.

- Then, for each component in turn, the cells of the mask that are inside the
  component are painted with the opaque or transparent value depending on
  whether the component is opaque or transparent.

- The parts of the mask overlapping the boundaries of the topmost components
  are set to an intermediate value between the opaque and transparent ones and
  (approximately) proportionally to the transparent fraction of the cell area.

Note that the order of the components of the mask is relevant: an aperture
component drills holes in the previously opaque parts while an obscuration
hides previously transparent parts.

Keyword `antialiasing` can be set to specify the number of sub-cells (per side)
to determine the transmission of grid cells partially overlapping the boundary
delimiting the mask components. By default, `antialiasing =
$default_antialiasing`. If `antialiasing ≤ 1`, a 50% transmission is assumed
for partially overlapping cells (sharp edges); otherwise, overlapping cells are
subdivided in `antialiasing × antialiasing` sub-cells to estimate their partial
transmission.

Keywords `opaque` and `transparent` can be used to specify the values of the
the respectively opaque and transparent parts of the mask. Values of partially
opaque/transparent parts will be interpolated between these.

Example to forge a mask representing the primary mirror of a telescope with its
spider arms:

    using TwoDimensional
    using Unitful: μm, mm, cm, m
    center = Point(0mm,0mm) # central position
    outer_radius = 1.8m
    inner_radius = outer_radius/3
    spider_thickness = 2.3cm # thickness of spider arms
    grid_step = 2.0mm # grid sampling step
    spider_length = 2*(outer_radius + grid_step) # length of spider arms
    margin = 5 # margin in pixels
    xmin = floor(Int, (center.x - outer_radius)/grid_step) - margin
    ymin = floor(Int, (center.y - outer_radius)/grid_step) - margin
    xmax = ceil(Int, (center.x + outer_radius)/grid_step) + margin
    ymax = ceil(Int, (center.y + outer_radius)/grid_step) + margin
    X = (xmin:xmax)*grid_step # coordinates along 1st dimension
    Y = (ymin:ymax)*grid_step # coordinates along 2nd dimension
    voff = Point(spider_thickness, spider_length)
    hoff = Point(spider_length, spider_thickness)
    mask = forge_mask(
        X, Y,
        circular_aperture(center, outer_radius), # aperture
        circular_obscuration(center, inner_radius), # central obscuration
        rectangular_obscuration(center - voff/2, center + voff/2),
        rectangular_obscuration(center - hoff/2, center + hoff/2))

"""
function forge_mask(X::AbstractVector,
                    Y::AbstractVector,
                    args...; kwds...)
    T = floating_point_type(promote_type(eltype(X), eltype(Y)))
    return forge_mask!(Matrix{T}(undef, length(X), length(Y)), X, Y, args...; kwds...)
end


"""
    TwoDimensional.forge_mask!(dst, X, Y, objs...; kwds...) -> dst

In-place version of [`TwoDimensional.forge_mask`](@ref), it overwrites the
destination array `dst` with the mask and returns it.

"""
function forge_mask!(dst::AbstractMatrix,
                     X::AbstractVector,
                     Y::AbstractVector,
                     args::Tuple{Vararg{MaskElement}};
                     kwds...)
    return forge_mask!(dst, X, Y, args...; kwds...)
end

function forge_mask!(dst::AbstractMatrix,
                     X::AbstractVector,
                     Y::AbstractVector,
                     args::MaskElement...;
                     antialiasing::Integer = default_antialiasing,
                     opaque = zero(eltype(dst)),
                     transparent = oneunit(eltype(dst)))
    # Check that coordinate units are compatible and determine a suitable unit
    # for all coordinates.
    T = float(promote_type(eltype(X), eltype(Y), map(coord_type, args)...))

    # Check arguments.
    I, J = axes(dst)
    axes(X) == (I,) || throw(DimensionMismatch("coordinates `X` have incompatible indices"))
    axes(Y) == (J,) || throw(DimensionMismatch("coordinates `Y` have incompatible indices"))
    antialiasing = as(Int , antialiasing)
    opaque = as(eltype(dst), opaque)
    transparent = as(eltype(dst), transparent)

    # Call the real method.
    return unsafe_forge_mask!(dst, convert_eltype(T, X), convert_eltype(T, Y),
                              map(Fix1(convert_coord_type, T), args)...;
                              antialiasing = antialiasing, opaque = opaque,
                              transparent = transparent)
end

function unsafe_forge_mask!(dst::AbstractMatrix{T},
                            X::AbstractVector{C},
                            Y::AbstractVector{C},
                            objs::MaskElement{C}...;
                            antialiasing::Int,
                            opaque::T,
                            transparent::T) where {T,C}
    partial = interpolate(opaque, transparent, 1//2)
    δx = grid_step(X)
    δy = grid_step(Y)

    # Determine whether a cell of the mask is fully or partially opaque of
    # transparent.
    initial = true
    for obj in objs
        unsafe_forge_mask!(dst, X, δx, Y, δy, obj, opaque, partial, transparent, initial)
        initial = false
    end
    if initial
        # No objects, assume mask is transparent.
        fill!(mask, transparent)
    end

    if antialiasing > 1
        # Refine mask at partially opaque cells.
        I, J = axes(dst)
        DX = subrange(antialiasing, δx)
        DY = subrange(antialiasing, δy)
        state = Array{Bool}(undef, antialiasing, antialiasing)
        @inbounds for j in J
            y = Y[j]
            for i in I
                if dst[i,j] == partial
                    # Cell is partially opaque/transparent.
                    x = X[i]
                    count = -1 # to trigger resetting of state array
                    for obj in objs
                        count = unsafe_forge_mask!(state, x, DX, y, DY, obj, count)
                    end
                    #@assert max(count, 0) == Base.count(state)
                    fraction = max(count, 0)//antialiasing^2
                    dst[i,j] = interpolate(opaque, transparent, fraction)
                end
            end
        end
    end

    return dst
end

function unsafe_forge_mask!(dst::AbstractMatrix{V},
                            X::AbstractVector{T}, δx::T,
                            Y::AbstractVector{T}, δy::T,
                            mask::MaskElement{T},
                            opaque::V,
                            partial::V,
                            transparent::V,
                            initial::Bool) where {T,V}
    if initial
        fill!(dst, is_opaque(mask) ? transparent : opaque)
    end
    return unsafe_forge_mask!(dst, X, δx, Y, δy, shape(mask),
                              is_opaque(mask) ? opaque : transparent,
                              partial)
end

function unsafe_forge_mask!(dst::AbstractMatrix{V},
                            X::AbstractVector{T}, δx::T,
                            Y::AbstractVector{T}, δy::T,
                            obj::ShapeElement{T},
                            value::V,
                            partial::V) where {T,V}
    box = grow(BoundingBox(obj), δx, δy)
    I, J = axes(dst)
    @inbounds for j in J
        y = Y[j]
        if (y < box.ymin)|(y > box.ymax)
            # Cell is outside the bounding-box of the mask.
            continue
        end
        for i in I
            x = X[i]
            if (x < box.xmin)|(x > box.xmax)
                # Cell is outside the bounding-box of the mask.
                continue
            end
            # FIXME We know that coordinates are in order...
            cell = Rectangle((x - δx/2, y - δy/2), (x + δx/2, y + δy/2))
            overlap = Overlap(cell, obj)
            if overlap == INSIDE
                dst[i,j] = value
            elseif overlap == PARTIAL
                dst[i,j] = partial
            end
        end
    end
    return dst
end

function unsafe_forge_mask!(state::AbstractMatrix{Bool},
                            x0::T, dx::AbstractVector{T},
                            y0::T, dy::AbstractVector{T},
                            obj::MaskElement{T},
                            count::Int) where {T}
    # FIXME: Skip all this if pixel is outside boundaries.
    I, J = axes(state)
    if count < 0
        # This is the first mask applied to this cell.
        count = 0
        @inbounds for j in J
            y = y0 + dy[j]
            for i in I
                x = x0 + dx[i]
                inside = Overlap(Point(x,y), obj) != OUTSIDE
                transp = xor(inside, is_opaque(obj))
                state[i,j] = transp
                count += transp
            end
        end
    elseif is_opaque(obj)
        @inbounds for j in J
            y = y0 + dy[j]
            for i in I
                x = x0 + dx[i]
                if Overlap(Point(x,y), obj) != OUTSIDE
                    # Point is considered as being inside boundaries of opaque
                    # mask.
                    if state[i,j]
                        # Sub-cell previously counted as transparent.
                        state[i,j] = false
                        count -= 1
                    end
                end
            end
        end
    else # element is transparent
        @inbounds for j in J
            y = y0 + dy[j]
            for i in I
                x = x0 + dx[i]
                if Overlap(Point(x,y), obj) != OUTSIDE
                    # Point is considered as being inside boundaries of transparent mask.
                    if !state[i,j]
                        # Sub-cell previously considered as opaque.
                        state[i,j] = true
                        count += 1
                    end
                end
            end
        end
    end
    return count
end

"""
    TwoDimensional.grid_step(x::AbstractVector) -> stp

yields the step along a vector of coordinates, throwing an error if the
increment between successive values of `x` is not positive or not uniform.

"""
function grid_step(rng::Union{AbstractRange#= FIXME ,AbstractCoordinates =#})
    stp = step(rng)
    stp > zero(stp) || throw(ArgumentError("grid step must be positive"))
    return stp
end

function grid_step(x::AbstractVector)
    len = length(x)
    len > 1 || throw(ArgumentError("zero-length vector of coordinates has no defined step"))
    stp = (last(x) - first(x))/(len - 1)
    stp > zero(stp) || throw(ArgumentError("grid step must be positive"))
    tol = sqrt(eps(floating_point_type(stp)))
    stpmin = stp*(one(tol) - tol)
    stpmax = stp*(one(tol) + tol)
    flag = true
    @inbounds for i ∈ firstindex(x):lastindex(x)-1
        flag &= (stpmin ≤ x[i+1] - x[i] ≤ stpmax)
    end
    flag || throw(ArgumentError("vector of coordinates has non-uniform steps"))
    return stp
end

"""
    TwoDimensional.interpolate(a, b, f) -> x

yields linearly interpolated value between `a` and `b` by a fraction `f`. If
`f` is a dimensionless factor, then the result is:

    x = a + f*(b - a)

otherwise, e.g. if `f` has units, the result is:

    x = (oneunit(f) - f)*a + f*b

In any case, `a` and `b` are promoted to the same type.

Hence, if `f` is a real, then `f = 0` yields `a` while `f = 1` yields `b`.

"""
interpolate(a, b, f::Number) = interpolate(promote(a, b)..., f)
interpolate(a::T, b::T, f::Real) where {T} = a + f*(b - a)
interpolate(a::T, b::T, f::Number) where {T} = (oneunit(f) - f)*a + f*b

"""
    TwoDimensional.Overlap(pxl, obj)

yields the overlapping of pixel `pxl` with shape object `obj`. `pxl` may be a
point, a bounding-box, or a rectangle. Returned value is:

- `INSIDE` if `pxl` is fully inside the boundaries of `obj`;

- `OUTSIDE` if `pxl` is fully outside the boundaries of `obj`;

- `PARTIAL` if `pxl` straddles some boundaries of `obj`.

"""
Overlap(pxl::Union{Point,Rectangle}, obj::MaskElement) = Overlap(pxl, shape(obj))

Overlap(pxl::BoundingBox, obj::GeometricObject) =
    isempty(pxl) ? OUTSIDE : Overlap(Rectangle(pxl), obj)

Overlap(pxl::Union{Point,Rectangle}, obj::Union{Rectangle,Circle,Polygon}) =
    # FIXME: Converting to the same coordinate type has some cost...
    Overlap(promote_coord_type(pxl, obj)...)

function Overlap(point::Point{T}, rect::Rectangle{T}) where {T}
    x, y = point
    (xmin, ymin), (xmax, ymax) = rect
    if (xmin < x < xmax) & (ymin < y < ymax)
        return INSIDE
    elseif (x < xmin) | (x > xmax) | (y < ymin) | (y > ymax)
        return OUTSIDE
    else
        return PARTIAL
    end
end

function Overlap(cell::Rectangle{T}, rect::Rectangle{T}) where {T}
    (x0, y0), (x1, y1) = cell
    (xmin, ymin), (xmax, ymax) = rect
    if (x0 ≥ xmin) & (x1 ≤ xmax) & (y0 ≥ ymin) & (y1 ≤ ymax)
        return INSIDE
    elseif (x0 > xmax) | (x1 < xmin) | (y0 > ymax) | (y1 < ymin)
        return OUTSIDE
    else
        return PARTIAL
    end
end

function Overlap(point::Point{T}, circle::Circle{T}) where {T}
    d² = abs2(point - circle.center)
    r² = abs2(circle.radius)
    if d² < r²
        return INSIDE
    elseif d² == r²
        return PARTIAL
    else
        return OUTSIDE
    end
end

function Overlap(cell::Rectangle{T}, circle::Circle{T}) where {T}
    # Get coordinates of cell corners relative to circle center and quickly
    # check whether cell is certainly outside circle.
    xc, yc = circle.center
    x0 = cell.x0 - xc
    x1 = cell.x1 - xc
    y0 = cell.y0 - yc
    y1 = cell.y1 - yc
    r = circle.radius
    if (x0 > r) | (x1 < -r) | (y0 > r) | (y1 < -r)
        return OUTSIDE
    end
    r² = r^2
    # Most distant point in cell from center of circle is one of the corners.
    # If this point is inside the circle, then the cell is fully inside the
    # circle.
    x²max = max(x0^2, x1^2)
    y²max = max(y0^2, y1^2)
    if x²max + y²max ≤ r²
        return INSIDE
    end
    # Otherwise, if closest point of one of the cell edge is inside the circle,
    # overlap is partial.
    if clamp(zero(T), x0, x1)^2 + y²max ≤ r²
        return PARTIAL
    end
    if x²max + clamp(zero(T), y0, y1)^2 ≤ r²
        return PARTIAL
    end
    # Otherwise, there is no overlap.
    return OUTSIDE
end

# FIXME: A point can only inside or outside.
Overlap(point::Point, polygon::Polygon) = point ∈ polygon ? INSIDE : OUTSIDE

# FIXME: Check that the following algorithm is correct for a cell having zero
#        width or height.
function Overlap(cell::Rectangle{T}, polygon::Polygon{T}) where {T}
    (x0, y0), (x1, y1) = cell # retrieve coordinates of cell vertices
    V = vertices(polygon) # retrieve the list of vertices of the polygon

    # If any polygon edges (strictly) crosses one of the 4 edges of the cell, then
    # there is some partial overlapping.
    A = last(V) # initialize A, the 1st point of edges in cyclic list of points
    @inbounds for B in V # loop over B, the 2nd point of edges
        xmin, xmax = minmax(A.x, B.x)
        ymin, ymax = minmax(A.y, B.y)
        if xmin < x0 < xmax
            # Left edge of cell may cross the segment AB.
            y = interpolate(A.y, B.y, (x0 - A.x)/(B.x - A.x))
            if y0 < y < y1
                return PARTIAL
            end
        end
        if xmin < x1 < xmax
            # Right edge of cell may cross the segment AB.
            y = interpolate(A.y, B.y, (x1 - A.x)/(B.x - A.x))
            if y0 < y < y1
                return PARTIAL
            end
        end
        if ymin < y0 < ymax
            # Bottom edge of cell may cross the segment AB.
            x = interpolate(A.x, B.x, (y0 - A.y)/(B.y - A.y))
            if x0 < x < x1
                return PARTIAL
            end
        end
        if ymin < y1 < ymax
            # Top edge of cell may cross the segment AB.
            x = interpolate(A.x, B.x, (y1 - A.y)/(B.y - A.y))
            if x0 < x < x1
                return PARTIAL
            end
        end
        A = B # update 1st point of next edge
    end

    # Otherwise, the cell is either fully inside or fully outside the polygon.
    # This can be tested for any point of the cell, its center for example.
    return Overlap(Point((x0 + x1)/2, (y0 + y1)/2), polygon)
end
