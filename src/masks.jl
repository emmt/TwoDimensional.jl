"""
    TwoDimensional.Mask(elems...)
    TwoDimensional.Mask{T}(elems...)

build a composite mask consisting in the ordered list of mask elements `elems...`.
Optional type parameter `T` is the coordinate type of the masks elements. All mask
elements are converted as needed to have this coordinate type. If `T` is not specified, it
is inferred by promoting the coordinate types of the mask elements.

A mask can be moved, scaled, rotated, etc., the coordinate type of its elements may be
converted to another type.

!!! warning
    Masks with a large number of elements should preferably be created with a vector, not
    a tuple, of mask elements.

"""
Mask(elems::MaskElement...) = Mask(elems)
Mask{T}(elems::MaskElement...) where {T} = Mask{T}(elems)
Mask(elems::List{<:MaskElement}) = Mask{coord_type(elems)}(elems)
Mask{T}(elems::List{<:MaskElement}) where {T} =
    Mask{T}(map(Fix1(convert_coord_type, T), elems))
Mask{T}(elems::List{<:MaskElement{T}}) where {T} = Mask{T,eltype(elems)}(elems)

# Extend abstract array API for masks.
Base.values(msk::Mask) = elements(msk)
Base.length(msk::Mask) = length(elements(msk))
Base.size(msk::Mask) = (length(msk),)
Base.axes(msk::Mask) = (Base.OneTo(length(msk)),)
Base.IndexStyle(::Type{<:Mask}) = IndexLinear()
@inline function Base.getindex(A::Mask, i::Int)
    @boundscheck checkbounds(A, i)
    return @inbounds getindex(elements(A), i)
end
@inline function Base.setindex!(A::Mask, x, i::Int)
    @boundscheck checkbounds(A, i)
    @inbounds setindex!(elements(A), x, i)
    return A
end

"""
    msk = MaskElement{T}(shape::ShapeElement; opaque)

builds an elementary mask whose shape is given by the elementary geometric object `shape`
and with keyword `opaque` set to `true` for an opaque mask (an *obscuration*) and to
`false` for a transparent mask (an *aperture*). `T` is the coordinate type which may be
omitted.

To change the opacity (by default the same opacity is kept) and/or the coordinate type:

    other_msk = MaskElement{T′}(msk; opaque=...)

All arithmetic operations on elementary geometric objects also apply for a mask: the
operation is applied to the shape of the mask leaving the opacity unchanged. To toggle the
opacity:

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
MaskElement(obj::MaskElement; opaque::Bool = is_opaque(obj)) =
    MaskElement(shape(obj); opaque = opaque)
MaskElement{T}(obj::MaskElement; opaque::Bool = is_opaque(obj)) where {T} =
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
collect_mask_elements(objs::MaskElement...) = MaskElement{coord_type(objs)}[objs...]

Base.first(obj::RectangularMask) = first(shape(obj))
Base.last(obj::RectangularMask) = last(shape(obj))

"""
    TwoDimensional.is_opaque(msk)

yields whether mask element `msk` is opaque.

"""
is_opaque(msk::MaskElement) = msk.opaque

"""
    TwoDimensional.is_transparent(msk)

yields whether mask element `msk` is transparent.

"""
is_transparent(msk::MaskElement) = !is_opaque(msk)

"""
    TwoDimensional.aperture(obj)

yields a transparent mask element of same shape as `obj`.

"""
aperture(obj::ShapeElement) = MaskElement(obj; opaque=false)
aperture(obj::MaskElement) = aperture(shape(obj))

"""
    TwoDimensional.obscuration(obj)

yields an opaque mask element of same shape as `obj`.

"""
obscuration(obj::ShapeElement) = MaskElement(obj; opaque=true)
obscuration(obj::MaskElement) = obscuration(shape(obj))

"""
    TwoDimensional.rectangular_aperture(args...; kwds...)

yields an elementary mask object representing a rectangular aperture defined by given
arguments `args...` and keywords `kwds...` and whose edges are aligned with the Cartesian
axes. See [`TwoDimensional.Rectangle`](@ref) constructor for possible arguments and
keywords. A rectangular aperture is a transparent rectangular mask.

""" rectangular_aperture

"""
    TwoDimensional.rectangular_obscuration(args...; kwds...)

yields an elementary mask object representing a rectangular obscuration defined by given
arguments `args...` and keywords `kwds...` and whose edges are aligned with the Cartesian
axes. See [`TwoDimensional.Rectangle`](@ref) constructor for possible arguments and
keywords. A rectangular obscuration is an opaque rectangular mask.

""" rectangular_obscuration

"""
    TwoDimensional.circular_aperture(args...; kwds...)

yields an elementary mask object representing a circular aperture defined by given
arguments `args...` and keywords `kwds...`. See [`TwoDimensional.Circle`](@ref)
constructor for possible arguments and keywords. A circular aperture is a transparent
circular mask.

""" circular_aperture

"""
    TwoDimensional.circular_obscuration(args...; kwds...)

yields an elementary mask object representing a circular obscuration defined by given
arguments `args...` and keywords `kwds...`. See [`TwoDimensional.Circle`](@ref)
constructor for possible arguments and keywords. A circular obscuration is an opaque
circular mask.

""" circular_obscuration

"""
    TwoDimensional.polygonal_aperture(args...; kwds...)

yields an elementary mask object representing a polygonal aperture defined by given
arguments `args...` and keywords `kwds...`. See [`TwoDimensional.Polygon`](@ref)
constructor for possible arguments and keywords. A polygonal aperture is a transparent
polygonal mask.

""" polygonal_aperture

"""
    TwoDimensional.polygonal_obscuration(args...; kwds...)

yields an elementary mask object representing a polygonal obscuration defined by given
arguments `args...` and keywords `kwds...`. See [`TwoDimensional.Polygon`](@ref)
constructor for possible arguments and keywords. A polygonal obscuration is an opaque
polygonal mask.

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

multiplies the values of the input 2-dimensional array `Aᵢₙ` by a mask defined by
arguments `args...` and keywords `kwds...` and returns the resulting output array `Aₒᵤₜ`.
The input array `Aᵢₙ` is left unmodified, method [`TwoDimensional.apply_mask!`](@ref) may
be used for in-place operation. See [`TwoDimensional.forge_mask`](@ref) for how to define
a mask.

"""
apply_mask(A::AbstractMatrix, args...; kwds...) = apply_mask!(copy(A), args...; kwds...)

"""
    TwoDimensional.apply_mask!(A, args...; kwds...) -> A

multiplies in-place the 2-dimensional array `A` by a mask defined by arguments `args...`
and keywords `kwds...` and returns `A`. See [`TwoDimensional.apply_mask`](@ref) for an
out-of-place version and for details.

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
    TwoDimensional.forge_mask(A::AbstractMatrix, msk; kwds...)
    TwoDimensional.forge_mask(A::AbstractMatrix, objs...; kwds...)

yield a 2-dimensional array with entries set to the transmission by the mask `msk` for an
array-like `A`. The mask may also be specified by the list `objs...` of elementary mask
objects. The coordinates of the mask are assumed to be given in fractional Cartesian
indices for `A`.

"""
forge_mask(A::AbstractMatrix, objs::MaskElement...; kwds...) = forge_mask(A, Mask(objs); kwds...)
forge_mask(A::AbstractMatrix, msk::Mask; kwds...) =
    forge_mask!(similar(A, floating_point_type(A)), msk; kwds...)

"""
    TwoDimensional.forge_mask([T,] X, Y, msk; kwds...) -> arr
    TwoDimensional.forge_mask([T,] X, Y, elems...; kwds...) -> arr

yield a 2-dimensional array filled with transmission values computed for the mask `msk` at
coordinates given by `X` and `Y` along the 1st and 2nd dimensions. The mask may also be
specified by combining elementary mask objects `elems...`. Optional argument `T` is to specify
the element type of the result.

The following *painting* algorithm is used:

- The mask array is initially filled with the transparent or opaque value depending on
  whether the first component is opaque or transparent.

- Then, for each component in turn, the cells of the mask array that are inside the
  component are painted with the opaque or transparent value depending on whether the
  component is opaque or transparent.

- The cells of the mask array overlapping the boundaries of the topmost components are set
  to an intermediate value between the opaque and transparent ones and (approximately)
  proportionally to the transparent fraction of the cell area.

Note that the order of the components of the mask is relevant: an aperture component
drills holes in the previously opaque parts while an obscuration hides previously
transparent parts.

Keyword `antialiasing` can be set to specify the number of sub-cells (per side) to
determine the transmission of grid cells partially overlapping the boundary delimiting the
mask components. By default, `antialiasing = $default_antialiasing`. If `antialiasing ≤
1`, a 50% transmission is assumed for partially overlapping cells (sharp edges);
otherwise, overlapping cells are subdivided in `antialiasing × antialiasing` sub-cells to
estimate their partial transmission.

Keywords `opaque` and `transparent` can be used to specify the values of the the
respectively opaque and transparent parts of the mask. Values of partially
opaque/transparent parts will be interpolated between these.

Keyword `multithreading` specifies whether to use multiple threads for the computations.

Example to forge a mask representing the primary mirror of a telescope with its spider
arms:

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
forge_mask(X::AbstractVector, Y::AbstractVector, elems::MaskElement...; kwds...) =
    forge_mask(X, Y, Mask(elems); kwds...)
forge_mask(::Type{T}, X::AbstractVector, Y::AbstractVector, elems::MaskElement...; kwds...) where {T} =
    forge_mask(T, X, Y, Mask(elems); kwds...)

forge_mask(X::AbstractVector, Y::AbstractVector, elems::List{<:MaskElement}; kwds...) =
    forge_mask(X, Y, Mask(elems); kwds...)
forge_mask(::Type{T}, X::AbstractVector, Y::AbstractVector, elems::List{<:MaskElement}; kwds...) where {T} =
    forge_mask(T, X, Y, Mask(elems); kwds...)

forge_mask(X::AbstractVector, Y::AbstractVector, msk::Mask; kwds...) =
    forge_mask(floating_point_type(promote_type(eltype(X), eltype(Y), coord_type(msk))),
               X, Y, msk; kwds...)
forge_mask(X::AbstractVector{T}, Y::AbstractVector{T}, msk::Mask{T}; kwds...) where {T} =
    forge_mask(floating_point_type(T), X, Y, msk; kwds...)

# NOTE: For offset-arrays, the in-place version `forge_mask!` must be used.
forge_mask(::Type{T}, X::AbstractVector, Y::AbstractVector, msk::Mask; kwds...) where {T} =
    forge_mask!(Matrix{T}(undef, length(X), length(Y)), X, Y, msk; kwds...)

"""
    TwoDimensional.forge_mask!(dst, [X, Y,] msk; kwds...) -> dst
    TwoDimensional.forge_mask!(dst, [X, Y,] elems...; kwds...) -> dst

In-place version of [`TwoDimensional.forge_mask`](@ref), it overwrites the destination
array `dst` with the mask and returns it. If coordinates `X` and `Y` along the axes of
`dst` are not specified, `(X, Y) = axes(dst)` is assumed.

"""
forge_mask!(A::AbstractMatrix, msk::Mask; kwds...) = forge_mask!(A, axes(A)..., msk; kwds...)

forge_mask!(A::AbstractMatrix, elems::MaskElement...; kwds...) =
    forge_mask!(A, Mask(elems); kwds...)

forge_mask!(A::AbstractMatrix, elems::List{<:MaskElement}; kwds...) =
    forge_mask!(A, Mask(elems); kwds...)

function forge_mask!(dst::AbstractMatrix, X::AbstractVector, Y::AbstractVector,
                     elems::MaskElement...; kwds...)
    return forge_mask!(dst, X, Y, Mask(elems); kwds...)
end

function forge_mask!(dst::AbstractMatrix, X::AbstractVector, Y::AbstractVector,
                     elems::List{<:MaskElement}; kwds...)
    return forge_mask!(dst, X, Y, Mask(elems); kwds...)
end

function forge_mask!(dst::AbstractMatrix, X::AbstractVector, Y::AbstractVector,
                     msk::Mask;
                     antialiasing::Integer = default_antialiasing,
                     opaque = zero(eltype(dst)),
                     transparent = oneunit(eltype(dst)),
                     multithreading::Bool = Threads.nthreads() > 1)
    # Check that coordinate units are compatible and determine a suitable unit for all
    # coordinates.
    T = float(promote_type(eltype(X), eltype(Y), coord_type(msk)))

    # Check arguments.
    I, J = axes(dst)
    axes(X) == (I,) || throw(DimensionMismatch("coordinates `X` have incompatible indices"))
    axes(Y) == (J,) || throw(DimensionMismatch("coordinates `Y` have incompatible indices"))
    antialiasing = as(Int , antialiasing)
    opaque = as(eltype(dst), opaque)
    transparent = as(eltype(dst), transparent)

    # Call the real method.
    return unsafe_forge_mask!(dst, convert_eltype(T, X), convert_eltype(T, Y),
                              convert_coord_type(T, msk);
                              antialiasing = antialiasing, opaque = opaque,
                              transparent = transparent,
                              multithreading = multithreading)
end

function unsafe_forge_mask!(dst::AbstractMatrix{V},
                            X::AbstractVector{T},
                            Y::AbstractVector{T},
                            msk::Mask{T};
                            antialiasing::Int,
                            opaque::V,
                            transparent::V,
                            multithreading::Bool) where {T,V}
    partial = interpolate(opaque, transparent, 1//2)
    δx = grid_step(X)
    δy = grid_step(Y)

    # Determine whether a cell of the mask is fully or partially opaque or transparent.
    if isempty(msk)
        # No objects, assume mask is transparent.
        fill!(mask, transparent)
    else
        initial = true
        for elem in msk
            unsafe_coarse_transmission!(dst, X, δx, Y, δy, elem,
                                        opaque, partial, transparent, initial)
            initial = false
        end
    end

    if antialiasing > 1
        # Refine mask at partially opaque cells.
        DX = subrange(antialiasing, δx)
        DY = subrange(antialiasing, δy)
        if multithreading && Threads.nthreads() > 1
            states = [Array{Bool}(undef, antialiasing, antialiasing) for k in 1:Threads.nthreads()]
            values = [Array{Tuple{Int,V}}(undef, 0) for k in 1:Threads.nthreads()]
            R = CartesianIndices(dst)
            # In the multi-threaded phase, only read `dst`.
            Threads.@threads for i in eachindex(IndexLinear(), dst)
                if dst[i] == partial
                    # Cell is partially opaque/transparent.
                    let k = Threads.threadid(), I = R[i], x = X[I[1]], y = Y[I[2]],
                        v = unsafe_sampled_transmission!(states[k], x, DX, y, DY,
                                                         msk, opaque, transparent)
                        push!(values[k], (i, v))
                    end
                end
            end
            # Override, coarse partial opaque values in `dst` by refined values computed
            # by the threads.
            for k in eachindex(values)
                for (i, v) in values[k]
                    dst[i] = v
                end
            end
        else
            state = Array{Bool}(undef, antialiasing, antialiasing)
            @inbounds for i in CartesianIndices(dst)
                if dst[i] == partial
                    # Cell is partially opaque/transparent.
                    x, y = X[i[1]], Y[i[2]]
                    dst[i] = unsafe_sampled_transmission!(state, x, DX, y, DY,
                                                          msk, opaque, transparent)
                end
            end
        end
    end

    return dst
end

"""
    TwoDimensional.unsafe_coarse_transmission!(dst, X, δx, Y, δy, elem,
                                               opaque, partial, transparent, initial)

updates 2-dimensional array `dst` with a coarse estimation of the transmission due to mask
element `elem`. `X` and `Y` give the coordinates of `dst` along its dimensions with
respective steps `δx` and `δy`. `opaque`, `partial`, and `transparent` are the three
possible values to store in `dst`. `initial` indicates whether `elem` is the first of the
mask elements.

The function is *unsafe* because it assumes without checking that `dst` is indexed as `X`
and `Y` along its first and second dimensions respectively.

"""
function unsafe_coarse_transmission!(dst::AbstractMatrix{V},
                                     X::AbstractVector{T}, δx::T,
                                     Y::AbstractVector{T}, δy::T,
                                     elem::MaskElement{T},
                                     opaque::V,
                                     partial::V,
                                     transparent::V,
                                     initial::Bool) where {T,V}
    if initial
        fill!(dst, is_opaque(elem) ? transparent : opaque)
    end
    return unsafe_coarse_transmission!(dst, X, δx, Y, δy, shape(elem),
                                       is_opaque(elem) ? opaque : transparent,
                                       partial)
end

"""
    TwoDimensional.unsafe_coarse_transmission!(dst, X, δx, Y, δy, obj, inside, partial)

updates 2-dimensional array `dst` with a coarse estimation of the transmission due to
geometrical shape object `obj`. `X` and `Y` give the coordinates of the cells of `dst`
along its dimensions with respective steps `δx` and `δy`. `inside` is the transmission for
cells fully inside the boundaries of `obj`, while `partial` is the transmission for cells
overlapping the boundaries of `obj`. Cells fully outside `obj` are left unchanged.

The function is *unsafe* because it assumes without checking that `dst` is indexed as `X`
and `Y` along its first and second dimensions respectively.

"""
function unsafe_coarse_transmission!(dst::AbstractMatrix{V},
                                     X::AbstractVector{T}, δx::T,
                                     Y::AbstractVector{T}, δy::T,
                                     obj::ShapeElement{T},
                                     inside::V,
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
            # FIXME We know that coordinates are in order so constructor could be faster...
            cell = Rectangle((x - δx/2, y - δy/2), (x + δx/2, y + δy/2))
            overlap = Overlap(cell, obj)
            if overlap == INSIDE
                dst[i,j] = inside
            elseif overlap == PARTIAL
                dst[i,j] = partial
            end
        end
    end
    return dst
end

"""
    TwoDimensional.unsafe_sampled_transmission!(state, x, sx, y, sy,
                                                msk, opaque, transparent) -> v

yields the transmission due to the elements in mask `msk` for a rectangular cell centered
at `(x,y)` by sampling the cell at offsets `sx` and `sy` along its dimensions. The
returned value `v` is interpolated in the range `[opaque,transparent]` of transmissions
from the fully opaque to the fully transparent parts of the mask.

Argument `state` is a workspace array used to sample the transmission. The function
is *unsafe* because it assumes without checking that `state` is indexed as `sx` and `sy`
along its first and second dimensions respectively.

"""
function unsafe_sampled_transmission!(state::AbstractMatrix{Bool},
                                      x::T, sx::AbstractVector{T},
                                      y::T, sy::AbstractVector{T},
                                      msk::Mask{T},
                                      opaque::V, transparent::V) where {T,V}
    X = x .+ sx
    Y = y .+ sy
    count = -1 # to trigger resetting of state array
    for elem in msk
        count = unsafe_sampled_transmission!(state, X, Y, elem, count)
    end
    return interpolate(opaque, transparent, max(count, 0)//length(state))
end

# Helper function to dispatch on the type of the mask element. Return the number of
# transparent samples. `count` is initially set to a negative value.
function unsafe_sampled_transmission!(state::AbstractMatrix{Bool},
                                      X::AbstractVector{T},
                                      Y::AbstractVector{T},
                                      obj::MaskElement{T},
                                      count::Int) where {T}
    # Operation can be much faster if cell is fully outside object boundaries. For complex
    # shaped object, it is cheaper to check whether the cell is outside the bounding-box
    # of the object.
    box = BoundingBox(obj)
    if box.xmax < first(X) || box.xmin > last(X) || box.ymax < first(Y) || box.ymin > last(Y)
        # The cell is outside the bounding-box of the mask element. The mask element has
        # thus no indicence unless it is the first of the stack.
        if count < zero(count)
            # This is the first mask element applied to this cell.
            if is_opaque(obj)
                fill!(state, true)
                count = oftype(count, length(state))
            else
                fill!(state, false)
                count = zero(count)
            end
        end
    else
        # The cell may overlap the boundaries of the object.
        I, J = axes(state)
        if count < zero(count)
            # This is the first mask element applied to this cell.
            count = zero(count)
            @inbounds for j in J
                y = Y[j]
                for i in I
                    x = X[i]
                    inside = Overlap(Point(x, y), obj) != OUTSIDE
                    transp = xor(inside, is_opaque(obj))
                    state[i,j] = transp
                    count += oftype(count, transp)
                end
            end
        elseif is_opaque(obj)
            @inbounds for j in J
                y = Y[j]
                for i in I
                    x = X[i]
                    if Overlap(Point(x, y), obj) != OUTSIDE
                        # Point is considered as being inside boundaries of opaque
                        # mask.
                        if state[i,j]
                            # Sub-cell previously counted as transparent.
                            state[i,j] = false
                            count -= oneunit(count)
                        end
                    end
                end
            end
        else # element is transparent
            @inbounds for j in J
                y = Y[j]
                for i in I
                    x = X[i]
                    if Overlap(Point(x, y), obj) != OUTSIDE
                        # Point is considered as being inside boundaries of transparent mask.
                        if !state[i,j]
                            # Sub-cell previously considered as opaque.
                            state[i,j] = true
                            count += oneunit(count)
                        end
                    end
                end
            end
        end
    end
    return count
end

"""
    TwoDimensional.grid_step(x::AbstractVector) -> stp

yields the step along a vector of coordinates, throwing an error if the increment between
successive values of `x` is not positive or not uniform.

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
    if len > 2
        tol = sqrt(eps(floating_point_type(stp)))
        stpmin = stp*(one(tol) - tol)
        stpmax = stp*(one(tol) + tol)
        flag = true
        @inbounds for i ∈ firstindex(x):lastindex(x)-1
            flag &= (stpmin ≤ x[i+1] - x[i] ≤ stpmax)
        end
        flag || throw(ArgumentError("vector of coordinates has non-uniform steps"))
    end
    return stp
end

"""
    TwoDimensional.interpolate(a, b, f) -> x

yields linearly interpolated value between `a` and `b` by a fraction `f`. If `f` is a
dimensionless factor, then the result is:

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

yields the overlapping of pixel `pxl` with shape object `obj`. `pxl` may be a point, a
bounding-box, or a rectangle. Returned value is:

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
    # Get coordinates of cell corners relative to circle center and quickly check whether
    # cell is certainly outside circle.
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
    # Most distant point in cell from center of circle is one of the corners. If this
    # point is inside the circle, then the cell is fully inside the circle.
    x²max = max(x0^2, x1^2)
    y²max = max(y0^2, y1^2)
    if x²max + y²max ≤ r²
        return INSIDE
    end
    # Otherwise, if closest point of one of the cell edge is inside the circle, overlap is
    # partial.
    if clamp(zero(T), x0, x1)^2 + y²max ≤ r²
        return PARTIAL
    end
    if x²max + clamp(zero(T), y0, y1)^2 ≤ r²
        return PARTIAL
    end
    # Otherwise, there is no overlap.
    return OUTSIDE
end

# FIXME: A point can only be inside or outside.
Overlap(point::Point, polygon::Polygon) = point ∈ polygon ? INSIDE : OUTSIDE

# FIXME: Check that the following algorithm is correct for a cell having zero width or
#        height.
function Overlap(cell::Rectangle{T}, polygon::Polygon{T}) where {T}
    (x0, y0), (x1, y1) = cell # retrieve coordinates of cell vertices
    V = vertices(polygon) # retrieve the list of vertices of the polygon

    # If any polygon edges (strictly) crosses one of the 4 edges of the cell, then there
    # is some partial overlapping.
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

    # Otherwise, the cell is either fully inside or fully outside the polygon. This can be
    # tested for any point of the cell, its center for example.
    return Overlap(Point((x0 + x1)/2, (y0 + y1)/2), polygon)
end
