module TwoDimensionalOffsetArraysExt
if isdefined(Base, :get_extension)
    using TwoDimensional, OffsetArrays
else
    using ..TwoDimensional, ..OffsetArrays
end

using Base: OneTo
using TypeUtils
import TwoDimensional: new_array

new_array(::Type{T}, rngs::ArrayAxes{N}) where {T,N} =
    OffsetArray(Array{T}(undef, as_array_size(rngs)), rngs)

end # module
