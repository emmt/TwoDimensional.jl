#
# basics.jl --
#
# Utilities for Julia interface to TAO.
#
#-------------------------------------------------------------------------------
#
# This file if part of the TAO software (https://github.com/emmt/TAO) licensed
# under the MIT license.
#
# Copyright (C) 2019, Éric Thiébaut.
#

"""

`dimensions(A)` yields the list of dimensions of array `A`.  The result is the
same as `size(A)` but it is also checked that `A` has standard indices
(starting at 1).

"""
function dimensions(A::AbstractArray{<:Any,N}) where {N}
    # The fastest code is LazyAlgebra.has_standard_indexing(A)
    LazyAlgebra.has_standard_indexing(A) ||
        error("array has non-standard indexing")
    return size(A)
end
