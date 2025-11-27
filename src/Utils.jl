"""
    vec(A)

Return a column-stacked vector of the entries of `A`, matching MATLAB's `(:)` behavior.
"""
function vec(A::AbstractArray)
    return reshape(copy(A), :)
end

"""
    unvec(v, m, n)

Reshape vector `v` into an `m × n` matrix using column-major ordering.
"""
function unvec(v::AbstractVector, m::Integer, n::Integer)
    length(v) == m * n || throw(ArgumentError("vector length does not match target dimensions"))
    return reshape(v, m, n)
end
