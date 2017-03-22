
function compensated_summation{T}(x::T, y::T, e::T)
    local err::T
    local res::T

    err = e + x
    res = err + y
    e   = err + (y - res)

    return (res, e)
end


function istriustrict{T}(A::Matrix{T})
    m, n = size(A)
    @inbounds for j = 1:min(n,m-1), i = j:m
        if A[i,j] ≠ zero(T)
            return false
        end
    end
    return true
end

function istrilstrict{T}(A::Matrix{T})
    m, n = size(A)
    @inbounds for j = 2:n, i = 1:min(j,m)
        if A[i,j] ≠ zero(T)
            return false
        end
    end
    return true
end

function L2norm(x)
    local l2::eltype(x) = 0
    for xᵢ in x
        l2 += xᵢ^2
    end
    l2
end

function l2norm(x)
    sqrt(L2norm(x))
end

function simd_scale!(x, a)
    @inbounds for i=1:length(x)
        x[i] *= a
    end
    nothing
end

function simd_copy!{T}(x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] = x[i]
    end
    nothing
end

function simd_copy!{T}(x::Matrix{T}, y::Matrix{T})
    @assert size(x) == size(y)
    @inbounds for j=1:size(y,2)
        for i=1:size(y,1)
            y[i,j] = x[i,j]
        end
    end
    nothing
end

function simd_copy_scale!(a, x, y)
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] = a*x[i]
    end
    nothing
end

"Copy the first dimension of a 2D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:length(x)
        x[i] = y[i,j]
    end
    nothing
end

"Copy the first dimension of a 3D array y into a 1D array x."
function simd_copy_xy_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:length(x)
        x[i] = y[i,j,k]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 2D array y."
function simd_copy_yx_first!(x, y, j)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:size(y, 1)
        y[i,j] = x[i]
    end
    nothing
end

"Copy a 1D array x into the second dimension of a 2D array y."
function simd_copy_yx_second!(x, y, j)
    @assert length(x) == size(y, 2)
    @inbounds for i=1:size(y, 2)
        y[j,i] = x[i]
    end
    nothing
end

"Copy a 1D array x into the first dimension of a 3D array y."
function simd_copy_yx_first!(x, y, j, k)
    @assert length(x) == size(y, 1)
    @inbounds for i=1:size(y, 1)
        y[i,j,k] = x[i]
    end
    nothing
end

function simd_copy_yx_first_last!(x, y, j)
    @assert size(x, 1) == size(y, 1)
    @assert size(x, 2) == size(y, 3)
    @inbounds for k=1:size(y, 3)
        for i=1:size(y, 1)
            y[i,j,k] = x[i,k]
        end
    end
    nothing
end

function simd_xpy!{T}(x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] += x[i]
    end
    nothing
end

function simd_xpy!{T}(x::Matrix{T}, y::Matrix{T})
    @assert size(x) == size(y)
    @inbounds for j=1:size(y,2)
        for i=1:size(y,1)
            y[i,j] += x[i,j]
        end
    end
    nothing
end

function simd_axpy!{T}(a, x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y)
    @inbounds for i=1:length(y)
        y[i] += a*x[i]
    end
    nothing
end

function simd_axpy!{T}(a, x::Vector{T}, y::Vector{T}, e::Vector{T})
    @assert length(x) == length(y) == length(e)

    local err::eltype(e)
    local ty::eltype(y)

    @inbounds for i=1:length(y)
        err = e[i] + a*x[i]
        ty  = y[i]
        y[i] = err + ty
        e[i] = err + (ty - y[i])
    end
    nothing
end

function simd_axpy!{T}(a, x::Matrix{T}, y::Matrix{T})
    @assert size(x) == size(y)
    @inbounds for j=1:size(y,2)
        for i=1:size(y,1)
            y[i,j] += a*x[i,j]
        end
    end
    nothing
end

function simd_wxpy!{T}(w::Vector{T}, x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y) == length(w)
    @inbounds for i=1:length(w)
        w[i] = x[i] + y[i]
    end
    nothing
end

function simd_wxpy!{T}(w::Matrix{T}, x::Matrix{T}, y::Matrix{T})
    @assert size(w) == size(x) == size(y)
    @inbounds for j=1:size(y,2)
        for i=1:size(y,1)
            w[i,j] = x[i,j] + y[i,j]
        end
    end
    nothing
end

function simd_waxpy!{T}(w::Vector{T}, a, x::Vector{T}, y::Vector{T})
    @assert length(x) == length(y) == length(w)
    @inbounds for i=1:length(w)
        w[i] = a*x[i] + y[i]
    end
    nothing
end

function simd_waxpy!{T}(w::Matrix{T}, a, x::Matrix{T}, y::Matrix{T})
    @assert size(w) == size(x) == size(y)
    @inbounds for j=1:size(y,2)
        for i=1:size(y,1)
            w[i,j] = a*x[i,j] + y[i,j]
        end
    end
    nothing
end

function simd_aXbpy!{T}(a, b::Vector{T}, X::Matrix{T}, y::Vector{T})
    @assert length(y) == size(X, 1)
    @assert length(b) == size(X, 2)
    local ty::T
    @inbounds for i=1:length(y)
        ty = zero(T)
        for j=1:length(b)
            ty += X[i,j] * b[j]
        end
        y[i] = a*ty
    end
    nothing
end

function simd_abXpy!{T}(a::T, b::Vector{T}, X::Matrix{T}, y::Vector{T})
    @assert length(y) == size(X, 2)
    @assert length(b) == size(X, 1)
    local ty::T
    @inbounds for i=1:length(y)
        ty = zero(T)
        for j=1:length(b)
            ty += b[j] * X[j,i]
        end
        y[i] += a*ty
    end
    nothing
end

function simd_mult!{T}(w::Vector{T}, X::Matrix{T}, y::Vector{T})
    @assert length(w) == size(X, 1)
    @assert length(y) == size(X, 2)
    local tw::T
    @inbounds for i=1:length(w)
        tw = zero(T)
        for j=1:length(y)
            tw += X[i,j] * y[j]
        end
        w[i] = tw
    end
    nothing
end

function simd_mult!{T}(w::Vector{T}, y::Vector{T}, X::Matrix{T})
    @assert length(w) == size(X, 2)
    @assert length(y) == size(X, 1)
    local tw::T
    @inbounds for i=1:length(w)
        tw = 0
        for j=1:length(y)
            tw += y[j] * X[j,i]
        end
        w[i] = tw
    end
    nothing
end

function simd_mult!{T}(W::Matrix{T}, X::Matrix{T}, Y::Matrix{T})
    @assert size(W,1) == size(X, 1)
    @assert size(W,2) == size(Y, 2)
    @assert size(X,2) == size(Y, 1)
    local tw::T
    @inbounds for j in 1:size(W,2)
        for i in 1:size(W,1)
            tw = 0
            for k in 1:size(X,2)
                tw += X[i,k] * Y[k,j]
            end
            W[i,j] = tw
        end
    end
    nothing
end
