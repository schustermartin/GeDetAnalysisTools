function gauss(x::T, p::Array{T, 1})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2))
end

function gauss(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ gauss(v, p) for v in x ]
end

function linear_polynom(x::T, p::Array{T, 1})::T where {T}
	offset::T = p[1]
	slope::T  = p[2]
	return @fastmath offset + slope * x
end
function linear_polynom(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
	return T[ linear_polynom(v, p) for v in x ]
end

function gauss_plus_first_order_polynom(x::T, p::Array{T, 1})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    cp0::T   = p[4] 
    cp1::T   = p[5]
    return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
end

function gauss_plus_first_order_polynom(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ gauss_plus_first_order_polynom(v, p) for v in x ]
end


function exponential_decay(x::T, p::Array{T, 1})::T where {T}
    return p[1] * exp( -x / p[2])
end    
function exponential_decay(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ exponential_decay(v, p) for v in x ]
end

function scaled_cauchy(x::T, p::Array{T, 1})::T where {T}
    scale = p[1]
    σ = p[2]
    μ  = p[3]
    return scale * σ / (π * (σ^2 + (x - μ)^2))
end
function scaled_cauchy(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ scaled_cauchy(v, p) for v in x ]
end


function linear_function_fixed_offset_at_zero(x::T, p::Array{T, 1})::T where {T}
    return p[1] * x
end
function linear_function_fixed_offset_at_zero(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ linear_function_fixed_offset_at_zero(v, p) for v in x ]
end