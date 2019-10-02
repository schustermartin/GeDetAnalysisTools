function gauss(x::T, p::Array{T, 1})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2))
end

function gauss(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ gauss(v, p) for v in x ]
end



function gauss_plus_cut_exp(x::T, p::Vector{T}) ::T where T
    exp_scale::T = p[1] # Left Skew Ampl
    exp_τ::T = p[2] # Left
    Scale::T = p[3] # Amplitude of Gaussian
    σ::T = p[4]
    μ::T = p[5]
    if x<μ
        myexp = exp_scale * exp((x-(μ-σ))/exp_τ)
    else
        myexp = 0.0
    end
    return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + myexp
end
function gauss_plus_cut_exp(x::Vector{T}, p::Vector{T}) ::Vector{T} where T
    return map(ix->gauss_plus_cut_exp(ix,p), x)
end


function gauss_plus_left_skew(x::T, p::Vector{T}) ::T where T
    LSScale::T = p[1] # Left Skew Ampl
    LSSlope::T = p[2] # Left
    Scale::T = p[3] # Amplitude of Gaussian
    σ::T = p[4]
    μ::T = p[5]
    return @fastmath Scale .* sqrt.(π) ./ 2  .* σ .* LSScale .* exp( (σ ./ (2 .*LSScale) )^2 .+ (x .- μ) ./ LSSlope ) .* erfc.(σ ./ (2 .* LSSlope) .+ (x .- μ) ./ σ )
end
function emg(x::T, p::Vector{T}) ::T where T
    scale::T = p[1] # Amplitude of Gaussian
    σ::T = p[2]
    μ::T = p[3]
    τ::T =p[4]# exponential decay
    return @fastmath scale .* σ ./ τ .* sqrt.(π ./ 2) .* exp( (1 ./ 2) .* (σ ./ τ) )^2 .- (x .- μ) ./ τ ) .* erfc.( (1 ./ sqrt.(2) ) .* (σ ./ τ .- (x .- μ) ./ σ ) )
end
function emg_right(x::Vector{T}, p::Vector{T}) ::Vector{T} where T
    scale::T = p[1] # Amplitude of Gaussian
    σ::T = p[2]
    μ::T = p[3]
    τ::T =p[4]# exponential decay
    return @fastmath scale .* σ ./ τ .* sqrt.(π ./ 2) .* exp( (1 ./ 2) .* (σ ./ τ) )^2 .- (x .- μ) ./ τ ) .* erfc.( (1 ./ sqrt.(2) ) .* (σ ./ τ .- (x .- μ) ./ σ ) )
end



function linear_polynom(x::T, p::Array{T, 1})::T where {T}
	offset::T = p[1]
	slope::T  = p[2]
	return @fastmath offset + slope * x
end
function linear_polynom(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
	return T[ linear_polynom(v, p) for v in x ]
end

function flat_offset(x::T, p::Array{T, 1})::T where {T}
	offset::T = p[1]
	return @fastmath offset
end
function flat_offset(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
	return T[ flat_offset(v, p) for v in x ]
end


function gauss_plus_first_order_polynom(x::T, p::Array{T, 1})::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    cp0::T   = p[4]
    cp1::T   = p[5]
    if scale < 0 || σ <= 0
        return -T(Inf)
    else
        return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
    end
end

function gauss_plus_first_order_polynom(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
    return T[ gauss_plus_first_order_polynom(v, p) for v in x ]
end


function gauss_plus_fixed_first_order_polynom(x::T, p::Array{T, 1}; offset = 0.0, slope = 0.0)::T where {T}
    scale::T = p[1]
    σ::T     = p[2]
    μ::T     = p[3]
    cp0::T   = offset
    cp1::T   = slope
    if scale < 0 || σ < 0
        return 0
    else
        return @fastmath scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
    end
end
function gauss_plus_fixed_first_order_polynom(x::Array{T, 1}, p::Array{T, 1}; offset = 0.0, slope = 0.0)::Array{T, 1} where {T}
    return T[ gauss_plus_first_order_polynom(v, p, offset = offset, slope = slope) for v in x ]
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
