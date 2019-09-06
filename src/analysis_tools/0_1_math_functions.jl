function gauss(x, p)
    scale = p[1]
    σ     = p[2]
    μ     = p[3]
    return @fastmath @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2))
end


function linear_polynom(x, p)
	offset = p[1]
	slope  = p[2]
	return @fastmath @. offset + slope * x
end


function flat_offset(x, p) 
	return @. p[1]
end

function gauss_plus_first_order_polynom(x, p)
    scale = p[1]
    σ    = p[2]
    μ     = p[3]
    cp0   = p[4] 
    cp1   = p[5]
    if scale < 0 || σ <= 0 
        return -Inf
    else
        return @fastmath @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
    end
end


function gauss_plus_fixed_first_order_polynom(x, p; offset = 0.0, slope = 0.0)
    scale = p[1]
    σ     = p[2]
    μ     = p[3]
    cp0   = offset
    cp1   = slope
    if scale < 0 || σ < 0 
        return 0
    else
        return @fastmath @. scale * exp(-0.5 * ((x - μ)^2) / (σ^2)) / (sqrt(2 * π * σ^2)) + cp0 + cp1 * (x - μ)
    end
end


function exponential_decay(x, p)
    return @. p[1] * exp( -x / p[2])
end    


function scaled_cauchy(x, p)
    scale = p[1]
    σ = p[2]
    μ  = p[3]
    return @. scale * σ / (π * (σ^2 + (x - μ)^2))
end


function linear_function_fixed_offset_at_zero(x, p)
    return @. p[1] * x
end
