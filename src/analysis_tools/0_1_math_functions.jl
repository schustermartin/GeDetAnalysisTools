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
