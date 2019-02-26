function one_boundary(x, par)
    A = par[1] # amplitude
    s = par[2] # `slope` 
    b = par[3] # boundary position
    o = par[4] # Offset
    return @. 0.5 * A * tanh( s * (x - b) ) + o + A * 0.5
end

function two_boundaries(x, par)
    A = par[1] # amplitude
    s = par[2] # `slope` 
    bl = par[3] # boundary position
    br = par[4] # boundary position
    o = par[5] # Offset
    return @. A * 0.5 * ( tanh( s * (x - bl) ) + tanh( -s * (x - br))) + o 
end

