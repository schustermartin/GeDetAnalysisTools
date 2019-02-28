function one_boundary(x, par)
    A = 1.0#par[1] # amplitude
    s = par[2] #par[2] # `slope` 
    b = par[3] # boundary position
    o = 0.0#par[4] # Offset
    return @. 0.5 * A * erf( s * (x - b) ) + o + A * 0.5
end

function two_boundaries(x, par)
    A = 1.0#par[1] # amplitude
    s = par[2] # `slope` 
    bl = par[3] # boundary position
    br = par[4] # boundary position
    o = 0.0 #par[5] # Offset
    return @. A * 0.5 * ( erf( s * (x - bl) ) + erf( -s * (x - br))) + o 
end

