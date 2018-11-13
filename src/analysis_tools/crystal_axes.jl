function crystal_axis_osc_func(x::T, par::Array{T, 1})::T where {T<:Real} 
	return par[1] * sind(4 * (x - par[2])) + par[3]
end
function crystal_axis_osc_func(x::Array{T, 1}, par::Array{T, 1})::Array{T, 1} where {T<:Real}  
	return T[ crystal_axis_osc_func(v, par) for v in x ]
end

risetime_crystal_axes_dependancy = MFunction("Rise time crystal axes dependancy", 
														 [20.0, 0.0, 250.],
														 ["Amplitude", "Shift", "Offset"],
														 crystal_axis_osc_func)
