"""
	The field `f` should have two methods: 
	- f(x::T, p::Array{T, 1})::T where {T}
	- f(x::Array{T, 1}, p::Array{T, 1})::Array{T, 1} where {T}
"""
mutable struct MFunction{T, N} 
	name::AbstractString
	par::AbstractArray{T, N}
	par_names::AbstractArray{<:AbstractString, N}
	f::Function

	function MFunction(name::AbstractString, par::AbstractArray{T, N}, par_names::AbstractArray{<:AbstractString, N}, f::Function) where {T, N} 
		new{eltype(par), ndims(par)}(name, par, par_names, f)
	end
end

function Base.println(io::IO, f::MFunction)
	println("Function: $(f.name)")
	for i in eachindex(f.par)
		println("\t$(f.par_names[i]): $(f.par[i])")
	end
end
