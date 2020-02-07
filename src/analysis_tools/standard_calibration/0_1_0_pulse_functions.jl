@inline @fastmath function fastsum(arr::AbstractArray{T, 1})::T where {T}
	sum = zero(T)
	@inbounds @simd for i in eachindex(arr)
		sum += arr[i]
	end
	return sum
end

@inline @fastmath function fastmean(arr::AbstractArray{T, 1}, inv_length::T)::T where {T<:AbstractFloat}
	sum = fastsum(arr)
	return sum * inv_length
end

@inline @fastmath function fastmean(arr::AbstractArray{T, 1})::AbstractFloat where {T<:AbstractFloat}
	sum = zero(T)
	@inbounds @simd for i in eachindex(arr)
		sum += arr[i]
	end
	return sum / length(arr)
end

@fastmath function baseline_substraction( pulse::AbstractArray{T, 1}, baseline_length::Int, inv_length::T)::AbstractArray{T, 1} where {T<:AbstractFloat}
	@inbounds return pulse .- fastmean( pulse[1:baseline_length], inv_length)
end
@fastmath function baseline_substraction!(pulse::AbstractArray{T, 1}, baseline_length::Int, inv_length::T)::Nothing where {T<:AbstractFloat}
	@inbounds pulse[:] = pulse .- fastmean( pulse[1:baseline_length], inv_length)
	return nothing
end

function baseline_substraction(pulses::AbstractArray{T, 2}, baseline_length::Int, inv_length::T)::AbstractArray{T, 2} where {T<:AbstractFloat}
	rp = Array{T, 2}(undef, size(pulses)) #zeros(T, size(pulses))
	@inbounds for ichn in axes(pulses, 2)
		rp[:, ichn] = baseline_substraction(pulses[:, ichn], baseline_length, inv_length)
	end
	return rp
end
function baseline_substraction!(pulses::AbstractArray{T, 2}, baseline_length::Int, inv_length::T)::Nothing where {T<:AbstractFloat}
	@inbounds for ichn in axes(pulses, 2)
		pulses[:, ichn] = baseline_substraction(pulses[:, ichn], baseline_length, inv_length)
	end
	return nothing
end


@fastmath function decay_correction(pulse::AbstractArray{T, 1}, decay_factor::T)::AbstractArray{T, 1} where {T<:AbstractFloat}
	@inbounds begin
		rp = Array{T, 1}(undef, length(pulse))
		rp[1] = pulse[1]
		for i in 2:length(rp)
			rp[i] = pulse[i] + rp[i-1] - pulse[i-1] * decay_factor
		end
		return rp
	end
end
@fastmath function decay_correction!(pulse::AbstractArray{T, 1}, decay_factor::T)::Nothing where {T<:AbstractFloat}
	@inbounds pulse[:] = decay_correction(pulse, decay_factor)
	return nothing
end


@fastmath function decay_correction(pulses::AbstractArray{T, 2}, decay_factors::Array{T, 1})::AbstractArray{T, 2} where {T<:AbstractFloat}
	@inbounds begin
		rp = Array{T, 2}(undef, size(pulses)...)
		for chn in eachindex(decay_factors)
			rp[1, chn] = pulses[1, chn]
			for i in 2:size(rp, 1)
				rp[i, chn] = pulses[i, chn] + rp[i-1, chn] - pulses[i-1, chn] * decay_factors[chn]
			end
		end
		return rp
	end
end
@fastmath function decay_correction!(pulses::AbstractArray{T, 2}, decay_factors::Array{T, 1})::Nothing where {T<:AbstractFloat}
	@inbounds for chn in eachindex(decay_factors)
		pulses[:, chn] = decay_correction(pulses[:, chn], decay_factors[chn])
	end
	return nothing
end


@fastmath function baseline_substraction_and_decay_correction(pulse::AbstractArray{T, 1}, baseline_length::Int, inv_length::T, decay_factor::T)::AbstractArray{T, 1} where {T<:AbstractFloat}
    return decay_correction(baseline_substraction(pulse, baseline_length, inv_length), decay_factor)
end
@fastmath function baseline_substraction_and_decay_correction!(pulse::AbstractArray{T, 1}, baseline_length::Int, inv_length::T, decay_factor::T)::Nothing where {T<:AbstractFloat}
	baseline_substraction!(pulse, baseline_length, inv_length)
    decay_correction!(pulse, decay_factor)
end


@fastmath function baseline_substraction_and_decay_correction(pulses::AbstractArray{T, 2}, baseline_length::Int, inv_length::T, decay_factors::Array{T, 1})::AbstractArray{T, 2} where {T<:AbstractFloat}
	return decay_correction(baseline_substraction(pulses, baseline_length, inv_length), decay_factors)	
end
@fastmath function baseline_substraction_and_decay_correction!(pulses::AbstractArray{T, 2}, baseline_length::Int, inv_length::T, decay_factors::Array{T, 1})::Nothing where {T<:AbstractFloat}
	baseline_substraction!(pulses, baseline_length, inv_length)
	decay_correction!(pulses, decay_factors)
end


@inline @fastmath function calibrate_pulses(pulses::AbstractArray{T, 2}, c::AbstractArray{T, 2})::AbstractArray{T, 2} where {T}
	@inbounds pulses * c 
end
@fastmath function calibrate_pulses!(pulses::AbstractArray{T, 2}, c::AbstractArray{T, 2})::Nothing where {T}
	@inbounds pulses[:] = calibrate_pulses(pulses, c)
	return nothing
end


