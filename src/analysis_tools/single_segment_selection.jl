@inline @fastmath function ss_event_rel(ce::T, se::T, d::T)::Bool where {T<:Real}
	abs(ce - se) * 2 / (ce + se) < d ? true : false
end
@inline @fastmath function ss_event_abs(ce::T, se::T, d::T)::Bool where {T<:Real}
	abs(ce - se) < d ? true : false
end

"""
	get_single_segment_channel_index_rel(energies::Array{T, 1}, d::T)::UInt8 where {T<:Real}

d::T = allowed relative difference in unit one.
"""
@fastmath function get_single_segment_channel_index_rel(energies::Array{T, 1}, d::T)::UInt8 where {T<:Real}
	@inbounds begin
		idx::UInt8 = 0
		counter::UInt8 = 0
		l::UInt8 = length(energies)
		for ichn in UInt8(2):l
			if ss_event_rel(energies[1], energies[ichn], d)
				idx = ichn
				counter += one(UInt8)
			end
		end
		counter == one(UInt8) ? idx : zero(UInt8)
	end
end


@fastmath function get_single_segment_channel_index_abs(energies::Array{T, 1}, d::T, d_high_energy::T = d)::UInt8 where {T<:Real}
	@inbounds begin
		idx::UInt8 = 0
		counter::UInt8 = 0
		l::UInt8 = length(energies)
		for ichn in UInt8(2):l
			if ss_event_abs(energies[1], energies[ichn], energies[1] < T(5000) ? d : d_high_energy)
				idx = ichn
				counter += one(UInt8)
			end
		end
		counter == one(UInt8) ? idx : zero(UInt8)
	end
end