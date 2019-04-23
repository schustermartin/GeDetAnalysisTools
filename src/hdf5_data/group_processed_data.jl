function get_energies(fn::AbstractString)::Array{<:Real,2}
    return h5read(fn, "Processed_data/energies")
end
function get_energies(m::Measurement)::Array{<:Real,2}
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = GeDetAnalysisTools.get_eltype_of_dataset(inputfiles[1], "Processed_data", "energies")
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    energies = Array{T, 2}(undef, n_channel, n_total_events)
    last_idx::Int = 0
    @inbounds for i in eachindex(inputfiles)
        en = get_energies(inputfiles[i])
        n_new_events::Int = size(en, 2)
        energies[:, last_idx+1:last_idx+n_new_events] = en
        last_idx += n_new_events
    end
    return energies
end

function get_measured_pulse_amplitudes(fn::AbstractString)::Array{<:Real,2}
    return h5read(fn, "Processed_data/measured_pulse_amplitudes")
end

function get_measured_pulse_amplitudes(m::Measurement)::Array{<:Real,2}
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = GeDetAnalysisTools.get_eltype_of_dataset(inputfiles[1], "Processed_data", "measured_pulse_amplitudes")
    n_total_events = get_number_of_events(m)
    n_channel = get_number_of_channel(m)
    mpas = Array{T, 2}(undef, n_channel, n_total_events)
    last_idx::Int = 0
    @inbounds for i in eachindex(inputfiles)
        en = get_measured_pulse_amplitudes(inputfiles[i])
        n_new_events::Int = size(en, 2)
        mpas[:, last_idx+1:last_idx+n_new_events] = en
        last_idx += n_new_events
    end
    return mpas
end

function get_single_segment_indices(fn::AbstractString)
    return h5read(fn,"Processed_data/single_segment_indices")
end
function get_single_segment_indices(m::Measurement)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    return vcat([get_single_segment_indices(f) for f in inputfiles]...)
end
