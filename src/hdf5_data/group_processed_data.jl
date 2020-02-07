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


function get_pulses(fn::AbstractString, event_index::Int, sampling_rate::T, tdcs::Vector{T}, bl, bl_inv, decay_factors, c)::AbstractArray{T, 2} where {T <: Real}
    daq_pulses = get_daq_pulse(fn, event_index)
    tdc_pulses::Array{T, 2} = baseline_substraction_and_decay_correction(T.(daq_pulses), bl, bl_inv, decay_factors);
    return calibrate_pulses(tdc_pulses, c)
end
function get_pulses(m::Measurement, event_index::Int, file_number::Int=1)::AbstractArray{<:Real, 2}
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = Float32
    bl::Int = Int(m.daq.baseline_length)
    bl_inv::T = T(1 / bl)
    decay_window_start_index = m.daq.n_samples - m.daq.decay_window_length + 1
    decay_window_length_inv = T(1 / m.daq.decay_window_length)
    sampling_rate = T(m.daq.sampling_rate)
    tdcs = T.(read_analysis_result_dataset(m, "init_tau_decay_constants"))
    decay_factors::Vector{T} = T[exp( - (1 / sampling_rate) / (tdc * (10^(-6.0)) ) ) for tdc in tdcs]
    c::Array{T, 2} = T.(read_analysis_result_dataset(m, "calibration_matrix"))
    return get_pulses(inputfiles[file_number], event_index, sampling_rate, tdcs, bl, bl_inv, decay_factors, c)
end

function get_pulses(m::Measurement, energy_range::Interval; single_channel_idx::Int = 1, max_nevents::Int = 10, single_channel_energy_atol = 5.0)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    T = Float32
    bl::Int = Int(m.daq.baseline_length)
    bl_inv::T = T(1 / bl)
    decay_window_start_index = m.daq.n_samples - m.daq.decay_window_length + 1
    decay_window_length_inv = T(1 / m.daq.decay_window_length)
    sampling_rate = T(m.daq.sampling_rate)
    tdcs = T.(read_analysis_result_dataset(m, "init_tau_decay_constants"))
    decay_factors::Vector{T} = T[exp( - (1 / sampling_rate) / (tdc * (10^(-6.0)) ) ) for tdc in tdcs]
    c::Array{T, 2} = T.(read_analysis_result_dataset(m, "calibration_matrix"))
    max_n_pulses = [ Array{T, 2}(undef, m.daq.n_samples, length(tdcs)) for evt in 1:max_nevents ]
    n_pulses_selected::Int = 0
    T_energy_range = T(energy_range.left)..T(energy_range.right)
    T_single_channel_energy_atol = T(single_channel_energy_atol)
    for ifn in inputfiles
        h5open(ifn, "r") do h5f
            g_daq = g_open(h5f, "DAQ_Data")
            g_pd = g_open(h5f, "Processed_data")
            d_energies = d_open(g_pd, "energies")
            d_daq_pulses = d_open(g_daq, "daq_pulses")
            chunk_size::NTuple{3, Int} = get_chunk(d_daq_pulses)
            n_samples::Int, n_channel::Int, n_events::Int = size(d_daq_pulses)
            evt_ranges::Vector{UnitRange{Int}} = event_range_iterator(n_events, chunk_size[3])

            for evt_range in evt_ranges
                chunk_pulses::Array{T, 3} = d_daq_pulses[: , :, evt_range]
                chunk_energies::Array{T, 2} = d_energies[:, evt_range]
                for i in 1:length(evt_range)
                    evt_energies::Vector{T} = chunk_energies[:, i]
                    if evt_energies[1] in T_energy_range
                        if isapprox(evt_energies[1], evt_energies[single_channel_idx], atol = T_single_channel_energy_atol)
                            waveforms::Array{T, 2} = baseline_substraction_and_decay_correction(chunk_pulses[:, :, i], bl, bl_inv, decay_factors); 
                            calibrate_pulses!(waveforms, c)
                            n_pulses_selected += 1
                            max_n_pulses[n_pulses_selected] = waveforms
                        end
                    end
                    if n_pulses_selected == max_nevents break end
                end
                if n_pulses_selected == max_nevents break end
            end
        end
        if n_pulses_selected == max_nevents break end
    end
    @info n_pulses_selected
    return max_n_pulses
end

function get_event_flags(fn::AbstractString)::Vector{EventFlag}
    return h5read(fn, "Processed_data/event_flags")
end
function get_event_flags(m::Measurement)::Vector{EventFlag}
    inputfiles::Vector{AbstractString} = gather_absolute_paths_to_hdf5_input_files(m)
    n_events::Int = get_number_of_events(m)
    flags::Vector{EventFlag} = Vector{EventFlag}(undef, n_events)
    last_evt_idx::Int = 0
    for fn in inputfiles
        n::Int = get_number_of_events(fn)
        flags[last_evt_idx+1:last_evt_idx+n] = get_event_flags(fn)
        last_evt_idx += n
    end
    return flags
end