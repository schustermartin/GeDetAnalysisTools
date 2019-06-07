function get_matching_index(dt::Vector{T}, dt1::T, dt2::T, dt3::T; acceptance::T = 5e-8)::Int where {T <: AbstractFloat}
    @fastmath @inbounds for i in 1:length(dt) - 2
        rms::T = abs(dt[i] - dt1) + abs(dt[i + 1] - dt2) + abs(dt[i + 2] - dt3)
        if rms <= acceptance
            return i  
        end
    end
    return -1
end

function merge(t1::Vector{T}, t2::Vector{T}) where {T <: AbstractFloat}
    t1 .-= t1[1]
    t2 .-= t2[1]
    Δt1 = diff(t1)
    Δt2 = diff(t2)

    f1_evt_indices = Int[]
    f2_evt_indices = Int[]

    j = -1
    i = 1
    no_match_in_this_chunk = false
    while j < 0 && !no_match_in_this_chunk
        j = get_matching_index(Δt2, Δt1[i], Δt1[i + 1], Δt1[i + 2])
        if j > 0 
            push!(f1_evt_indices, i)
            push!(f2_evt_indices, j)
        else
            i += 1     
            if i > length(Δt1)
                no_match_in_this_chunk = true
            end
        end
    end
    @show j, i

    if !no_match_in_this_chunk
        @showprogress 1 "Event merging..." for i in 2:length(Δt1) - 2
            j_last::Int = f2_evt_indices[end]
            j = @inbounds get_matching_index(Δt2[j_last + 1:end], Δt1[i:i+2]...)
            if j > 0
                push!(f1_evt_indices, i)
                push!(f2_evt_indices, j + j_last)
            end
        end
    end
    
    return f1_evt_indices, f2_evt_indices
end

function combine_two_hdf5_files(fn1::AbstractString, fn2::AbstractString, ofn::AbstractString; overwrite::Bool=false, chunk_n_events::Int = 100)
    h5i1 = h5open(fn1, "r")
    h5i2 = h5open(fn2, "r")
    g1 = g_open(h5i1, "DAQ_Data")
    g2 = g_open(h5i2, "DAQ_Data")
    d_event_number_1 = d_open(g1, "event_number")
    d_event_number_2 = d_open(g2, "event_number")
    d_daq_energy_1 = d_open(g1, "daq_energies")
    d_daq_energy_2 = d_open(g2, "daq_energies")
    d_daq_time_1 = d_open(g1, "daq_time")
    d_daq_time_2 = d_open(g2, "daq_time")

    waveform_format = exists(g1, "daq_pulses") ? :integers : :none
    if waveform_format == :integers
        d_daq_pulses_1 = d_open(g1, "daq_pulses")
        d_daq_pulses_2 = d_open(g2, "daq_pulses")
    end

    g_general_1 = g_open(h5i1, "INFO")
    g_general_2 = g_open(h5i2, "INFO")
    a_total_n_events_1 = a_open(g_general_1, "N_events_in_binary_file")
    a_total_n_events_2 = a_open(g_general_2, "N_events_in_binary_file")
    @info "File 1: $(read(a_total_n_events_1)) total events"
    @info "File 2: $(read(a_total_n_events_2)) total events"
    a_N_events_1 = a_open(g_general_1, "N_events")
    a_N_events_2 = a_open(g_general_2, "N_events")
    @info "File 1: $(read(a_N_events_1)) written events"
    @info "File 2: $(read(a_N_events_2)) written events"
    a_N_corrupted_events_1 = a_open(g_general_1, "N_corrupted_events")
    a_N_corrupted_events_2 = a_open(g_general_2, "N_corrupted_events")
    @info "File 1: $(read(a_N_corrupted_events_1)) corrupted events"
    @info "File 2: $(read(a_N_corrupted_events_2)) corrupted events"

    n_events_1 = size(d_event_number_1, 1)
    n_events_2 = size(d_event_number_2, 1)
    n_max_events = min(n_events_1, n_events_2)
    daq_n_channels_1 = size(d_daq_energy_1, 1)
    daq_n_channels_2 = size(d_daq_energy_2, 1)
    daq_n_channels = daq_n_channels_1 + daq_n_channels_2
    if waveform_format == :integers
        daq_n_samples = size(d_daq_pulses_1, 1)
    end

    t1 = read(d_daq_time_1)[1, :]
    t2 = read(d_daq_time_2)[1, :]
    T::Type = eltype(t1)

    n1=length(t1)
    n2=length(t2)
    nmax = min(n1, n2)
       
    f1_evt_indices, f2_evt_indices = merge(t1, t2)
    
	n_matches = length(f1_evt_indices)
    # n_matches = 2000 # debugging
    @info "$(n_matches) matches"
	@info "$(round(1 - n_matches / nmax, sigdigits = 4) * 100) % loss due event merging"

    @info "Now writing them to new file:"

    chunk_size_pulses = chunk_n_events
    chunk_size_others = 5000

    chunk_n_events = n_matches < chunk_n_events ? n_matches : chunk_n_events
    chunk_size_others = n_matches < chunk_size_others ? n_matches : chunk_size_others
    h5o = h5open(ofn, "w")
    try
        g_daq = g_create(h5o, "DAQ_Data")
        d_event_number  = d_create(g_daq, "event_number", Int32, ((n_matches,),(n_matches,)), "chunk", (chunk_size_others,))
        d_daq_time      = d_create(g_daq, "daq_time",   T, ((2, n_matches),(2,n_matches)), "chunk", (2, chunk_size_others))
        d_daq_energy    = d_create(g_daq, "daq_energies", Int32, ((daq_n_channels,n_matches),(daq_n_channels,n_matches)), "chunk", (daq_n_channels,chunk_size_others))
        if waveform_format == :integers
            d_daq_pulses    = d_create(g_daq, "daq_pulses", Int32, ((daq_n_samples,daq_n_channels,n_matches),(daq_n_samples,daq_n_channels,n_matches)), "chunk", (daq_n_samples,daq_n_channels,chunk_size_pulses), "shuffle", (), "deflate", 3 )
        end

        buf_iterator = 1:chunk_n_events:n_matches
        prog = ProgressMeter.Progress(n_matches, 1)
        for ievt in buf_iterator
            evt_range = ievt:ievt+chunk_n_events-1
            if ievt+chunk_n_events-1 > n_matches
                evt_range = ievt:n_matches
            end
            evt_range_1 = f1_evt_indices[evt_range]
            evt_range_2 = f2_evt_indices[evt_range]

            daq_energies_1 = d_daq_energy_1[:, evt_range]
            daq_energies_2 = d_daq_energy_2[:, evt_range]

            if waveform_format == :integers
                daq_pulses_1 = d_daq_pulses_1[:, : , evt_range]
                daq_pulses_2 = d_daq_pulses_2[:, : , evt_range]
            end

            buf_daq_times::Array{T, 2} = zeros(T, 2, length(evt_range))
            buf_daq_times[1, :] = t1[evt_range_1]
            buf_daq_times[2, :] = t2[evt_range_2]

            evt_range_1_chunk = evt_range_1 .- (first(evt_range) - 1)
            evt_range_2_chunk = evt_range_2 .- (first(evt_range) - 1)

            buf_daq_energies::Array{Int32, 2} = zeros(Int32, daq_n_channels, length(evt_range))
            if waveform_format == :integers
                buf_daq_pulses::Array{Int32, 3} = zeros(Int32, daq_n_samples, daq_n_channels, length(evt_range))
            end
            l = length(evt_range)
            for i in eachindex(evt_range) # this is slow because indicidual reading for each event
                inchunk = evt_range_1_chunk[i] <= l && evt_range_2_chunk[i] <= l
                if inchunk
                    buf_daq_energies[:, i] = vcat(daq_energies_1[:, evt_range_1_chunk[i]], daq_energies_2[:, evt_range_2_chunk[i]])
                else
                    buf_daq_energies[:, i] = vcat(d_daq_energy_1[:, evt_range_1[i]], d_daq_energy_2[:, evt_range_2[i]])
                end
                if waveform_format == :integers
                    if inchunk
                        buf_daq_pulses[:, :, i] = hcat( daq_pulses_1[:, :, evt_range_1_chunk[i]], daq_pulses_2[:, :, evt_range_2_chunk[i]] )
                    else
                        buf_daq_pulses[:, :, i] = hcat( d_daq_pulses_1[:, :, evt_range_1[i]], d_daq_pulses_2[:, :, evt_range_2[i]] )
                    end
                end
                next!(prog)
            end

            d_event_number[evt_range]   = collect(evt_range)
            d_daq_time[:, evt_range]    = buf_daq_times
            d_daq_energy[:, evt_range]  = buf_daq_energies
            if waveform_format == :integers
                d_daq_pulses[:, :, evt_range] = buf_daq_pulses
            end

            # next!(prog)
        end
        finish!(prog)


        g_general = g_create(h5o, "INFO")
        attrs(g_general)["N_events_in_binary_file"] = [read(a_total_n_events_1),read(a_total_n_events_2)]
        attrs(g_general)["N_events"] = n_matches
        attrs(g_general)["N_corrupted_events"] = read(a_total_n_events_1)-n_matches


        close(h5o)
        close(h5i1)
    catch err
        close(h5o)
        close(h5i1)
        close(h5i2)
        error(err)
    end

    run(`chmod g+rw $ofn`)
    return ofn
end

function two_sis3316_to_hdf5(fn1::AbstractString, fn2::AbstractString; evt_merge_window::AbstractFloat = 100e-9, waveform_format = :waveform_format,
                                                                        overwrite = false, chunk_n_events::Int=1000, keep_individual_hdf5_files::Bool=false, waveform_type::DataType = Int32)
    reg_adc_unit = r"adc1-"
    occursin(reg_adc_unit, fn1) ? nothing : error("filename '$fn1' does not contain 'adc1-'")
    tmpidx = match(reg_adc_unit, fn1).offset
    ofn = join([fn1[1:tmpidx-1],fn1[tmpidx+5:end]])
    if endswith(ofn, "bz2")
        ofn = "$(ofn[1:end-4][1:first(findlast(".", ofn[1:end-4]))-1]).hdf5"
    else
        ofn = "$(ofn[1:first(findlast(".", ofn))-1]).hdf5"
    end
    if !overwrite && isfile(ofn)
        @info "Already converted. Skipping $(ofn)."
        return ofn
    end
    @info "convert adc1 file:"
    ofn1 = sis3316_to_hdf5(fn1, evt_merge_window=evt_merge_window, waveform_format=waveform_format, compress=false, overwrite=overwrite, use_true_event_number=true, chunk_n_events=chunk_n_events, waveform_type = waveform_type)
    @info "convert adc2 file:"
    ofn2 = sis3316_to_hdf5(fn2, evt_merge_window=evt_merge_window, waveform_format=waveform_format, compress=false, overwrite=overwrite, use_true_event_number=true, chunk_n_events=chunk_n_events, waveform_type = waveform_type)
    @info "merge them:"
    combine_two_hdf5_files(ofn1, ofn2, ofn, overwrite=overwrite, chunk_n_events=chunk_n_events)
    if !keep_individual_hdf5_files
        rm(ofn1)
        rm(ofn2)
    end
    @info "merging done"
    return ofn
end

function twostrucks_convert_all_data_files_in_raw_data_folder(raw_dir=pwd(); overwrite=false, evt_merge_window::AbstractFloat=100e-9, 
                                        waveform_format=:integers, chunk_n_events::Int=100, keep_individual_hdf5_files::Bool = false, compress_raw_data::Bool = true,
                                        waveform_type::DataType = Int32)
    current_dir = pwd()
    cd(raw_dir)

    if !isdir("../conv_data") mkdir("../conv_data") end

    all_files = readdir(raw_dir)
    adc1_files = filter(x -> occursin("adc1", x), all_files)
    adc2_files = filter(x -> occursin("adc2", x), all_files)

    fn1_dat_files = String[]
    fn2_dat_files = String[]
    for fn1 in adc1_files
        endswith(fn1, ".bz2") ? tmp_name = fn1[1:end-4] : tmp_name = fn1
        if !in(tmp_name, fn1_dat_files) push!(fn1_dat_files, tmp_name) end
    end
    for fn2 in adc2_files
        endswith(fn2, ".bz2") ? tmp_name = fn2[1:end-4] : tmp_name = fn2
        if !in(tmp_name, fn2_dat_files) push!(fn2_dat_files, tmp_name) end
    end
    fn1_dat_files = filter(x->endswith(x, ".dat"), fn1_dat_files)
    fn2_dat_files = filter(x->endswith(x, ".dat"), fn2_dat_files)

    if length(fn1_dat_files) != length(fn2_dat_files)
        error("Different number of adc1 ($(length(fn1_dat_files))) and adc2 ($(length(fn2_dat_files))) files.")
    end

    function process_file_idx(i)
        cd(raw_dir)
        ifn1 = fn1_dat_files[i]
        ifn2 = fn2_dat_files[i]
        ofn = joinpath("../conv_data", get_conv_data_hdf5_filename(ifn1))
        if !isfile(ofn) || overwrite
            f1_is_compressed::Bool = false
            f2_is_compressed::Bool = false
            if !isfile(ifn1) && isfile(ifn1 * ".bz2")
                @info "Now on $(myid()): Decompressing file `$(ifn1 * ".bz2")`"
                f1_is_compressed = true
                decompress_file(ifn1 * ".bz2", overwrite=true, keep_input_files=true)
            end
            if !isfile(ifn2) && isfile(ifn2 * ".bz2")
                @info "Now on $(myid()): Decompressing file `$(ifn2 * ".bz2")`"
                f2_is_compressed = true
                decompress_file(ifn2 * ".bz2", overwrite=true, keep_input_files=true)
            end

            @info "Now on $(myid()): converting to $(ofn))"
            ofn_tmp = two_sis3316_to_hdf5(  ifn1, ifn2,
                                            evt_merge_window=evt_merge_window,
                                            waveform_format=waveform_format,
                                            waveform_type = waveform_type,
                                            overwrite=overwrite,
                                            chunk_n_events=chunk_n_events,
                                            keep_individual_hdf5_files=keep_individual_hdf5_files)
            mv(ofn_tmp, ofn, force = true )

            if compress_raw_data
                @info "Now on $(myid()): Compression dat files."
                if (f1_is_compressed && isfile(ifn1)) rm(ifn1) else compress_file(ifn1, keep_input_files=false) end
                if (f2_is_compressed && isfile(ifn1)) rm(ifn2) else compress_file(ifn2, keep_input_files=false) end
            end

            @info "Now on $(myid()): Finished with $(ofn)."
        else
            @info "Skipping $ifn1. Already converted."
        end
    end

    pmap( process_file_idx,  1:length(fn1_dat_files) )

    # run(`chmod -R ug+rw ../`)

    cd(current_dir)
    return nothing
end
