mktemp_custom(parent=tempdir(), format="tmpXXXXXX") = begin
    b = joinpath(parent, format)
    p = ccall(:mkstemp, Int32, (Ptr{UInt8},), b) # modifies b
    systemerror(:mktemp, p == -1)
    return (b, fdio(p, true))
end

function sis3316_get_nchannel(ifn::AbstractString)::Int
    input_io = open(CompressedFile(ifn), "r")
    n_channel::Int = 0
    for unsorted in eachchunk(input_io, SIS3316.UnsortedEvents)
        for (ch, events) in unsorted
            if !isempty(events) n_channel += 1 end
        end
        break   
    end
    return n_channel
end
function sis3316_get_nsamples(ifn::AbstractString; evt_merge_window::AbstractFloat = 100e-9)::Int
    input_io = open(CompressedFile(ifn), "r")
    n_channel::Int = 0
    n_samples::Int = 0
    for unsorted in eachchunk(input_io, SIS3316.UnsortedEvents)
        sorted = sortevents(unsorted, merge_window = evt_merge_window)
        evtv = Vector{Pair{Int64, SIS3316.RawChEvent}}()
        for evt in sorted
            n_samples = length(evt[1].samples)
            break
        # for (ch, events) in unsorted

        #     if !isempty(events) 
        #         display(events)
        #         n_channel += 1 
        #     end
        end
        break   
    end
    close(input_io)
    return n_samples
end

function sis3316_get_nchunks(ifn::AbstractString)::Int
    input_io = open(CompressedFile(ifn), "r")
    n_chunks::Int = 0
    for unsorted in eachchunk(input_io, SIS3316.UnsortedEvents)
        n_chunks += 1
    end
    return n_chunks
end

function get_conv_data_hdf5_filename(fn::AbstractString)::String
    ofn::String = ""
    if occursin("-adc", fn)
        if occursin(r"adc1-", fn) 
            tmpidx = match(r"adc1-", fn).offset
            ofn = join([fn[1:tmpidx-1],fn[tmpidx+5:end]])
        end
        if occursin(r"adc2-", fn) 
            tmpidx = match(r"adc2-", fn).offset
            ofn = join([fn[1:tmpidx-1],fn[tmpidx+5:end]])
        end
        if endswith(ofn, ".bz2")
            ofn = "$(ofn[1:end-4][1:first(findlast(".", ofn[1:end-4]))-1]).hdf5"
        else
            ofn = "$(ofn[1:first(findlast(".", ofn))-1]).hdf5"
        end
    else
        if endswith(fn, "bz2")
            ofn = "$(fn[1:end-4][1:first(findlast(".", fn[1:end-4]))-1]).hdf5"
        else
            ofn = "$(fn[1:first(findlast(".", fn))-1]).hdf5"
        end
    end
    return ofn
end

function sis3316_to_hdf5(ifn::AbstractString;   evt_merge_window::AbstractFloat = 100e-9, waveform_format = :none, 
                                                compress=true, overwrite=false, use_true_event_number=false,              
                                                chunk_n_events::Int=1000, waveform_type::DataType = Int32) 
    # if endswith(ifn, "bz2")
    #     ofn = "$(ifn[1:end-4][1:first(findlast(".", ifn[1:end-4]))-1]).hdf5"
    # else
    # 	ofn = "$(ifn[1:first(findlast(".", ifn))-1]).hdf5"
    # end
    ofn = get_conv_data_hdf5_filename(ifn)

    if overwrite || !isfile(ofn)
        n_channel::Int = sis3316_get_nchannel(ifn)
        n_samples::Int = sis3316_get_nsamples(ifn, evt_merge_window = evt_merge_window)
        @info "creating $ofn"
        @info "Detected $(n_channel) active channel in this file."
        @info "Detected $(n_samples) samples per channel in this file."
        output_tmpname, tmpio = mktemp_custom(pwd(), "$(ofn).tmp-XXXXXX")
        close(tmpio)
        input_io = open(CompressedFile(ifn), "r")
        h5f = h5open(output_tmpname, "w")
        sis3316_to_hdf5(input_io, h5f, n_channel=n_channel, n_samples_per_channel = n_samples; evt_merge_window = evt_merge_window, waveform_format = waveform_format, compress=compress, 
                                                            use_true_event_number=use_true_event_number, chunk_n_events=chunk_n_events, waveform_type = waveform_type)

        mv(output_tmpname, ofn, force = overwrite)
        run(`chmod g+rw $ofn`) # give read+write rights for the group
        @info "Single file $(ofn) done!"
    else
        @info "Output file \"$(ofn)\" already exists, skipping \"$(ifn)\"."
    end
    return ofn
end


function sis3316_to_hdf5(input_io::IO, output_hdf5_file;    n_channel=16, n_samples_per_channel = 5000, evt_merge_window::AbstractFloat = 100e-9, 
                                                            waveform_format = :none, compress=true, use_true_event_number=false,
                                                            chunk_n_events::Int=100, waveform_type::DataType = Int32)
  	_time(x::Pair{Int64, SIS3316.RawChEvent}) = time(x.second)

    info_idx = Ref{Int32}(0)
    info_time =  Ref{Float64}(0)

    raw_pp_ch = Int32[]
    raw_pp_mca = Int32[]
    raw_pp_trig_max = Int32[]
    raw_pp_peak_pos =  Int32[]
    raw_pp_peak_height =  Int32[]
    # raw_pp_acc = map(i -> bindings[Symbol("raw_pp_acc_$i")] = Int32[], 1:8)

    raw_trig_ch =  Int32[]
    raw_trig_trel = Float64[]
    raw_trig_pileup = Int32[]
    raw_trig_overflow =  Int32[]

    raw_wf_ch = Int32[]
    raw_wf_smpl_n = Int32[]
    raw_wf_smpl_v = Int32[]


    ch_sized_vecs = Vector[raw_pp_ch, raw_pp_mca, raw_pp_trig_max,
        raw_pp_peak_pos, raw_pp_peak_height,
        raw_trig_ch, raw_trig_trel, raw_trig_pileup, raw_trig_overflow
    ]

    for v in ch_sized_vecs sizehint!(v, n_channel) end  

    reader_1 = eachchunk(input_io, SIS3316.UnsortedEvents)

    evtno = 0
    number_of_corrupted_events = 0
    myevtno = 0

    #dirty hardcode but works for the moment
    n_max_events = -1 #from 300000
    start_n_events = 10000000
    daq_n_channels = n_channel
    daq_n_samples = n_samples_per_channel
    if daq_n_samples == 0
        @warn "Detected no pulse samples in input files. Setting `waveform_format` to `:none`"
        waveform_format = :none
    end
    @info "chunk_n_events: $chunk_n_events"

    chunk_size_pulses = chunk_n_events
    chunk_size_others = 5000

    g_daq = g_create(output_hdf5_file, "DAQ_Data")
    d_event_number  = d_create(g_daq, "event_number", Int32, ((start_n_events,),(n_max_events,)), "chunk", (chunk_size_others,))
    d_daq_time      = d_create(g_daq, "daq_time",   Float64, ((1, start_n_events,),(1, n_max_events,)), "chunk", (1, chunk_size_others,))
    d_daq_energy    = d_create(g_daq, "daq_energies", waveform_type, ((daq_n_channels,start_n_events),(daq_n_channels,n_max_events)), "chunk", (daq_n_channels,chunk_size_others))
    if waveform_format == :integers
        if compress
            d_daq_pulses = d_create(g_daq, "daq_pulses", Float32, ((daq_n_samples,daq_n_channels,start_n_events),(daq_n_samples,daq_n_channels,n_max_events)), 
                                            "chunk", (daq_n_samples,daq_n_channels,chunk_size_pulses), "shuffle", (), "deflate", 3  )
        else
            d_daq_pulses = d_create(g_daq, "daq_pulses", Int32, ((daq_n_samples,daq_n_channels,start_n_events),(daq_n_samples,daq_n_channels,n_max_events)), "chunk", (daq_n_samples,daq_n_channels,chunk_size_pulses) )
        end
    end

    temp_event_number::Array{Int32, 1}  = zeros(Int32, chunk_n_events)
    temp_daq_energy::Array{Float64, 2}  = zeros(Int32, daq_n_channels, chunk_n_events)
    temp_daq_time::Array{Float64, 2} = zeros(Float64, 1, chunk_n_events)
    temp_daq_pulses::Array{waveform_type, 3} = zeros(waveform_type, daq_n_samples, daq_n_channels, chunk_n_events)
    
    i_sub_buffer::Int = 1
    evt_sub_no::Int = 0

    chunk_i = 0    
    for unsorted in reader_1
    	chunk_i += 1
        # if chunk_i > 1500 break end # debugging
    # 	if chunk_i % 100 == 0
			 # @info "Now processing chunk: $(chunk_i)"
    # 	end
        sorted = sortevents(unsorted, merge_window = evt_merge_window)
        evtv = Vector{Pair{Int64, SIS3316.RawChEvent}}()
        timestamps = Vector{Float64}()

        energynull = SIS3316.EnergyValues(0, 0)
        mawnull    = SIS3316.MAWValues(0, 0, 0)
        psanull    = SIS3316.PSAValue(0, 0)
        flagsnull  = SIS3316.EvtFlags(false, false, false, false)

        for evt in sorted
            evtno += 1
            if evtno%10000==0
              @info "Now prcessing event: $(evtno+1)..."
            end

            resize!(evtv, length(evt))
            evtv[:] = [Pair(key, evt[key]) for key in keys(evt)]
            sort!(evtv, by = first)
            resize!(timestamps, length(evtv))
            map!(_time, timestamps, evtv)
            starttime = isempty(timestamps) ? zero(Float64) : minimum(timestamps)

            info_idx.x = evtno
            info_time.x = starttime

            for v in ch_sized_vecs empty!(v) end

            empty!(raw_wf_ch)
            empty!(raw_wf_smpl_n)
            empty!(raw_wf_smpl_v)


            for (ch, chevt) in evtv
                push!(raw_pp_ch, ch)
                push!(raw_pp_mca, chevt.energy.maximum)

                
                push!(raw_pp_trig_max, chevt.trig_maw.maximum)
                push!(raw_pp_peak_pos, chevt.peak_height.index)
                push!(raw_pp_peak_height, chevt.peak_height.value)


                push!(raw_trig_ch, ch)
                push!(raw_trig_trel, time(chevt) - starttime)
                push!(raw_trig_pileup, chevt.pileup_flag + 2 * chevt.flags.pileup +  4 * chevt.flags.repileup)
                push!(raw_trig_overflow, 1 * chevt.flags.overflow +  2 * chevt.flags.underflow)

                if !isempty(chevt.samples)
                    if waveform_format == :integers
                        push!(raw_wf_ch, ch)
                        push!(raw_wf_smpl_n, length(chevt.samples))
                        append!(raw_wf_smpl_v, chevt.samples)
                    end
                end
            end

            myevtno = evtno-number_of_corrupted_events

           
            if length(raw_pp_mca) == daq_n_channels  # == 5
                evt_sub_no += 1
                # d_daq_energy[:,myevtno] = raw_pp_mca[:]
                temp_daq_energy[:,evt_sub_no] = raw_pp_mca[:]
                if use_true_event_number
                    # d_event_number[1,myevtno] = evtno
                    # d_event_number[myevtno] = evtno
                    temp_event_number[evt_sub_no] = evtno
                else
                    # d_event_number[1,myevtno] = myevtno
                    # d_event_number[myevtno] = myevtno
                    temp_event_number[evt_sub_no]= myevtno
                end
                # d_daq_time[1,myevtno] = starttime
                # d_daq_time[myevtno] = starttime
                temp_daq_time[:,evt_sub_no] .= starttime
				#println("test1")
                if waveform_format == :integers
                    if length(raw_wf_smpl_v) == daq_n_samples*daq_n_channels  # ==25000
                        @inbounds for ich in 1:daq_n_channels # == 5
                            temp_daq_pulses[:,ich,evt_sub_no] = waveform_type.(raw_wf_smpl_v[(ich-1)*daq_n_samples+1:ich*daq_n_samples])
                        end
                    end
                end

                if evt_sub_no == chunk_n_events # now sub buffer is full -> write to file
                    evt_range = (i_sub_buffer-1)*evt_sub_no+1:(i_sub_buffer)*evt_sub_no
                    d_event_number[evt_range] = temp_event_number
                    d_daq_time[:, evt_range] = temp_daq_time
                    d_daq_energy[:,evt_range] = temp_daq_energy
                    if waveform_format == :integers
                        d_daq_pulses[:,:,evt_range] = temp_daq_pulses
                    end
                    i_sub_buffer += 1
                    evt_sub_no = 0
                    temp_event_number = zeros(Int32, chunk_n_events)
                    temp_daq_energy = zeros(Int32, daq_n_channels, chunk_n_events)
                    temp_daq_time = zeros(Float64, 1, chunk_n_events)
                    temp_daq_pulses = zeros(waveform_type, daq_n_samples, daq_n_channels, chunk_n_events)
                end
            else
                number_of_corrupted_events+=1
            end
            # @info "evt_sub_no: $evt_sub_no"
        end
  #       # if (i_sub_buffer == 3 &&  evt_sub_no == 315) break  end
    end

    # now write last unfull buffer to file
    idxzero = findfirst(iszero, temp_event_number)
    if idxzero > 1
        tmp_range = 1:idxzero-1
        evt_range = (i_sub_buffer-1)*chunk_n_events+1:(i_sub_buffer-1)*chunk_n_events+1+idxzero-1-1
        d_event_number[evt_range] = temp_event_number[tmp_range]
        d_daq_time[:, evt_range] = temp_daq_time[tmp_range]
        d_daq_energy[:,evt_range] = temp_daq_energy[:,tmp_range]
        if waveform_format == :integers
            d_daq_pulses[:,:,evt_range] = temp_daq_pulses[:,:,tmp_range]
        end
    end
    n_events_written = (i_sub_buffer-1)*chunk_n_events+1+idxzero-1-1

    # set_dims!(d_event_number, (1,myevtno))
    # set_dims!(d_daq_time,     (1,myevtno))
    @info n_events_written
    set_dims!(d_event_number, (n_events_written,))
    set_dims!(d_daq_time,     (size(d_daq_time, 1), n_events_written))
    set_dims!(d_daq_energy,   (daq_n_channels, n_events_written))
    if waveform_format == :integers
        set_dims!(d_daq_pulses, (daq_n_samples,daq_n_channels,n_events_written))
    end

    g_general = g_create(output_hdf5_file, "INFO")
    attrs(g_general)["N_events_in_binary_file"] = evtno
    attrs(g_general)["N_events"] = n_events_written
    attrs(g_general)["N_corrupted_events"] = number_of_corrupted_events

    close(output_hdf5_file)
    return nothing
end
