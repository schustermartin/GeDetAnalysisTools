function convert_data_files_sis3316_to_hdf5( raw_dir = pwd(); 
                                             overwrite=false, evt_merge_window::AbstractFloat=100e-9, waveform_format=:integers, 
                                             chunk_n_events::Int=100, keep_individual_hdf5_files::Bool = false, compress_raw_data::Bool = true, 
                                             waveform_type::DataType = Int32)
    current_dir = pwd()
    cd(raw_dir)
    if !isdir("../conv_data") mkdir("../conv_data") end

    all_files = readdir(raw_dir)

    sis_files = filter(x -> occursin(".dat", x), all_files)

    function process_file( fn::AbstractString )
        cd(raw_dir)
        ofn = joinpath("../conv_data", get_conv_data_hdf5_filename(fn))
        @show ofn
        if !isfile(ofn) || overwrite
            file_is_compressed::Bool = endswith(fn, ".bz2")
            if file_is_compressed
                @info "Now on $(myid()): Decompressing file `fn`"
                file_is_compressed = true
                decompress_file(fn, overwrite = true, keep_input_files = true)
            end

            @info "Now on $(myid()): converting to $(ofn))"
            ofn_tmp = sis3316_to_hdf5(  fn,
                                        evt_merge_window=evt_merge_window,
                                        waveform_format=waveform_format,
                                        overwrite=overwrite,
                                        chunk_n_events=chunk_n_events,
                                        waveform_type = waveform_type)
            mv(ofn_tmp, ofn, force = true )

            if compress_raw_data
                @info "Now on $(myid()): Compression dat files."
                if (file_is_compressed && isfile(fn)) rm(fn) else compress_file(fn, keep_input_files=false) end
            end

            @info "Now on $(myid()): Finished with $(ofn)."
        else
            @info "Skipping $fn. Already converted."
        end
    end

    pmap( process_file, sis_files )

    # run(`chmod -Rf ug+rw ../`)

    cd(current_dir)
    return nothing
end
