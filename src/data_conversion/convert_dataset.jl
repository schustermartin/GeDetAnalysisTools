function convert_data_files_sis3316_to_hdf5( raw_dir = pwd();
                                             overwrite=false, evt_merge_window::AbstractFloat=100e-9, waveform_format=:integers,
                                             chunk_n_events::Int=100, keep_individual_hdf5_files::Bool = false, compress_raw_data::Bool = true,
                                             waveform_type::DataType = Int32)
    current_dir = pwd()
    cd(raw_dir)
    if !isdir("../conv_data") mkdir("../conv_data") end

    all_files = readdir(raw_dir)

    if any(fn -> occursin("-adc", fn), all_files)
        @info "Detected multi STRUCK dataset"
        twostrucks_convert_all_data_files_in_raw_data_folder(raw_dir, overwrite = overwrite, evt_merge_window = evt_merge_window, waveform_format = waveform_format,
                                                                waveform_type = waveform_type, chunk_n_events = chunk_n_events, keep_individual_hdf5_files = keep_individual_hdf5_files,
                                                                compress_raw_data = compress_raw_data )
    else
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
                    @info "Now on $(myid()): Compression dat files. Skipping..."
                    # if (file_is_compressed && isfile(fn)) rm(fn) else compress_file(fn, keep_input_files=false) end
                end

                @info "Now on $(myid()): Finished with $(ofn)."
            else
                @info "Skipping $fn. Already converted."
            end
        end

        pmap( process_file, sis_files )
    end
    # run(`chmod -Rf ug+rw ../`)
    cd(current_dir)
    return nothing
end


function convert_data_files_sis3316_to_hdf5( d::Dataset;    overwrite=false, evt_merge_window::AbstractFloat=100e-9, waveform_format=:integers,
                                                            chunk_n_events::Int=100, keep_individual_hdf5_files::Bool = false, compress_raw_data::Bool = true,
                                                            waveform_type::DataType = Int32)
    convert_data_files_sis3316_to_hdf5( get_path_to_raw_data(d),   overwrite = overwrite, evt_merge_window = evt_merge_window, waveform_format = waveform_format,
                                                                    waveform_type = waveform_type, chunk_n_events = chunk_n_events, keep_individual_hdf5_files = keep_individual_hdf5_files,
                                                                    compress_raw_data = compress_raw_data )
end


export convert_all_files
function convert_all_files(;directory = pwd(),process_pulses::Bool=true,compressed::Bool=false,compress_after_conversion::Bool=true,move_converted_file::Bool=true,use_true_event_number::Bool=false, new_pulse_format=true)
  if compressed==false
   list_of_files_to_convert = filter(x->endswith(x,".dat"),readdir(directory))
  # elseif compressed==false
  #  list_of_files_to_convert = filter(x->( startswith(x,"bkg_") && endswith(x,"Z.dat")), readdir(directory))
  elseif compressed ==true
   list_of_files_to_convert = filter(x->endswith(x,"-raw.dat.bz2"),readdir(directory))
  end
  @info("Total $(length(list_of_files_to_convert)) files to convert..")

 function process_one_file(file)
    if compressed==true
      run(`bzip2 -d $file`)
      file = file[1:searchindex(file,match(r"bz2",file).match)-2]
    end
    output_filename = "$(file[1:end-4]).hdf5"
    # cd("$directory")
    if !ispath("../conv_data/$output_filename")
      @info("now processing file $file")
      if process_pulses == true
        sis3316_to_hdf5(file,evt_merge_window=100e-9,waveform_format = :integers,use_true_event_number=use_true_event_number)
        # sis3316_to_hdf5(file,evt_merge_window=100e-9,n_channel=5,waveform_format = :integers,use_true_event_number=use_true_event_number, new_pulse_format=new_pulse_format)
      elseif process_pulses == false
        sis3316_to_hdf5(file,n_channel=5,use_true_event_number=use_true_event_number, new_pulse_format=new_pulse_format)
      end
      if move_converted_file==true
        @info("moving $output_filename...")
        mv(output_filename, "../conv_data/$output_filename")
      end
      if compress_after_conversion==true
        @info("compressing binary file... ")
        run(`bzip2 $file`)
        @info("finished compression")
      end
    else
      @info("file $output_filename already exists")
      nothing
    end
   end

   map(process_one_file,list_of_files_to_convert)

  if basename(pwd())=="raw_data"
      chmod.("../conv_data/".*filter(x -> endswith(x, ".hdf5"), readdir("../conv_data/")), 0o774)
  end

end

function convert_data_files_sis3316_to_hdf5( d::Dataset; kwargs...)
    convert_data_files_sis3316_to_hdf5( get_path_to_raw_data(d); kwargs... )
end
