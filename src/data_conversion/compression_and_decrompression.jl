function compress_file(filename::AbstractString; overwrite=false, keep_input_files=false)::Nothing
    isfile(filename) ? nothing : error("File not found: '$filename'")
    current_dir = pwd()
    if isabspath(filename) cd(dirname(filename)) end
    options=""
    if overwrite
        keep_input_files ? options=` -fk` : options=` -f`
    else
        keep_input_files ? options=` -k` : options=``
    end
    if !isfile("$filename.bz2") || overwrite
        length(options) > 0 ? run(`pbzip2 $options $filename`) : run(`pbzip2 $filename`)
    elseif !keep_input_files
        rm(filename)
    end
    cd(current_dir)
    return nothing
end

function compress_files(filenames::Array{<:AbstractString}=readdir(); overwrite=false, keep_input_files=false)::Nothing
    for ubf in filenames
        compress_file(ubf, overwrite=overwrite, keep_input_files=keep_input_files)
    end
    return nothing
end

function decompress_file(filename::AbstractString; overwrite=false, keep_input_files=false)::Nothing
    endswith(filename, ".bz2") ? nothing : error("File is not a '.bz2' file.")
    isfile(filename) ? nothing : error("File not found: '$filename'")
    current_dir = pwd()
    if isabspath(filename) cd(dirname(filename)) end
    options=""
    if overwrite
        keep_input_files ? options=`-dfk` : options=`-df`
    else
        keep_input_files ? options=`-dk` : options=`-d`
    end
    if ((!in(basename(filename)[1:end-4], readdir())) || overwrite)
        @info "decompress file $filename"
        length(options) > 0 ? run(`pbzip2 $options $filename`) : run(`pbzip2 $filename`)
    end
    cd(current_dir)
    return nothing
end

function decompress_files(filenames::Array{<:AbstractString}=readdir(); overwrite=false, keep_input_files=true)::Nothing
    for cbf in filenames
        decompress_file(cbf, overwrite=overwrite, keep_input_files=keep_input_files)
    end
    return nothing
end