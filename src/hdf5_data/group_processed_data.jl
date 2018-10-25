function get_energies(fn::AbstractString)::Array{<:Real,2}
    h = h5open(fn, "r+")
    g = g_open(h,"Processed_data")
    d = d_open(g,"energies")
    energies = read(d)
    close(h)
    return energies
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