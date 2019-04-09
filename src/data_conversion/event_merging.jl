# using Plots; pyplot()
@everywhere using GeDetAnalysisTools, HDF5

raw_dir = pwd()
all_files = readdir(raw_dir)
adc1_files = filter(x -> occursin("adc1", x), all_files)
adc2_files = filter(x -> occursin("adc2", x), all_files)
adc1_files = filter(x -> endswith(x, ".hdf5"), adc1_files)
adc2_files = filter(x -> endswith(x, ".hdf5"), adc2_files)

@info raw_dir
@sync @distributed for i in eachindex(adc1_files)
	ifn1 = adc1_files[i]
	ifn2 = adc2_files[i]
	ofn = joinpath("../conv_data", GAT.get_conv_data_hdf5_filename(ifn1))
	GeDetAnalysisTools.combine_two_hdf5_files(ifn1, ifn2, ofn, overwrite = true, chunk_n_events = 100 )
end