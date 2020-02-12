mutable struct Dataset <: AbstractVector{Measurement}
	name::AbstractString
	data_dir_name::AbstractString
	measurements::Vector{Measurement}
	plotcolor::Union{String, Int, Symbol}

	Dataset() = new("unnamed", "unknown", Measurement[], :red)

	Dataset(name::AbstractString, data_dir_name::AbstractString, measurements::Array{Measurement, 1}) = new(name, data_dir_name, measurements, :red)
end

function Dataset(name::AbstractString, data_dir_name::AbstractString)::Dataset
	Dataset(name, data_dir_name, data_set_from_conv_data(data_dir_name))
end

function Dataset(name::AbstractString, data_dir_name::AbstractString, detector::Detector, daq::DAQ)::Dataset
	ms = data_set_from_conv_data(data_dir_name)
	for m in ms
		m.daq = daq
		m.detector = detector
	end
	Dataset(name, data_dir_name, ms)
end

function getindex(d::Dataset, i::Int)::Measurement
	return d.measurements[i]
end
function size(d::Dataset)
	return size(d.measurements)
end
function length(d::Dataset)
	return length(d.measurements)
end

function Base.push!(d::Dataset, m::Measurement)::Nothing
	push!(d.measurements, m)
end


function println(io::IO, d::Dataset)
	println(io, d.name)
	println(io, "Contains $(length(d.measurements)) measurements")
	for (i, m) in enumerate(d.measurements)
		println(i, "\t", m.motor_pos_r, "\t", m.motor_pos_phi, "\t", m.motor_pos_z)
	end
end
function print(io::IO, d::Dataset)
	print(io, d.name, " - Contains $(length(d.measurements)) measurements")
end
display(io::IO, d::Dataset) = println(io, d)
show(io::IO, d::Dataset) = println(io, d)

function show(io::IO, ::MIME"text/plain", d::Dataset)
    show(io, d)
end

function sort(d::Dataset; kwargs...)::Dataset
	ms_sorted = sort(d.measurements; kwargs...)
	new_ds = Dataset()
	new_ds.name = d.name
	new_ds.measurements = ms_sorted
	@show new_ds[1].motor_pos_phi
	new_ds.plotcolor = d.plotcolor
	return new_ds
end

function sort!(d::Dataset; kwargs...)::Nothing
	ms_sorted = sort(d.measurements; kwargs...);
	d.measurements[:] = ms_sorted[:];
	return nothing
end

function get_path_to_raw_data(d::Dataset)::AbstractString
	return joinpath(GAT.USER_DATA_PATH, d.data_dir_name[1:4], d.data_dir_name, "raw_data")
end