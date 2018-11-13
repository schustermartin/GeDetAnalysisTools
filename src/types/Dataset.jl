mutable struct Dataset
	name::AbstractString
	measurements::Array{Measurement, 1}

	Dataset() = new("unnamed", Measurement[])

	Dataset(name::AbstractString, measurements::Array{Measurement, 1}) = new(name, measurements)
end

function Dataset(name::AbstractString, data_dir_name::AbstractString)::Dataset
	Dataset(name, data_set_from_conv_data(data_dir_name))
end

function Dataset(name::AbstractString, data_dir_name::AbstractString, detector::Detector, daq::DAQ)::Dataset
	ms = data_set_from_conv_data(data_dir_name)
	for m in ms
		m.daq = daq
		m.detector = detector
	end
	Dataset(name, ms)
end

function Base.getindex(d::Dataset, i::Int)::Measurement
	return d.measurements[i]
end
function Base.size(d::Dataset)
	return size(d.measurements)
end
function Base.length(d::Dataset)
	return length(d.measurements)
end
function Base.iterate(d::Dataset, state=1)
 	if state > length(d)
 		return nothing
 	else
 		return d[state]
 	end
end

function Base.push!(d::Dataset, m::Measurement)::Nothing
	push!(d.measurements, m)
end


function Base.println(io::IO, d::Dataset)
	println(io, d.name)
	println(io, "Contains $(length(d.measurements)) measurements")
	for (i, m) in enumerate(d.measurements)
		println(i, "\t", m.motor_pos_r, "\t", m.motor_pos_phi, "\t", m.motor_pos_z)
	end
end
function Base.print(io::IO, d::Dataset)
	print(io, d.name, " - Contains $(length(d.measurements)) measurements")
end
Base.display(io::IO, d::Dataset) = println(io, d)
Base.show(io::IO, d::Dataset) = println(io, d)
