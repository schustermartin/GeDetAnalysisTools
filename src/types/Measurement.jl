mutable struct Measurement
	name::AbstractString
	data_set_name::AbstractString
	path_to_raw_data::AbstractString
	path_to_conv_data::AbstractString
	results_path::AbstractString
	datetime::DateTime
	side_source::AbstractString
	top_source::AbstractString
	outside_source::AbstractString
	motor_pos_z::Float64
	motor_pos_r::Float64
	motor_pos_phi::Float64
	sidesourcefacingsegment::Int
	topsourcefacingsegment::Int
	detector::Detector
	measurement_duration::Int
	daq::DAQ
	daq_livetime::Float64
	temperature::Float64
	pressure::Float64
	new_data_structure::Bool
	r::Function
	phi_side::Function
	phi_top::Function
	z::Function
	plotcolor::Union{String, Int, Symbol}

	function Measurement()
		new("unknown",  # name
			"unknown",  # data_set_name
			"unknown",  # path_to_raw_data
			"unknown",  # path_to_raw_data
			"unknown",  # results_path
			DateTime(0),  # date: 0000-01-01T00:00:00
			"unknown",  # side_source
			"unknown",  # top_source
			"unknown",  # outside_source
			NaN64,      # motor_pos_z
			NaN64,      # motor_pos_r
			NaN64,      # motor_pos_phi
			-1,         # sidesourcefacingsegment
			-1,         # topsourcefacingsegment
			Detector(), # Detector
			-1,         # measurement_duration
			DAQ(),      # daq
			-1,         # daq_livetime
			NaN64,      # temperature
			NaN64,      # pressure
			true,
			identity,
			identity,
			identity,
			identity,
			:red,
			)
	end

	function Base.show(io::IO, m::Measurement)
		for n in fieldnames(Measurement)
			println(io, "$n: $(getfield(m,n))")
		end
	end
end

import Base.println
import Base.print
import Base.show
import Base.display
function println(io::IO, dataset::Array{Measurement, 1})
	if length(dataset) > 0
		m = dataset[1]
		println(io, m.data_set_name)
		println(io, "contains $(length(dataset)) measurements")
	else
		"empty dataset"
	end
end
print(io::IO, dataset::Array{Measurement, 1}) = println(io, dataset)
show(io::IO, dataset::Array{Measurement, 1}) = println(io, dataset)
display(io::IO, dataset::Array{Measurement, 1}) = println(io, dataset)


r(m::Measurement)   = m.r(m.motor_pos_r)
phi_side(m::Measurement) = m.phi_side(m.motor_pos_phi)
phi_top(m::Measurement) = m.phi_top(m.motor_pos_phi)
z(m::Measurement)   = m.z(m.motor_pos_z)

############# Obtaining Information of a Measurement ##########
# new_data_structure -> Struck data storage structure: no subfolders in raw_data folders
# true -> new / Struck
# false -> OLD / PIXI
function get_measurement_name_from_compressed_data_file_name(filename::AbstractString; new_data_structure=true)
	if new_data_structure == true
		# measurement_name = filename
		measurement_name = filename[1:match(r"-\d{8}T\d{6}Z.*",filename).offset-1]
		return measurement_name
	else
		if endswith(filename, ".bin.bz2")
			measurement_name = filename[1:end-12]
		else
			error("File is not a compressed data file.")
		end
		return measurement_name
	end
end

function get_datetime_from_measurement_name(m::Measurement; new_data_structure=true)::DateTime
	files = filter(x -> startswith(x, m.name), readdir(m.path_to_raw_data))
	dt = if new_data_structure == true
		df = DateFormat("yyyymmddTHHMMSSZ")
		dts = match(r"-\d{8}T\d{6}Z.*", files[1]).match
		dts = dts[2:end-4]
		if endswith(dts, ".dat")
			dts = dts[1:end-4]
		end
		if endswith(dts, "-raw")
			dts = dts[1:end-4]
		end
		DateTime(dts, df)
	else
		df = DateFormat("yyyymmdd")
		ds = files[1][1:8]
		DateTime(ds, df)
	end
	return dt
end

function get_pressure_from_measurement_name(m::Measurement)
	reg = r"_p_\d{1}.\d{1}e-*\d{0,3}mbar"
	if occursin(reg, m.name)
		mt = match(reg, m.name)
		pressure = parse(Float64, mt.match[4:end-4])
		return pressure
	else
		return NaN
	end
end

function get_motor_pos_r_from_measurement_name(m::Measurement; new_data_structure=true)
	motor_pos_r = NaN64
	if new_data_structure == true
		reg = r"R_[0-9]*.{0,1}[0-9]*mm"
		if occursin(reg, m.name)
			motor_pos_r = parse(Float64, match( reg, m.name ).match[3:end-2])
		end
	else
		reg = r"pos_[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*"
		if occursin(reg, m.name)
			pos_string = match( reg, m.name ).match
			value_reg = r"[0-9]*[.][0-9]*"
			r_string = matchall(value_reg, pos_string)[2]
			motor_pos_r = parse(Float64, r_string)
		end
	end
	return motor_pos_r
end

function get_motor_pos_z_from_measurement_name(m::Measurement; new_data_structure=true)
	motor_pos_z = NaN64
	if new_data_structure == true
		reg = r"Z_[0-9]*.{0,1}[0-9]*mm"
		if occursin(reg, m.name)
			motor_pos_z = parse(Float64, match( reg, m.name ).match[3:end-2])
		end
	else
		reg = r"pos_[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*"
		if occursin(reg, m.name)
			pos_string = match( reg, m.name ).match
			value_reg = r"[0-9]*[.][0-9]*"
			z_string = matchall(value_reg, pos_string)[1]
			motor_pos_z = parse(Float64, z_string)
		end
	end
	return motor_pos_z
end

function get_motor_pos_phi_from_measurement_name(m::Measurement; new_data_structure=true)
	motor_pos_phi = NaN64
	if new_data_structure == true
		reg = r"Phi_[0-9]*.{0,1}[0-9]*deg"
		if occursin(reg, m.name)
			motor_pos_phi = parse(Float64, match( reg, m.name ).match[5:end-3])
		end
	else
		reg = r"pos_[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*-{0,1}_{0,1}[0-9]*.{0,1}[0-9]*"
		if occursin(reg, m.name)
			pos_string = match( reg, m.name ).match
			value_reg = r"[0-9]*[.][0-9]*"
			phi_string = matchall(value_reg, pos_string)[3]
			motor_pos_phi = parse(Float64, phi_string)
		end
	end
	return motor_pos_phi
end

function get_measurement_duration_from_measurement_names(m::Measurement; new_data_structure=true)
	measurement_duration = 0 # in sec
	if new_data_structure == true
		files = filter(x -> startswith(x, m.name), readdir(joinpath(m.path_to_raw_data, "../conv_data" ) ))
		reg1 = r"measuretime_[0-9]*sec"
		reg2 = r"t_[0-9]*sec"
		for file in files
			if occursin(reg1, m.name)
				measurement_duration += parse(Int, match( reg1, m.name ).match[13:end-3])
			elseif occursin(reg2, m.name)
				measurement_duration += parse(Int, match( reg2, m.name ).match[3:end-3])
			end
		end
		measurement_duration == 0 ? measurement_duration=-1 : nothing
	else
		reg = r"dur_[0-9]*s"
		if occursin(reg, m.data_set_name)
			measurement_duration = parse(Int, match( reg, m.data_set_name ).match[5:end-1])
		end
		if occursin(reg, m.name)
			measurement_duration = parse(Int, match( reg, m.name ).match[5:end-1])
		end
	end
	return measurement_duration
end

function get_temperature_from_measurement_name(m::Measurement; new_data_structure=true)
	temperature = NaN64
	if new_data_structure == true
		reg = r"T_[0-9]*.{0,1}[0-9]*K"
		if occursin(reg, m.name)
			temperature = parse(Float64, match( reg, m.name ).match[3:end-1])
		end
	end
	return temperature
end

function get_side_source_from_measurement_name(m::Measurement)
	ss = "unknown"
	reg = r"ssrc_[A-z]*[0-9]*"
	if occursin(reg, m.data_set_name)
		ss = match(reg, m.data_set_name).match[6:end]
	end
	if occursin(reg, m.name)
		ss = match(reg, m.data_set_name).match[6:end]
	end
	return ss
end
function get_top_source_from_measurement_name(m::Measurement)
	ts = "unknown"
	reg = r"tsrc_[A-z]*[0-9]*"
	if occursin(reg, m.data_set_name)
		ts = match(reg, m.data_set_name).match[6:end]
	end
	if occursin(reg, m.name)
		ts = match(reg, m.data_set_name).match[6:end]
	end
	return ts
end
function get_outside_source_from_measurement_name(m::Measurement)
	os = "unknown"
	reg = r"osrc_[A-z]*[0-9]*"
	if occursin(reg, m.data_set_name)
		os = match(reg, m.data_set_name).match[6:end]
	end
	if occursin(reg, m.name)
		os = match(reg, m.data_set_name).match[6:end]
	end
	return os
end

function Measurement(data_set_name::AbstractString, name::AbstractString; new_data_structure=true) #Give any(or the single one, if there is only one) complete filename of that measurement as 'name' argument.
	m = Measurement()
	# m.name = name
	# new_data_structure ? m.name = name[1:match(r"-\d{8}T\d{6}Z.*", name).offset-1] : m.name = name
	m.name = name
	m.data_set_name = data_set_name
	m.new_data_structure = new_data_structure
	if m.new_data_structure == true
		m.path_to_raw_data  = joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data" )
		m.path_to_conv_data = joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "conv_data" )
	else
		m.path_to_raw_data  = joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data",  m.name )
		m.path_to_conv_data = joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "conv_data", m.name )
	end
	m.results_path = results_path(m)
	m.datetime = get_datetime_from_measurement_name(m, new_data_structure=new_data_structure)
	m.motor_pos_r = get_motor_pos_r_from_measurement_name(m, new_data_structure=new_data_structure)
	m.motor_pos_z = get_motor_pos_z_from_measurement_name(m, new_data_structure=new_data_structure)
	m.motor_pos_phi = get_motor_pos_phi_from_measurement_name(m, new_data_structure=new_data_structure)
	m.side_source = get_side_source_from_measurement_name(m)
	m.top_source = get_top_source_from_measurement_name(m)
	m.outside_source = get_outside_source_from_measurement_name(m)
	m.measurement_duration = get_measurement_duration_from_measurement_names(m,new_data_structure=new_data_structure)
	m.temperature = get_temperature_from_measurement_name(m,new_data_structure=new_data_structure)
	m.pressure = get_pressure_from_measurement_name(m)
	return m
end

######### Creating Data Sets ################################

function scan_conv_data_folder_for_measurements(data_set_name::AbstractString; new_data_structure=true)
	if new_data_structure==true
		files = readdir(joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "conv_data") )
		files = filter( x -> endswith(x,".hdf5") , files)
		list_of_measurements = get_measurement_name_from_compressed_data_file_name.(files)
		list_of_measurements = union(list_of_measurements)
	else
		@info "not yet implemented"
	end
end
function data_set_from_conv_data(data_set_name::AbstractString; new_data_structure=true)
	data_set = Measurement[]
	if new_data_structure == true
		list_of_measurements = scan_conv_data_folder_for_measurements(data_set_name, new_data_structure=new_data_structure)
		for mn in list_of_measurements
			m = Measurement(data_set_name, mn)
			# println(in(m, data_set))
			push!(data_set, m) 
		end
	else
		@info "not yet implemented"
	end
	return data_set
end
function scan_raw_data_folder_for_measurements(data_set_name::AbstractString; new_data_structure=true)
	if new_data_structure==true
		files = readdir(joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data") )
		files = filter( x -> endswith(x,".log") , files)
		list_of_measurements = get_measurement_name_from_compressed_data_file_name.(files)
		list_of_measurements = union(list_of_measurements)
	else
		@info "not yet implemented"
	end
end
function data_set_from_raw_data(data_set_name::AbstractString; new_data_structure=true)
	data_set = Measurement[]
	if new_data_structure==true
		list_of_measurements = scan_raw_data_folder_for_measurements(data_set_name, new_data_structure=new_data_structure)
		for m in list_of_measurements
			println(m)
			push!(data_set, Measurement(data_set_name, m))
		end
		return data_set
	else
		@info "not yet implemented"
		return data_set
	end
end

function scan_data_set_for_measurement_names(data_set_name::AbstractString; new_data_structure=true)
	list_of_measurements = String[]
	if new_data_structure==true
		files = readdir(joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data") )
		files = filter( x -> (endswith(x,"raw.dat") || endswith(x,"dat.bz2")), files)
		list_of_measurements = get_measurement_name_from_compressed_data_file_name.(files)
		list_of_measurements = union(list_of_measurements)
	else
		measurement_subfolders = readdir(joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data") )
		list_of_measurements = (filter(x -> isdir(joinpath(USER_DATA_PATH, data_set_name[1:4], data_set_name, "raw_data", x)), measurement_subfolders ))
	end
	return list_of_measurements
end

function Data_set(data_set_name::AbstractString; new_data_structure=true)
	data_set = Measurement[]
	if new_data_structure==true
		list_of_measurements = scan_conv_data_folder_for_measurements(data_set_name)
		for m in list_of_measurements
			push!(data_set, Measurement(data_set_name, m))
		end
	else
		list_of_measurements = scan_data_set_for_measurement_names(data_set_name, new_data_structure=new_data_structure)
		for m in list_of_measurements
			push!(data_set, Measurement(data_set_name, m, new_data_structure=new_data_structure))
		end
	end
	return data_set
end
###############################################################
import Base.sort
function sort(data_set, keyword = "phi")
	sorted_dataset = deepcopy(data_set)
	keyword_array = []
	for m in data_set[:]
		keyword == "phi" ? append!(keyword_array, m.motor_pos_phi) : nothing
		keyword == "z"   ? append!(keyword_array, m.motor_pos_z)   : nothing
		keyword == "r"   ? append!(keyword_array, m.motor_pos_r)   : nothing
		keyword == "T"   ? append!(keyword_array, m.temperature)   : nothing
	end
	sorted_keyword_array = sort(keyword_array)
	icount = 1
	for entry in sorted_keyword_array[:]
		for m in data_set[:]
			compare_key = 0.
			keyword == "phi" ? compare_key = m.motor_pos_phi : nothing
			keyword == "z"   ? compare_key = m.motor_pos_z   : nothing
			keyword == "r"   ? compare_key = m.motor_pos_r   : nothing
			keyword == "T"   ? compare_key = m.temperature   : nothing
			if entry == compare_key
				sorted_dataset[icount]=m
			end
		end
		icount += 1
	end
	return sorted_dataset
end
###############################################################

function gather_absolute_paths_to_hdf5_input_files(m::Measurement)
	if m.new_data_structure==true
		# hdf5_data_path = joinpath(USER_DATA_PATH, "$(Dates.year(m.date))", m.data_set_name, "conv_data")
		# hdf5_data_path = m.path_to_raw_data[1:end-8] * "conv_data"
		hdf5_data_path = m.path_to_conv_data
		input_files = readdir(hdf5_data_path)
		input_files = filter(x -> endswith(x, ".hdf5"), input_files)
		input_files = filter(x -> startswith(x, m.name), input_files)
		input_files = [joinpath.(hdf5_data_path, file) for file in input_files]
	else
		hdf5_data_path = joinpath(USER_DATA_PATH, "$(Dates.year(m.date))", m.data_set_name, "conv_data", m.name)
		input_files = readdir(hdf5_data_path)
		input_files = filter(x -> endswith(x, ".hdf5"), input_files)
		input_files = [joinpath.(hdf5_data_path, file) for file in input_files]
	end
	return input_files
end

function get_path_to_data_summary_file(m::Measurement)
	cal_data_path = joinpath(USER_DATA_PATH, "$(Dates.year(m.date))", m.data_set_name, "cal_data")
	files = readdir(cal_data_path)
	files = filter(x -> endswith(x, ".hdf5"), files)
	files = filter(x -> startswith(x, m.name), files)
	files = [joinpath.(cal_data_path, file) for file in files]
	if length(files) == 1
		return files[1]
	elseif length(files) > 1
		error("Strange. More than one data summary file...")
	else
		error("No data summary file found for measurement $(m.name)")
	end
end

function results_path(m::Measurement)::AbstractString
	return joinpath(USER_OUTPUT_PATH, m.data_set_name)
end


function select(ms::Array{Measurement}, dimension::Symbol, value; approx::Bool=false)::Array{Measurement}
	new_ms = Measurement[]
	for m in ms
		if dimension == :r
			for m in ms
				pos = approx ? round(m.motor_pos_r, digits=0) : m.motor_pos_r
				value == pos ? push!(new_ms, m) : nothing
			end
		elseif dimension == :z
			for m in ms
				pos = approx ? round(m.motor_pos_z, digits=0) : m.motor_pos_z
				value == pos ? push!(new_ms, m) : nothing
			end
		elseif dimension == :phi
			for m in ms
				pos = approx ? round(m.motor_pos_phi, digits=0) : m.motor_pos_phi
				value == pos ? push!(new_ms, m) : nothing
			end
		else
			error("Wrong input for dimension. Given '$dimension', but only ':r',':z' or ':phi' are valid.")
		end
	end
	return new_ms
end
