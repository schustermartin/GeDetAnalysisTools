mutable struct DAQ
  name::AbstractString
  n_channels::Int
  n_samples::Int
  sampling_rate::Float64  # in secs
  baseline_length::Int
  decay_window_length::Int
  function DAQ()
    new("noname", -1, -1, NaN, -1, -1 )
  end
end

function PIXI()
    pixi = DAQ()
    pixi.name = "PIXI"
    pixi.n_channels = 20
    pixi.n_samples = 1023
    pixi.sampling_rate = 75e6
    pixi.baseline_length = 338  # 2016
    pixi.decay_window_length = 500  #2016
    return pixi
end

function STRUCK()
    struck = DAQ()
    struck.name = "STRUCK"
    struck.n_channels = 16
    struck.n_samples = 5000
    struck.sampling_rate = 250e6
    struck.baseline_length = 2000
    struck.decay_window_length = 5000-2750
    return struck
end

function STRUCK_CZT_SETTINGS()
    struck = DAQ()
    struck.name = "STRUCK"
    struck.n_channels = 16
    struck.n_samples = 600
    struck.sampling_rate = 250e6
    struck.baseline_length = 150
    struck.decay_window_length = 120
    return struck
end

function TWO_STRUCKS()
    struck = DAQ()
    struck.name = "STRUCK"
    struck.n_channels = 20
    struck.n_samples = 5000
    struck.sampling_rate = 250e6
    struck.baseline_length = 1200
    struck.decay_window_length = 5000-2500
    return struck
end

export get_time_array_of_pulses
function get_time_array_of_pulses(daq::DAQ)::Array{Float64,1}
    return collect( range(0, step= 1 / daq.sampling_rate, length=daq.n_samples)) # in s
end
function get_time_array_of_pulses(::Type{T}, daq::DAQ)::Array{T,1} where T
    return collect(T, range(0, step= 1 / daq.sampling_rate, length=daq.n_samples)) # in s
end
