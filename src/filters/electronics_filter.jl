"""
    the_filter(x, p::NamedTuple[, samplingtime::Float64=4e-9])
...
# Arguments
- `x`: pulse in array.
- `p::NamedTuple`: one set of parameters for the filter.
- `samplingtime::Float64=4e-9`: sampling time of the pulse.
...
"""
function the_filter(x, p::NamedTuple;
        samplingtime::Float64=4e-9)
    # Filter 1: Basic RC filter:
    #-------------------------------
    K  = 2*p.risetime/samplingtime
    Fc = 1/(2*pi*p.risetime)


    myfilter = PolynomialRatio([1,1], [1+K, 1-K])
    filtered = filt(myfilter, x)

    if p.window != 0
    # Filter 2: Basic Hamming Window:
    #-------------------------------
        responsetype = Lowpass(Fc; fs=1/samplingtime)
        designmethod = FIRWindow(hamming(round(Int, p.window)))
        filtered = filt(digitalfilter(responsetype, designmethod), filtered)
    end

    if p.butter != 0
    # Filter 3: Basic Butterworth:
    #-------------------------------
        responsetype = Bandpass(1, Fc; fs=1/samplingtime)
        designmethod = Butterworth(round(Int, p.butter))
        filtered = filt(digitalfilter(responsetype, designmethod), filtered)
    end

    # x^-2 folding for a longer/smoother rise (especially for the segments)
    f(x) = x .^(-2)
    vec = f(0.1 : p.fct : 3)
    vec = vec ./ sum(vec)

    return conv(filtered, vec)[1:length(x)]
end


"""
    electronics_filter(x, p[, samplingtime::Float64=4e-9])
...
# Arguments
- `x`: pulse in array or array of pulses.
- `p`: set of parameters for the filter or array of parameter sets.
- `samplingtime::Float64=4e-9`: sampling time of the pulse.
...
"""
function electronics_filter(x, p; sampling_time::Float64=4e-9)

    if length(x[1]) == 1

        if typeof(p) == NamedTuple{(:window, :butter, :fct, :risetime),Tuple{Int64,Int64,Float64,Float64}}
            return [the_filter(x,p, samplingtime=sampling_time)]

        elseif typeof(p[1]) == NamedTuple{(:window, :butter, :fct, :risetime),Tuple{Int64,Int64,Float64,Float64}}
            filtered_pulses = []
            i = 1
            while i <= length(p)
                push!(filtered_pulses, the_filter(x, p[i], samplingtime=sampling_time))
                i+= 1
            end
            return filtered_pulses
        end


    elseif length(x) != length(p)
        if typeof(p) != NamedTuple{(:window, :butter, :fct, :risetime),Tuple{Int64,Int64,Float64,Float64}}
            if length(x) == 1
                filtered_pulses = []
                i = 1
                while i <= length(p)
                    push!(filtered_pulses, the_filter(x[1], p[i], samplingtime=sampling_time))
                    i+= 1
                end
                return filtered_pulses
            else

                error("Check the dimensions! length(pulses) = "*string(length(x))*", length(p) = "*string(length(p)))
            end


        else
            filtered_pulses = []
            i = 1
            while i <= length(x)
                push!(filtered_pulses, the_filter(x[i], p, samplingtime=sampling_time))
                i+= 1
            end
            return filtered_pulses
        end


    else
        filtered_pulses = []
        i = 1
        while i <= length(x)
            push!(filtered_pulses, the_filter(x[i], p[i], samplingtime=sampling_time))
            i+= 1
        end
        return filtered_pulses

    end

end
