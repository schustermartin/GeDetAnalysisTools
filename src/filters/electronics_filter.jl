function the_filter(x, p::NamedTuple; 
        Samplingtime::Float64=4e-9)
    # Filter 1: Basic RC filter:
    #-------------------------------
    K  = 2*p.Rise_Time/Samplingtime
    Fc = 1/(2*pi*p.Rise_Time)


    myfilter = PolynomialRatio([1,1], [1+K, 1-K])
    filtered = filt(myfilter, x)

    if p.Window != 0
    # Filter 2: Basic Hamming Window:
    #-------------------------------
        responsetype = Lowpass(Fc; fs=1/Samplingtime)
        designmethod = FIRWindow(hamming(round(Int, p.Window)))
        filtered = filt(digitalfilter(responsetype, designmethod), filtered)
    end

    if p.Butter != 0
    # Filter 3: Basic Butterworth:
    #-------------------------------
        responsetype = Bandpass(1, Fc; fs=1/Samplingtime)
        designmethod = Butterworth(round(Int, p.Butter))
        filtered = filt(digitalfilter(responsetype, designmethod), filtered)
    end

    # x^-2 folding for a longer/smoother rise (especially for the segments)
    f(x) = x .^(-2)
    vec = f(0.1 : p.fct : 3)
    vec = vec ./ sum(vec)

    return conv(filtered, vec)[1:length(x)]
end



function electronics_filter(x, p; Sampling_Time::Float64=4e-9)
    
    if length(x[1]) == 1
        
        if typeof(p) == NamedTuple{(:Window, :Butter, :fct, :Rise_Time),Tuple{Int64,Int64,Float64,Float64}}
            return [the_filter(x,p, Samplingtime=Sampling_Time)]
            
        elseif typeof(p[1]) == NamedTuple{(:Window, :Butter, :fct, :Rise_Time),Tuple{Int64,Int64,Float64,Float64}}
            Filtered_Pulses = []
            i = 1
            while i <= length(p)
                push!(Filtered_Pulses, the_filter(x, p[i], Samplingtime=Sampling_Time))
                i+= 1
            end
            return Filtered_Pulses
        end
                
        
    elseif length(x) != length(p)
        if typeof(p) != NamedTuple{(:Window, :Butter, :fct, :Rise_Time),Tuple{Int64,Int64,Float64,Float64}}
            if length(x) == 1
                Filtered_Pulses = []
                i = 1
                while i <= length(p)
                    push!(Filtered_Pulses, the_filter(x[1], p[i], Samplingtime=Sampling_Time))
                    i+= 1
                end
                return Filtered_Pulses
            else
                
                error("Check the dimensions! length(pulses) = "*string(length(x))*", length(p) = "*string(length(p)))
            end
            
               
        else
            Filtered_Pulses = []
            i = 1
            while i <= length(x)
                push!(Filtered_Pulses, the_filter(x[i], p, Samplingtime=Sampling_Time))
                i+= 1
            end
            return Filtered_Pulses
        end
        
        
    else
        Filtered_Pulses = []
        i = 1
        while i <= length(x)
            push!(Filtered_Pulses, the_filter(x[i], p[i], Samplingtime=Sampling_Time))
            i+= 1
        end
        return Filtered_Pulses
            
    end
        
end
