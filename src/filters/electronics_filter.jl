#=
Manual:
-----------------------------------------
Required input:
x                     -> your simulated pulse

Optional input:
Segment::String       -> Core (default), Seg1, Seg2, Seg3, Seg4
- This takes the best parameters found so far for the chosen response function
p::NamedTuple         -> If you want to use your own parameters.
Samplingtime::Float64 -> Default 4e-9

=#

function electronics_filter(x ;
        Segment::String="Core",
        p::NamedTuple=(Rise_Time=9.42e-9, Window=6, butter=1, fct=0.2), 
        Samplingtime::Float64=4e-9)
    
    # if p is set to Default, the function takes the input of 'Segment'.
    # And if 'Segment' is not entered, it will use the Core parameters
    
    if p == (Rise_Time=9.42e-9, Window=6, butter=1, fct=0.2)
        if Segment == "Core"
            p = (Rise_Time=9.42e-9, Window = 6, butter = 1, fct = 0.2);
        elseif Segment == "Seg1"
            p = (Rise_Time=9.42e-9, Window = 60, butter = 1, fct = 0.015);
        elseif Segment == "Seg2"
            p = (Rise_Time=9.42e-9, Window = 60, butter = 1, fct = 0.0145);
        elseif Segment == "Seg3"
            p = (Rise_Time=10.42e-9, Window = 60, butter = 1, fct = 0.011);
        elseif Segment == "Seg4"
            p = (Rise_Time=9.42e-9, Window = 80, butter = 1, fct = 0.006);
        else # Default is "Core"
            p = (Rise_Time=9.42e-9, Window = 6, butter = 1, fct = 0.2);
        end
    end
    
    
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
    
    if p.butter != 0
    # Filter 3: Basic Butterworth:
    #-------------------------------
        responsetype = Bandpass(1, Fc; fs=1/Samplingtime)
        designmethod = Butterworth(round(Int, p.butter))
        filtered = filt(digitalfilter(responsetype, designmethod), filtered)
    end
    
    # x^-2 folding for a longer/smoother rise (especially for the segments)
    f(x) = x .^(-2)
    vec = f(0.1 : p.fct : 3)
    vec = vec ./ sum(vec)

    return conv(filtered, vec)[1:length(x)]
end