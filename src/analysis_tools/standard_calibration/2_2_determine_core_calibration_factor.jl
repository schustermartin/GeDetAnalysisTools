function determine_core_calibration_factor_with_mpas(d::Dataset, c_precal::Real; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], edges=0:1:3000, create_plots=true)
    inputfiles::Vector{AbstractString}=[]
    for m in d
        push!(inputfiles,gather_absolute_paths_to_hdf5_input_files(m))
    end
    return determine_core_calibration_factor_with_mpas(d[1], c_precal, photon_lines=photon_lines, edges=edges, create_plots = create_plots, inputfiles=inputfiles)
end

function determine_core_calibration_factor_with_mpas(m::Measurement, c_precal::Real; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], edges=0:1:3000, create_plots=true, inputfiles::Union{Vector{AbstractString},Missing} = missing)

    ismissing(inputfiles) ? inputfiles = gather_absolute_paths_to_hdf5_input_files(m) : nothing
    core::Int = 1

    precal_energies::Array{Float32, 1} = c_precal .* get_measured_pulse_amplitudes(m)[1, :]
    n_events = length(precal_energies)

    T = Float64

    c::T = 1.0
    h0 = fit(Histogram, precal_energies, edges, closed=:left)

    peak_fits = RadiationSpectra.FitFunction[ RadiationSpectra.FitFunction{T}( gauss_plus_first_order_polynom, 1, 5 ) for ichn in 1:length(photon_lines) ]
    # peak_fits = GeDetSpectrumAnalyserTmp.Fit[]
    for (i, pl) in enumerate(photon_lines)
        line::T = pl
        # fitrange = (line - 20 ):(line + 20 ) # +- 20 keV
        fitrange = (line - 20, line + 20)
        set_fitranges!(peak_fits[i], (fitrange, ))
        first_bin = StatsBase.binindex(h0, fitrange[1])
        last_bin  = StatsBase.binindex(h0, fitrange[2])
        p0_sigma = 1.0  # 1keV
        p0_sigma = line / 900.
        p0_scale = 0.5 * (maximum(h0.weights[first_bin:last_bin]) - (h0.weights[first_bin] + h0.weights[last_bin]) / 2) * 2 * p0_sigma
        p0_mean = midpoints(h0.edges[1])[findmax(h0.weights)[2]] #
        idx_max = findmax(h0.weights[StatsBase.binindex(h0, line - 20):StatsBase.binindex(h0, line + 20)])[2] + StatsBase.binindex(h0, line - 20)
        # p0_mean = line
        p0_mean = h0.edges[1][idx_max]
        p0_bg_offset = (h0.weights[first_bin] + h0.weights[last_bin]) / 2
        p0_bg_slope = (h0.weights[last_bin] - h0.weights[first_bin]) / (fitrange[2] - fitrange[1])
        p0 = T[ p0_scale, p0_sigma, p0_mean, p0_bg_offset, p0_bg_slope ]
        set_initial_parameters!(peak_fits[i], p0)
        RadiationSpectra.lsqfit!(peak_fits[i], h0)
    end

    fitted_peak_positions = [ fr.fitted_parameters[3] for fr in peak_fits ] ./ c_precal
    # fitted_peak_positions_err = [ fr.uncertainties[3] for fr in peak_fits ] ./ c_precal

    c_fit = RadiationSpectra.FitFunction{T}( linear_function_fixed_offset_at_zero, 1, 1)
    set_initial_parameters!(c_fit, [c_precal])
    set_fitranges!(c_fit, ((minimum(photon_lines), maximum(photon_lines)),))
    RadiationSpectra.lsqfit!( c_fit, photon_lines, fitted_peak_positions) #, fitted_peak_positions_err )
    c = inv(c_fit.fitted_parameters[1])

    if create_plots
        peak_fit_plots = []
        for (ipl, pl) in enumerate(photon_lines)
            line::Float64 = pl
            fitrange = (line - 20 ):(line + 20 ) # +- 20 keV
            first_bin = StatsBase.binindex(h0, first(fitrange))
            last_bin  = StatsBase.binindex(h0, last(fitrange))
            pfp = plot( h0.edges[1][first_bin]:step(h0.edges[1]):h0.edges[1][last_bin], h0.weights[first_bin:last_bin], st=:step, label="Data (precalibrated)" )
            plot!(peak_fits[ipl], label="LSQ Fit")
            plot!(peak_fits[ipl], label="LSQ Fit - Init", use_initial_parameters = true)
            plot!([pl], st=:vline, color=:green, label="Photon line")
            push!(peak_fit_plots, pfp)
        end
        p_fits = plot( peak_fit_plots..., size=(1920,1080));

        p_factor_fit = plot(photon_lines, fitted_peak_positions, st=:scatter, legend=false, xlabel="E / keV", ylabel="MPA / precalibrated")
        # p_factor_fit = plot(photon_lines, fitted_peak_positions, yerr=fitted_peak_positions_err, st=:scatter, legend=false, xlabel="E / keV", ylabel="MPA / precalibrated")
        plot!(c_fit)
        l = @layout [a{0.75h}
                            b]

        p = plot(p_fits, p_factor_fit, layout=l, size=(1920,1080))
        savefig(m, p, "2_2_core_calibration", "core_calibration_peak_fits", fmt=:png )
        p = 0; peak_fit_plots=0; p_fits=0; p_factor_fit=0;
    end

    return c, h0, peak_fits, c_fit
end
