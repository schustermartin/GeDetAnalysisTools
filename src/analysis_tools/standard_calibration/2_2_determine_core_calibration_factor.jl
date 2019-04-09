function determine_core_calibration_factor_with_mpas(m::Measurement, c_precal::Real; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], edges=0:1:3000, create_plots=true)
    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)
    core::Int = 1

    precal_energies::Array{Float32, 1} = c_precal .* get_measured_pulse_amplitudes(m)[1, :]
    n_events = length(precal_energies)

    c = 1.0
    h0 = fit(Histogram, precal_energies, edges, closed=:left)

    peak_fits = RadiationSpectra.FitFunction[ RadiationSpectra.FitFunction( gauss_plus_first_order_polynom  ) for ichn in 1:length(photon_lines) ]
    # peak_fits = GeDetSpectrumAnalyserTmp.Fit[]
    for (i, pl) in enumerate(photon_lines)
        line::Float64 = pl
        # fitrange = (line - 20 ):(line + 20 ) # +- 20 keV
        fitrange = (line - 20, line + 20)
        peak_fits[i].fitrange = fitrange
        first_bin = StatsBase.binindex(h0, fitrange[1])
        last_bin  = StatsBase.binindex(h0, fitrange[2])
        p0_sigma = 1.0  # 1keV
        p0_scale = (maximum(h0.weights[first_bin:last_bin]) - (h0.weights[first_bin] + h0.weights[last_bin]) / 2) * 2 * p0_sigma
        p0_mean = midpoints(h0.edges[1])[findmax(h0.weights)[2]] #line
        p0_bg_offset = (h0.weights[first_bin] + h0.weights[last_bin]) / 2
        p0_bg_slope = (h0.weights[last_bin] - h0.weights[first_bin]) / (fitrange[2] - fitrange[1])
        p0 = Float64[ p0_scale, p0_sigma, p0_mean, p0_bg_offset, p0_bg_slope ]
        peak_fits[i].initial_parameters = p0
        RadiationSpectra.lsqfit!(peak_fits[i], h0, estimate_uncertainties=false)
        # fr = GeDetSpectrumAnalyserTmp.fit(h0, fitrange, gauss_plus_first_order_polynom, p0, σ=1.0, estimate_uncertainties=false)
        # fr.uncertainties = GeDetSpectrumAnalyserTmp.estimate_uncertainties(fr, 1.0)
        # push!(peak_fits, fr)
    end

    fitted_peak_positions = [ fr.parameters[3] for fr in peak_fits ] ./ c_precal
    # fitted_peak_positions_err = [ fr.uncertainties[3] for fr in peak_fits ] ./ c_precal

    c_fit = RadiationSpectra.FitFunction( linear_function_fixed_offset_at_zero)
    c_fit.initial_parameters = [c_precal]
    RadiationSpectra.lsqfit!( c_fit, photon_lines, fitted_peak_positions, estimate_uncertainties=false) #, fitted_peak_positions_err )
    # c_fit = GeDetSpectrumAnalyserTmp.LSQFIT(photon_lines, fitted_peak_positions, fitted_peak_positions_err, linear_function_fixed_offset_at_zero, [c_precal] )
    # c_fit.uncertainties = GeDetSpectrumAnalyserTmp.estimate_uncertainties(c_fit, 1.0)
    # c = inv(c_fit.parameters[1])
    c_fit.uncertainties = zeros(Float64, length(c_fit.parameters))
    c_fit.uncertainties .= -1
    c = inv(c_fit.parameters[1])

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
