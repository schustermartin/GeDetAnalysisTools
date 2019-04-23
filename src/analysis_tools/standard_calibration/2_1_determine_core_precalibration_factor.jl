

function determine_core_precalibration_factor_with_mpas(m::Measurement; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], overwrite=false,
    min_n_peaks = length(photon_lines), max_n_peaks = 4 * length(photon_lines), nbins=6000, peak_threshold=10., peak_sigma=4.0, 
    averWindow=3, deconIterations=3, create_plots=true, α=0.01, rtol=1e-3, inputfiles = missing)

    ismissing(inputfiles) ? inputfiles = gather_absolute_paths_to_hdf5_input_files(m) : nothing

    core::Int = 1

    mpas::Array{Float32, 1} = get_measured_pulse_amplitudes(m)[1, :]
    n_events = length(mpas)

    h_core = fit(Histogram, mpas, nbins=nbins, closed=:left)

    # npeaks = 0
    # fPositionX = Float64[]
    # h_peaks = Histogram(0:1:1, :left)
    # while npeaks < min_npeaks
    #     h_peaks, fPositionX = RadiationSpectra.peakfinder(h_core, σ = peak_sigma, threshold = peak_threshold, deconIterations = deconIterations, averWindow = averWindow)
    #     npeaks = length(fPositionX)
    #     peak_threshold = 0.9 * peak_threshold
    # end

    # c0_pre = RadiationSpectra.determine_calibration_constant_through_peak_ratios(fPositionX, photon_lines, α = α, rtol = rtol )
    c_precal, h_deconv, peakPositions, threshold = RadiationSpectra.determine_calibration_constant_through_peak_ratios(  h_core, photon_lines, 
                                                                                                        min_n_peaks=min_n_peaks, max_n_peaks = max_n_peaks, 
                                                                                                        threshold = peak_threshold, α=α, σ=peak_sigma, rtol=rtol)
   

    if create_plots
        p_core  = plot(h_core, st=:step, size=(1200,600), label="core mpa spectrum (uncalibrated)");
        p_peaks = plot(h_deconv, st=:step, size=(1200,600), label="");
        plot!(p_peaks, [peak_threshold], st=:hline, color=:red, label="threshold")
        # p_pcf  = plot(h_pcf, st=:step, size=(1200,600), label="");
        # plot!([c_precal], st=:vline, label="precalibration factor");
        h_cal = Histogram(-200:1:6000, :left)
        append!(h_cal, mpas .* c_precal)
        p_core_cal = plot(h_cal, st=:step, size=(1200,600), label="core energy spectrum (rough calibrated)");
        plot!(p_core_cal, photon_lines, st=:vline, label="photon lines used for calibration")
        # p = plot(p_core, p_peaks, p_pcf, p_core_cal, layout=(4,1), size=(1920,1080));
        p = plot(p_core, p_peaks, p_core_cal, layout=(3,1), size=(1920,1080));
        savefig(m, p, "2_1_core_precalibration", "core_energy_spectrum_precalibration", fmt=:png )
        p_core = 0; p_peaks = 0; p_pcf = 0; p=0;
    end

    return c_precal, h_core, h_deconv#, h_pcf
end
