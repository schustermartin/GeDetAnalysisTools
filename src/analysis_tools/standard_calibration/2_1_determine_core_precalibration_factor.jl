function determine_core_precalibration_factor_with_mpas(m::Measurement; photon_lines=[609.312, 911.204, 1120.287, 1460.830, 1764.494, 2614.533], overwrite=false,
    min_npeaks=10, nbins=6000, peak_threshold=10., peak_sigma=4., averWindow=3, deconIterations=3, alpha=0.005, create_plots=true, α = 0.01)

    inputfiles = gather_absolute_paths_to_hdf5_input_files(m)

    core::Int = 1

    mpas::Array{Float32, 1} = get_measured_pulse_amplitudes(m)[1, :]
    n_events = length(mpas)

    h_core = fit(Histogram, mpas, nbins=nbins, closed=:left)

    npeaks = 0
    fPositionX = Float64[]
    h_peaks = Histogram(0:1:1, :left)
    while npeaks < min_npeaks
        h_peaks, fPositionX = RadiationSpectra.peakfinder(h_core, σ = peak_sigma, threshold = peak_threshold, deconIterations = deconIterations, averWindow = averWindow)
        npeaks = length(fPositionX)
        peak_threshold = 0.9 * peak_threshold
    end

    c0_pre = RadiationSpectra.determine_calibration_constant_through_peak_ratios(fPositionX, photon_lines, α = α )

    if create_plots
        p_core  = plot(h_core, st=:step, size=(1200,600), label="core mpa spectrum (uncalibrated)");
        p_peaks = plot(h_peaks, st=:step, size=(1200,600), label="");
        plot!(p_peaks, [peak_threshold], st=:hline, color=:red, label="threshold")
        # p_pcf  = plot(h_pcf, st=:step, size=(1200,600), label="");
        # plot!([c0_pre], st=:vline, label="precalibration factor");
        h_cal = Histogram(-200:1:6000, :left)
        append!(h_cal, mpas .* c0_pre)
        p_core_cal = plot(h_cal, st=:step, size=(1200,600), label="core energy spectrum (rough calibrated)");
        plot!(p_core_cal, photon_lines, st=:vline, label="photon lines used for calibration")
        # p = plot(p_core, p_peaks, p_pcf, p_core_cal, layout=(4,1), size=(1920,1080));
        p = plot(p_core, p_peaks, p_core_cal, layout=(3,1), size=(1920,1080));
        savefig(m, p, "2_1_core_precalibration", "core_energy_spectrum_precalibration", fmt=:png )
        p_core = 0; p_peaks = 0; p_pcf = 0; p=0;
    end

    return c0_pre, h_core, h_peaks#, h_pcf
end
