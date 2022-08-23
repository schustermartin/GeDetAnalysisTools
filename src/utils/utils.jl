
export EChnVsChnHist
EChnVsChnHist(e::Array{<:Real, 2}, xchn, ychn; xedges = -1000:4:6000, yedges = -1000:4:6000) = fit(Histogram, (e[xchn, :], e[ychn, :]), (xedges, yedges))

EChnVsChnHist(m::Measurement, args...; kwargs...) = EChnVsChnHist(GAT.get_energies(m), args...; kwargs...)


function get_alpha_peak_position_and_sigma(m, to_int = true; channel = 1)
    h_e_alpha = fit(Histogram,get_energies(m)[channel,:], 5100:1:5600)
    fit_alpha = RadiationSpectra_beforeBAT.fit_single_peak_histogram_w_gaussian_refined(h_e_alpha)
    if to_int == true
        return round(Int,fit_alpha.fitted_parameters[3]), round(Int,fit_alpha.fitted_parameters[2]), h_e_alpha, fit_alpha
    else
        return fit_alpha.fitted_parameters[3], fit_alpha.fitted_parameters[2], h_e_alpha, fit_alpha
    end
end
function get_single_segment_energies(m,ss)
    energies = get_energies(m)
    ss_idcs = get_single_segment_indices(m)
    idcs_ss = findall(x->x==UInt8(ss+1),ss_idcs)
    return energies[:,idcs_ss]
end

function get_LEgamma_peak_position_and_sigma(m, to_int = true; channel = 1)
    h_e_leg = fit(Histogram,get_energies(m)[channel,:], 20:1:90)
    fit_leg = RadiationSpectra_beforeBAT.fit_single_peak_histogram_refined(h_e_leg, (45.0,63.0), fit_function =:Gauss_pol1)
    if to_int == true
        return round(Int,fit_leg.fitted_parameters[3]), round(Int,fit_leg.fitted_parameters[2]), h_e_leg, fit_leg
    else
        return fit_leg.fitted_parameters[3], fit_leg.fitted_parameters[2], h_e_leg, fit_leg
    end
end

# function get_LEgamma_peak_position_and_sigma(m, ss = missing, to_int = true; channel = 1)
#     h_e_leg = fit(Histogram,get_single_segment_energies(m,ss)[channel,:], 20:1:90)
#     fit_leg = RadiationSpectra_beforeBAT.fit_single_peak_histogram_w_gaussian_refined(h_e_leg)

#     h_e_leg_ss = fit(Histogram,get_single_segment_energies(m,ss)[ss+1,:], 20:1:90)
#     fit_leg_ss = RadiationSpectra_beforeBAT.fit_single_peak_histogram_w_gaussian_refined(h_e_leg_ss)

#     if to_int == true
#         return round(Int,fit_leg_ss.fitted_parameters[3]), round(Int,fit_leg.fitted_parameters[2]), h_e_leg, fit_leg,
#         round(Int,fit_leg_ss.fitted_parameters[3]), round(Int,fit_leg_ss.fitted_parameters[2]), h_e_leg_ss, fit_leg_ss 
#     else
#         return fit_leg.fitted_parameters[3], fit_leg.fitted_parameters[2], h_e_leg, fit_leg,
#         fit_leg_ss.fitted_parameters[3], fit_leg_ss.fitted_parameters[2], h_e_leg_ss, fit_leg_ss
#     end
# end

# function get_LEgamma_peak_position_and_sigma_passivation(m, ss, to_int = true; channel = 1)
#     h_e_leg = fit(Histogram,get_energies(m)[channel,:], 20:1:90)
#     fit_leg = RadiationSpectra_beforeBAT.fit_single_peak_histogram_w_gaussian_refined(h_e_leg)
#     if to_int == true
#         return round(Int,fit_leg.fitted_parameters[3]), round(Int,fit_leg.fitted_parameters[2]), h_e_lge, fit_leg
#     else
#         return fit_leg.fitted_parameters[3], fit_leg.fitted_parameters[2], h_e_leg, fit_leg
#     end
#     h_e_leg = fit(Histogram,get_energies(m)[channel,:], 20:1:90)
#     fit_leg = RadiationSpectra_beforeBAT.fit_single_peak_histogram_w_gaussian_refined(h_e_leg)
# end