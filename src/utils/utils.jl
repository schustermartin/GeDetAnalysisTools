
export EChnVsChnHist
EChnVsChnHist(e::Array{<:Real, 2}, xchn, ychn; xedges = -1000:4:6000, yedges = -1000:4:6000) = fit(Histogram, (e[xchn, :], e[ychn, :]), (xedges, yedges))

EChnVsChnHist(m::Measurement, args...; kwargs...) = EChnVsChnHist(GAT.get_energies(m), args...; kwargs...)