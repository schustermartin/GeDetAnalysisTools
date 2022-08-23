# import Plots.xlims
# import Plots.ylims
# function rel_coord(x, y)
#     ax = x*(xlims()[2] - xlims()[1]) + xlims()[1]
#     ay = y*(ylims()[2] - ylims()[1]) + ylims()[1]
#     return ax, ay
# end
import Plots.savefig
function savefig(m::Measurement, plt, subdir, name; fmt=:pdf)
    output_path = joinpath(m.results_path, subdir)
    !ispath(output_path) ? mkpath(output_path) : nothing
    if fmt == :pdf
        pdf(plt, joinpath(m.results_path, subdir, "$(name)_$(m.name)"))
    elseif fmt == :png
        # png(plt, joinpath(m.results_path, subdir, "$(name)_$(m.name)"))
        savefig(plt, joinpath(m.results_path, subdir, "$(name)_$(m.name)"*".png"))
    elseif fmt == :eps
        eps(plt, joinpath(m.results_path, subdir, "$(name)_$(m.name)"))
    else
        pdf(plt, joinpath(m.results_path, subdir, "$(name)_$(m.name)"))
    end
    nothing
end
