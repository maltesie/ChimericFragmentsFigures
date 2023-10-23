function coverage_plot!(ax::Axis, coverage::Coverage, interval::Interval, color)
    vals = values(coverage, interval)
    strand(interval) == '-' && reverse!(vals)
    band!(ax, collect(0:(rightposition(interval)-leftposition(interval))), zeros(rightposition(interval)-leftposition(interval)+1), vals, color=color)
end

function ligation_points_plot!(ax::Axis, interact::Interactions, n1::String, n2::String,  max_fdr::Float64)
    idx = findfirst(interact.nodes.name .== n)
    cds = interact.nodes.cds[idx] == 0 ? (interact.nodes.strand[idx] == '-' ? interact.nodes.right[idx] : interact.nodes.left[idx]) : interact.nodes.cds[idx]
    st = interact.nodes.strand[idx]

end

function plot_figure_6(assets_folder::String, interact::Interactions)
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 800*resfactor))

    colors = Makie.wong_colors()
    colors = [colors[1], RGBAf(95/255, 158/255, 160/255, 191/255)]

    gb = fig[1:3,1] = GridLayout()

    ax_graph = Axis(gb[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_graph = rotr90(load(joinpath(assets_folder, "vc0715:vc0719_targets.png")))
    hidedecorations!(ax_graph)
    image!(ax_graph, img_graph, aspect = DataAspect())

    ax_ligation_points = Axis(gb[2, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_ligation_points = rotr90(load(joinpath(assets_folder, "ligation_points_with_ascii_plot_rect.png")))
    hidedecorations!(ax_ligation_points)
    image!(ax_ligation_points, img_ligation_points, aspect = DataAspect())

    ga = fig[1:2,2] = GridLayout()

    ax_genomic_context = Axis(ga[1, 1:4], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img_genomic_context = rotr90(load(joinpath(assets_folder, "genomic_context.png")))
    hidedecorations!(ax_genomic_context)
    image!(ax_genomic_context, img_genomic_context, aspect = DataAspect())

    ax1_1_1 = Axis(ga[2, 1], ylabel="count", xticks=([1, 950], ["1", "950"]))#, title="basepairing predictions")

    #limits_low, limits_high = 2150, 2270
    limits_low, limits_high = 1, 1000

    #aggregation_plot!(ax1_1_1, interact, "VC0715:VC0719", 0.3)
    aggregation_plot!(ax1_1_1, interact, "VC0715", 0.3)
    hidexdecorations!(ax1_1_1, grid = false)
    hidespines!(ax1_1_1, :l, :t, :r)
    xlims!(ax1_1_1, limits_low, limits_high)

    #idx = findfirst(interact.nodes.name .== "VC0715:VC0719")
    idx = findfirst(interact.nodes.name .== "VC0715")
    ref = interact.nodes.ref[idx]
    st = interact.nodes.strand[idx]
    interval = Interval(ref, interact.nodes.left[idx], interact.nodes.right[idx], st)
    coverage_tex = Coverage(joinpath(assets_folder, "trimmed_tex_01_1_forward.bw"), joinpath(assets_folder, "trimmed_tex_01_1_reverse.bw"))
    coverage_term = Coverage(joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_forward.bw"), joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_reverse.bw"))

    ax1_2_1 = Axis(ga[3,1], ylabel="count", xticks=([1, 950], ["1", "950"]))#, title="dRNA-seq")
    coverage_plot!(ax1_2_1, coverage_tex, interval, colors[1])
    hidexdecorations!(ax1_2_1, grid = false)
    hidespines!(ax1_2_1, :l, :t, :r)
    xlims!(ax1_2_1, limits_low, limits_high)

    ax1_3_1 = Axis(ga[4,1], ylabel="count", xticks=([1, 950], ["1", "950"]))#, xlabel="coordinate in annotation", title="TERM-seq")
    coverage_plot!(ax1_3_1, coverage_term, interval, colors[1])
    hidespines!(ax1_3_1, :l, :t, :r)
    xlims!(ax1_3_1, limits_low, limits_high)

    linkxaxes!(ax1_1_1, ax1_2_1, ax1_3_1)

    ax1_1_2 = Axis(ga[2, 2:3], title="basepairing predictions", xticks=([2150, 2250], ["3100", "3200"]))

    limits_low, limits_high = 2135, 2270
    #limits_low, limits_high = 0, 2000

    aggregation_plot!(ax1_1_2, interact, "VC0715:VC0719", 0.3)
    #axislegend(ax1_1_2)
    hidexdecorations!(ax1_1_2, grid = false)
    hideydecorations!(ax1_1_2, grid = false)
    hidespines!(ax1_1_2, :l, :t, :r)
    xlims!(ax1_1_2, limits_low, limits_high)

    idx = findfirst(interact.nodes.name .== "VC0715:VC0719")
    ref = interact.nodes.ref[idx]
    st = interact.nodes.strand[idx]
    interval = Interval(ref, interact.nodes.left[idx], interact.nodes.right[idx], st)
    #coverage_tex = Coverage(joinpath(assets_folder, "trimmed_tex_01_1_forward.bw"), joinpath(assets_folder, "trimmed_tex_01_1_reverse.bw"))
    #coverage_term = Coverage(joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_forward.bw"), joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_reverse.bw"))

    ax1_2_2 = Axis(ga[3,2:3], title="dRNA-seq", xticks=([2150, 2250], ["3100", "3200"]))
    coverage_plot!(ax1_2_2, coverage_tex, interval, colors[1])
    hidexdecorations!(ax1_2_2, grid = false)
    hideydecorations!(ax1_2_2, grid = false)
    hidespines!(ax1_2_2, :l, :t, :r)
    xlims!(ax1_2_2, limits_low, limits_high)

    ax1_3_2 = Axis(ga[4,2:3], xlabel=rich("coordinate relative to ", rich("vc0715"; font=:italic)), title="TERM-seq", xticks=([2150, 2250], ["3100", "3200"]))
    coverage_plot!(ax1_3_2, coverage_term, interval, colors[1])
    hideydecorations!(ax1_3_2, grid = false)
    hidespines!(ax1_3_2, :l, :t, :r)
    xlims!(ax1_3_2, limits_low, limits_high)

    linkxaxes!(ax1_1_2, ax1_2_2, ax1_3_2)

    ax1_1_3 = Axis(ga[2, 4], xticks=([100, 900], ["3600", "4400"]))

    #limits_low, limits_high = 2150, 2270
    limits_low, limits_high = 1, 1000

    #aggregation_plot!(ax1_1_3, interact, "VC0715:VC0719", 0.3)
    aggregation_plot!(ax1_1_3, interact, "VC0719", 0.3)
    hidexdecorations!(ax1_1_3, grid = false)
    hideydecorations!(ax1_1_3, grid = false)
    hidespines!(ax1_1_3, :l, :t, :r)
    xlims!(ax1_1_3, limits_low, limits_high)

    elements = [PolyElement(polycolor = c) for c in (RGBAf(0.388, 0.431, 0.98, 0.6), RGBAf(0.937, 0.333, 0.231, 0.6))]
    #axislegend(ax5, elements, labels, position=:rt)
    Legend(ga[2, 4], elements, ["RNA1", "RNA2"])

    #idx = findfirst(interact.nodes.name .== "VC0715:VC0719")
    idx = findfirst(interact.nodes.name .== "VC0719")
    ref = interact.nodes.ref[idx]
    st = interact.nodes.strand[idx]
    interval = Interval(ref, interact.nodes.left[idx], interact.nodes.right[idx], st)
    #coverage_tex = Coverage(joinpath(assets_folder, "trimmed_tex_01_1_forward.bw"), joinpath(assets_folder, "trimmed_tex_01_1_reverse.bw"))
    #coverage_term = Coverage(joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_forward.bw"), joinpath(assets_folder, "trimmed_12898-mock-rep1_S33_R1_001_reverse.bw"))

    ax1_2_3 = Axis(ga[3,4], xticks=([100, 900], ["3600", "4400"]))#, ylabel="count", title="dRNA-seq")
    coverage_plot!(ax1_2_3, coverage_tex, interval, colors[1])
    hidexdecorations!(ax1_2_3, grid = false)
    hideydecorations!(ax1_2_3, grid = false)
    hidespines!(ax1_2_3, :l, :t, :r)
    xlims!(ax1_2_3, limits_low, limits_high)

    ax1_3_3 = Axis(ga[4,4], xticks=([100, 900], ["3600", "4400"]))#, ylabel="count", xlabel="position in IGR", title="TERM-seq")
    coverage_plot!(ax1_3_3, coverage_term, interval, colors[1])
    hideydecorations!(ax1_3_3, grid = false)
    hidespines!(ax1_3_3, :l, :t, :r)
    xlims!(ax1_3_3, limits_low, limits_high)

    linkxaxes!(ax1_1_3, ax1_2_3, ax1_3_3)

    linkyaxes!(ax1_1_1, ax1_1_2, ax1_1_3)

    linkyaxes!(ax1_2_1, ax1_2_2, ax1_2_3)

    linkyaxes!(ax1_3_1, ax1_3_2, ax1_3_3)

    gc = fig[1:3,3] = GridLayout()

    img4 = rotr90(load(joinpath(assets_folder, "nb.png")))
    ax4 = Axis(gc[1:6, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())


    df_plate = DataFrame(CSV.File(joinpath(assets_folder, "platereader.csv")))

    ax7 = Axis(fig[3,2], ylabel="AphA-FLAG levels [AU]", xlabel="OD at 600nm",
        xticks = (1:nrow(df_plate), ["$n" for n in df_plate.name]))
    groups = vcat(fill(1,nrow(df_plate)), fill(2,nrow(df_plate)))

    ctrl = Matrix(df_plate[:, [:ctrl1, :ctrl2, :ctrl3]])
    std_one = std(ctrl[1,:]) / mean(ctrl[1,:])
    spot = Matrix(df_plate[:, [:exp1, :exp2, :exp3]])
    sd_spot_one = std(spot[1,:] ./ mean(ctrl[1,:]))

    factor = mean(spot[1,:]) ./ mean(ctrl[1,:])
    sd_factor = factor * sqrt(sd_spot_one^2 + std_one^2)
    ctrl = (ctrl' ./ ctrl[1,:])'
    mean_ctrl = vec(mean(ctrl; dims=2))

    sd_ctrl = vec(std(ctrl; dims=2))
    sd_ctrl[1] = std_one

    spot = (spot' ./ spot[1,:])'
    mean_spot_before = vec(mean(spot; dims=2))
    sd_spot_before = vec(std(spot; dims=2))

    spot .*= factor
    mean_spot = vec(mean(spot; dims=2))
    sd_spot = mean_spot .* sqrt.((sd_spot_before ./ mean_spot_before).^2 .+ (sd_factor/factor)^2)
    sd_spot[1] = sd_spot_one

    barplot!(ax7, vcat(1:nrow(df_plate), 1:nrow(df_plate)), vcat(mean_ctrl, mean_spot), dodge=groups, color=colors[groups])
    scatter!(ax7, collect(1:nrow(df_plate)) .- 0.35, vec(ctrl[:, 1]), color="black", markersize=5)
    scatter!(ax7, collect(1:nrow(df_plate)) .- 0.2, vec(ctrl[:, 2]), color="black", markersize=5)
    scatter!(ax7, collect(1:nrow(df_plate)) .- 0.05, vec(ctrl[:, 3]), color="black", markersize=5)
    errorbars!(ax7, collect(2:nrow(df_plate)) .- 0.2, mean_ctrl[2:end], sd_ctrl[2:end], whiskerwidth=5)
    scatter!(ax7, collect(1:nrow(df_plate)) .+ 0.35, vec(spot[:, 1]), color="black", markersize=5)
    scatter!(ax7, collect(1:nrow(df_plate)) .+ 0.2, vec(spot[:, 2]), color="black", markersize=5)
    scatter!(ax7, collect(1:nrow(df_plate)) .+ 0.05, vec(spot[:, 3]), color="black", markersize=5)
    errorbars!(ax7, collect(2:nrow(df_plate)) .+ 0.2, mean_spot[2:end], sd_spot[2:end], whiskerwidth=5)
    labels = ["control", "NetX"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    hidexdecorations!(ax7, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)

    axislegend(ax7, elements, labels, position=:rt)

    Label(gb[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[2,1, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    Label(ga[1,1, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(fig[3,2, TopLeft()], "e", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    Label(gc[1,1, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    #Label(gc[7,1, TopLeft()], "F", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    save("figure_6.svg", fig)
    save("figure_6.png", fig, px_per_unit = 2)
end