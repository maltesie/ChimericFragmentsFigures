function plot_figure_4(assets_folder::String)

    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 700*resfactor))

    Label(fig[1,1, TopLeft()], "a", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)
    Label(fig[2,1, TopLeft()], "b", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    Label(fig[1,2, TopLeft()], "c", fontsize = 26,font = :bold,padding = (0, 5, 5, 0), halign = :right)

    df_old = DataFrame(CSV.File(joinpath(assets_folder, "hfq_lcd_old.csv")))
    df_new = DataFrame(CSV.File(joinpath(assets_folder, "hfq_lcd_new.csv")))

    tested_pairs = [
        ("VCA0946", "Spot42"),
        ("aphA", "Spot42"),
        ("cyaA", "Spot42"),
        ("VC0972", "Spot42"),
        ("Spot42", "VC0972"),
        ("VCA0181", "Spot42"),
        ("VC1043", "Spot42"),
        ("VC0966", "Spot42"),
        ("VC2030", "Spot42"),
        ("VCA0166", "Spot42"),
        ("VC1904", "Spot42"),
        ("VC2602", "Spot42"),
    ]
    old_pairs = [(row["name1"], row["name2"]) for row in eachrow(df_old)]

    tested_index = collect(findall(row->(row["name1"], row["name2"]) in tested_pairs, eachrow(df_new)))
    old_index = collect(findall(row->(row["name1"], row["name2"]) in old_pairs, eachrow(df_new)))
    new_index = vcat(tested_index, [i for i in 1:nrow(df_new) if !((i in old_index) | (i in tested_index))])

    df_total_old = DataFrame(CSV.File(joinpath(assets_folder, "table_s1_lcd.csv")))
    df_total_new = DataFrame(CSV.File(joinpath(assets_folder, "interactions_hfq_lcd_new.csv")))

    types = (["sRNA"], ["IGR"], ["CDS", "5UTR", "3UTR", "CDS_UTRS"])

    counter = zeros(Float64, 12)
    for row in eachrow(df_total_old)
        (ismissing(row.type1) || ismissing(row.type2)) && continue
        for (i, t) in enumerate(types)
            counter[i] += Int(row.type1 in t)
            counter[i] += Int(row.type2 in t)
        end
    end
    for row in eachrow(df_total_new)
        for (i, t) in enumerate(types)
            if row.nb_ints >= 20
                counter[i+6] += Int(row.type1 in t)
                counter[i+6] += Int(row.type2 in t)
                if row.fisher_fdr <= 0.05
                    counter[i+3] += Int(row.type1 in t)
                    counter[i+3] += Int(row.type2 in t)
                end
            end
            if row.bp_fdr <= 0.25
                counter[i+9] += Int(row.type1 in t)
                counter[i+9] += Int(row.type2 in t)
            end
        end
    end
    counter ./= 2
    ax1 = Axis(fig[1,1], ylabel="count", title="annotation types", xticks = (1:4, ["previous", "improved\nmapping", "no FDR,\nreads >= 20", "FDR < 0.25,\nreads >= 3"]))
    colors = Makie.wong_colors()
    groups = [1,2,3,1,2,3,1,2,3,1,2,3]

    barplot!(ax1, vcat(fill(1, 3), fill(2, 3), fill(3, 3), fill(4, 3)), counter, stack=groups, color=colors[groups])
    labels = ["sRNA", "IGR", "mRNA"]
    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    axislegend(ax1, elements, labels, position=:lt)
    hidexdecorations!(ax1, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)

    axb2 = Axis(fig[2,1], ylabel="-log10(complementarity FDR)", xlabel="log10(reads count)", title="Spot 42 interactions")

    scatter!(axb2, log10.(df_new[new_index, :nb_ints]), -1 .* log10.(df_new[new_index, :bp_fdr]),
        label="Fisher FDR > 0.05", markersize=9, color=RGBAf(colors[1].r, colors[1].g, colors[1].b, 0.8))


    scatter!(axb2, log10.(df_new[old_index, :nb_ints]), -1 .* log10.(df_new[old_index, :bp_fdr]),
        label="Fisher FDR <= 0.05", markersize=9, color=RGBAf(colors[2].r, colors[2].g, colors[2].b, 0.8)) #color=RGBAf(0.4, 0.25, 0.8, 0.7))

    p_big = decompose(Point2f, Circle(Point2f(0), 0.8))
    p_small = decompose(Point2f, Circle(Point2f(0), 0.4))
    scatter!(axb2, log10.(df_new[tested_index, :nb_ints]), -1 .* log10.(df_new[tested_index, :bp_fdr]), label="tested",
        marker=Polygon(p_big, [p_small]), markersize=7, color=RGBAf(0.8, 0.1, 0.0, 1.0))
    axislegend(axb2, position=:lt)
    scatter!(axb2, log10.(df_new[tested_index, :nb_ints]), -1 .* log10.(df_new[tested_index, :bp_fdr]),  color=colors[1], markersize=9)
    scatter!(axb2, log10.(df_new[tested_index, :nb_ints]), -1 .* log10.(df_new[tested_index, :bp_fdr]), label="tested",
        marker=Polygon(p_big, [p_small]), markersize=7, color=RGBAf(0.8, 0.1, 0.0, 1.0))

    df_plate = DataFrame(CSV.File(joinpath(assets_folder, "platereader.csv")))
    xlabels_latex = [rich(String(cl); font=:italic) for cl in df_plate.name]

    axb4 = Axis(fig[1:2,2], title="Spot 42 reporter assay", xlabel="relative fluorescence [AU]", yticks = (1:nrow(df_plate), xlabels_latex))

    groups = vcat(fill(1,nrow(df_plate)), fill(2,nrow(df_plate)))
    ctrl = Matrix(df_plate[:, [:ctrl1, :ctrl2, :ctrl3]])
    mean_ctrl = vec(mean(ctrl; dims=2))
    ctrl ./= mean_ctrl
    sd_ctrl = vec(std(ctrl; dims=2))
    spot = Matrix(df_plate[:, [:exp1, :exp2, :exp3]])
    spot ./= mean_ctrl
    mean_spot = vec(mean(spot; dims=2))
    sd_spot = vec(std(spot; dims=2))
    mean_ctrl = ones(nrow(df_plate))
    barplot!(axb4, 1:nrow(df_plate), mean_spot, direction=:x, color=colors[1])#, dodge=groups, color=colors[groups])
    #scatter!(axb4, collect(1:nrow(df_plate)) .- 0.35, vec(ctrl[:, 1]), color="black", markersize=5)
    #scatter!(axb4, collect(1:nrow(df_plate)) .- 0.2, vec(ctrl[:, 2]), color="black", markersize=5)
    #scatter!(axb4, collect(1:nrow(df_plate)) .- 0.05, vec(ctrl[:, 3]), color="black", markersize=5)
    #errorbars!(collect(1:nrow(df_plate)) .- 0.2, mean_ctrl, sd_ctrl, whiskerwidth=5)
    scatter!(axb4, vec(spot[:, 1]), collect(1:nrow(df_plate)) .- 0.2, color="black", markersize=5)
    scatter!(axb4, vec(spot[:, 2]), collect(1:nrow(df_plate)), color="black", markersize=5)
    scatter!(axb4, vec(spot[:, 3]), collect(1:nrow(df_plate)) .+ 0.2, color="black", markersize=5)
    errorbars!(mean_spot, collect(1:nrow(df_plate)), sd_spot, whiskerwidth=5, direction=:x)
    #labels = ["control", "Spot 42"]
    #elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    #axislegend(axb4, elements, labels, position=:rt)
    hideydecorations!(axb4, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = false, minorticks = false)

    save( "figure_4.svg", fig)
    save( "figure_4.png", fig, px_per_unit = 2)

end