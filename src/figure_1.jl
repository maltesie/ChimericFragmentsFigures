function plot_figure_1(assets_dir::String)
    resfactor = 1.
    fig = Figure(resolution=(1200*resfactor, 800*resfactor))

    ga = fig[1:2,1] = GridLayout()

    ax1 = Axis(ga[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img1 = rotr90(load(joinpath(assets_dir, "1A.png")))
    hidedecorations!(ax1)
    image!(ax1, img1, aspect = DataAspect())

    ax2 = Axis(ga[1, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img2 = rotr90(load(joinpath(assets_dir, "1B.png")))
    hidedecorations!(ax2)
    image!(ax2, img2, aspect = DataAspect())

    ax3 = Axis(ga[1, 3], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img3 = rotr90(load(joinpath(assets_dir, "1C.png")))
    hidedecorations!(ax3)
    image!(ax3, img3, aspect = DataAspect())

    gb = fig[3:5,1] = GridLayout()

    ax4 = Axis(gb[1, 1], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img4 = rotr90(load(joinpath(assets_dir, "1D.png")))
    hidedecorations!(ax4)
    image!(ax4, img4, aspect = DataAspect())

    ax5 = Axis(gb[1, 2], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img5 = rotr90(load(joinpath(assets_dir, "1E.png")))
    hidedecorations!(ax5)
    image!(ax5, img5, aspect = DataAspect())

    ax6 = Axis(gb[1, 3], leftspinevisible = false, rightspinevisible = false, bottomspinevisible = false, topspinevisible = false, aspect = DataAspect())
    img6 = rotr90(load(joinpath(assets_dir, "1F.png")))
    hidedecorations!(ax6)
    image!(ax6, img6, aspect = DataAspect())

    Label(ga[1,1, TopLeft()], "a", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,2, TopLeft()], "b", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(ga[1,3, TopLeft()], "c", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,1, TopLeft()], "d", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,2, TopLeft()], "e", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)
    Label(gb[1,3, TopLeft()], "f", fontsize = 26, font = :bold, padding = (0, 5, 5, 0), halign = :right)

    #titlelayout = GridLayout(fig[0, 1], halign = :left, tellwidth = false)
    #Label(titlelayout[1, 1], "Fig. 1", halign = :left, fontsize=30)

    save("figure_1.pdf", fig)
    save("figure_1.png", fig, px_per_unit = 2)
end
