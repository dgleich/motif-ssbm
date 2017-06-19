##
using Plots
using PlotRecipes
using MatrixNetworks
using MatrixDepot
##
# matrixdepot("Pothen/commanche_dual", :get)
##
A = matrixdepot("Pothen/commanche_dual", :r)
xy = reshape(
      matrixdepot("Pothen/commanche_dual", :r, meta=true)["commanche_dual_coord"],
      size(A,1),3)



##
#A,xy = lollipop_graph(5, Val{true})
pyplot(size=(1200,600))
graphplot(findnz(A)[1:2]..., x =xy[:,1], y=xy[:,2], z=xy[:,3],
  dim=3,
  curves=false,
  markersize=1,
  linewidth=0.75, markerstrokewidth=1, markercolor=:white,
  markerstrokecolor=:blue,
  linealpha=0.5,
  linecolor=:grey, foreground_color_axis=:red, grid=false, markeralpha=0.5,
  ticks=nothing)
#Plots.PyPlotBackend.gca()
gui()
Plots.PyPlot.gca()[:set_axis_off]()
Plots.PyPlot.savefig("commanche-figure.png", dpi=300)
#savefig("commanche-figure.pdf")
#savefig("commanche-figure.png")
