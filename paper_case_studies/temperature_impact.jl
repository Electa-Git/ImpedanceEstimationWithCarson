##################################################
### This script analyzes the impact of 
### temperature on steady-state self-resistances
### with Carson's equations, and generates the
### picture featured in our journal paper
###################################################
import StatsPlots as _SP
using LaTeXStrings

T = 40:1:90
T65 = 65

# materials
materials = Dict(
   "Cu"      => Dict("ρ" => 17.77e-9, 
                     "α" => 0.00381),
   "Al-1350" => Dict("ρ" => 2.83e-8, 
                     "α" => 0.00403)
)


# cross-sections
A_cu = [141.03, 18.09, 28.27]
A_al = [265.90, 216.42, 18.09, 277.59, 18.86, 15.27]
gmr_cu = exp(-1/4)*sqrt.(A_cu./pi)
gmr_al = exp(-1/4)*sqrt.(A_al./pi)

# y-axis
R_cu(T::Int, mat::Dict, A::Float64) = 0.049348 + mat["Cu"]["ρ"]/(A/1E6) * (1+(mat["Cu"]["α"]*(T-20)))
R_al(T::Int, mat::Dict, A::Float64) = 0.049348 + mat["Al-1350"]["ρ"]/(A/1E6) * (1+(mat["Al-1350"]["α"]*(T-20)))
# X(A::Float64) = 0.062832*(log(1/( exp(-1/4)*sqrt(A/pi)*3.28084*1E-3 ))+8.0252)

# Rdiff(R)    = (R65-R)/R65*100  
# Zdiff(X, R) = (Z65-R)/R65*100

Rdiff_cu(A,T) = (R_cu(T, materials, A)-R_cu(65, materials, A))/R_cu(65, materials, A)*100
Rdiff_al(A,T) = (R_al(T, materials, A)-R_al(65, materials, A))/R_al(65, materials, A)*100

lower_cu = minimum.([[Rdiff_cu(A,t) for A in A_cu] for t in T])
upper_cu = maximum.([[Rdiff_cu(A,t) for A in A_cu] for t in T])

lower_al = minimum.([[Rdiff_al(A,t) for A in A_al] for t in T])
upper_al = maximum.([[Rdiff_al(A,t) for A in A_al] for t in T])

p = _SP.plot(xticks = (40:5:90, [L"%$i" for i in 40:5:90]), yticks = (-0.3:0.1:0.3, [L"%$i" for i in -0.3:0.1:0.3]), xlabel = L"\textrm{Temperature} \, \, [°\textrm{C}]", 
            ylabel = L"(\textrm{R}^{65}-\textrm{R})/\textrm{R}^{65} \times 100 \, \, [\%]", legend = :topleft,
            xtickfontsize=12, ytickfontsize=12, xlabelfontsize=15, labelfontsize=15, legendfontsize=13)

_SP.plot!(T, lower_cu, fillrange = upper_cu, color = :black, label = L"\textrm{Cu}", fillalpha = 0.2)
_SP.plot!(T, lower_al, fillrange = upper_al, color = :grey, label = L"\textrm{Al-1350}", fillalpha = 0.2)
_SP.plot!(T, upper_cu, color = :black, label = "", fillalpha = 0.3)
_SP.plot!(T, upper_al, color = :grey, label = "", fillalpha = 0.2)

font = text("").font

font.rotation = 1
annotate!(50, 0.01, text(L"Al, 277.59 mm^2", font, :grey, 14))

font.rotation = 35
annotate!(50, -0.25, text(L"Al, 15.27 mm^2", font, :grey, 14))

font.rotation = 20
annotate!(80, 0.14, text(L"Cu, 18.09 mm^2", font, :black, 14))

font.rotation = 2
annotate!(80, -0.015, text(L"Cu, 141.03 mm^2", font, :black, 14))

_SP.savefig(p, "temperature_impact.png")
_SP.savefig(p, "temperature_impact.pdf")
# for A in A_al
#     _SP.plot!(T, [Rdiff_al(A,t) for t in T], linestyle = :dot, label = "$A")
# end