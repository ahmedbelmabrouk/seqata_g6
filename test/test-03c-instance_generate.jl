lg1("Génération d'une instance... ")
inst = instance_build_mini10()
if lg2()
    # println("\n", to_s_alpx(inst))
    println("\n", to_s_alp(inst))
    print("Génération d'une instance... ")
end
lg1("=> ok\n")

lg1("Statistiques sur l'instance $(inst.name)... ")
@test inst.nb_planes == length(inst.planes)
@test inst.nb_planes == 10
@test inst.nb_kinds == 2
@test inst.name == "mini10"

# sommes des coûts précalculés pour l'avion p10
p10 = inst.planes[10]
cumul_cost_p10 = sum(p10.costs[p10.lb:p10.ub])
# @show p10.costs
# println("lb=", p10.lb)
# println("ub=", p10.ub)
# for (t,cost) in enumerate(p10.costs)
#     println("   p10.costs[$t]=$cost")
# end
# @show cumul_cost_p10
@test cumul_cost_p10 == 3.42639e6 # valeur entière donc ok

lg1("=> ok\n")
