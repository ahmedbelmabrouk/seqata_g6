
# EXEMPLE DE LIGNE DE COMMANDE TESTÉE :
# ./bin/run.jl val -t lp -i data/01.alp -s sols/alp_01_p10_k3=700.0.sol

# ===========
# Construction de l'instance mini
inst = instance_build_mini10()
@test inst.nb_planes == 10

# Création de la solution vide
sol = Solution(inst, algo = :earliest) # sera refait pour chaque test
# @show join(sol.planes, ",")

lg1("1. La solution initiale de alp10 doit être faisable... ")
sol = Solution(inst, algo = :earliest)
nbviols, violtxt = get_viol_description(sol)
nblines_in_violtxt = length(findall(r".+", violtxt))
@test nbviols == 0
@test nbviols == nblines_in_violtxt
@test length(violtxt) == 0
ln1("fait.")

lg1("2. La solution viole les bornes lb et ub... ")
sol = Solution(inst, algo = :earliest)
# on impose l'avion x[3] atterrit trop tot et x[5] trop tard 
sol.x[3] = 0
sol.x[5] = 1000
nbviols, violtxt = get_viol_description(sol)
expected = "
p3 atterrit trop tôt  x=0 < lb=90
p5 atterrit trop tard ub=556 < x=1000
p2->p3 : écart insuffisant x[p2]=196->x[p3]=0 => sep(p2,p3)=-196 au lieu de 15
p1->p3 : écart insuffisant x[p1]=130->x[p3]=0 => sep(p1,p3)=-130 au lieu de 15
p5->p6 : écart insuffisant x[p5]=1000->x[p6]=235 => sep(p5,p6)=-765 au lieu de 8
p5->p7 : écart insuffisant x[p5]=1000->x[p7]=243 => sep(p5,p7)=-757 au lieu de 8
p5->p8 : écart insuffisant x[p5]=1000->x[p8]=251 => sep(p5,p8)=-749 au lieu de 8
p5->p9 : écart insuffisant x[p5]=1000->x[p9]=259 => sep(p5,p9)=-741 au lieu de 8
p5->p10 : écart insuffisant x[p5]=1000->x[p10]=267 => sep(p5,p10)=-733 au lieu de 8
"
if lg3()
    println("\nexpected=\n", expected)
    println("\nvioltxt=\n", violtxt)
end
nblines_in_violtxt = length(findall(r".+", violtxt))
@test nbviols == 9
@test nbviols == nblines_in_violtxt
ln1("fait.")

lg1("3. La solution viole les temps de séparation... ")
sol = Solution(inst, algo = :earliest)
# on impose l'avion x[3] à la valeur de x[5] => p4 et p5 deviennet faux
sol.x[3] = sol.x[5]
nbviols, violtxt = get_viol_description(sol)
nblines_in_violtxt = length(findall(r".+", violtxt))
@test nbviols == 2
@test nbviols == nblines_in_violtxt
if lg3()
    println("\nvioltxt=\n", violtxt)
end

ln1(" fait.")
