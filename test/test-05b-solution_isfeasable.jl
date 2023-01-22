# Test de la méthode is_feasable(sol::Solution)

# PARTIE 1
# Test sur l'instance instance_build_mini10 (instance alp01)
# Comme elle reste faisable quelque soit la permutation, on force
# l'infaisablilités en modifiant les caractéristiques d'un avion 
# dans l'instance
# 

lg1("faisabilité d'une solution PARTIE 1 : test sur alp01... ")
ln2("\nGénération d'une instance instance_build_mini10.")

inst = instance_build_mini10()
@test inst.nb_planes == 10

# Création de la solution vide
sol = Solution(inst)

ln2("La solution initiale de alp10 doit être faisable.")
@test is_feasable(sol) == true

if lg3()
    println(to_s_long(sol))
end

###############################
# on rend la sol infaisable en imposant le UB du troisième avion égale 
# au LB du premier avion (i.e ON MODIFIE L'INSTANCE) !
ln2("On modifie alp10 par p3.ub := p1.lb pour rendre infaisable.")

# Mais les avions sont immutables, on doit donc créer un nouvel avion `bad_plane`
# par copie de planes[3] qui remplecera l'avion d'origine
p = sol.planes[3]
bad_ub = sol.planes[1].lb
bad_plane = Plane(
    id = p.id,
    name = p.name,
    kind = p.kind,
    at = p.at,
    lb = p.lb,
    target = p.target,
    ub = bad_ub,   # UB est modifié !
    ep = p.ep,
    tp = p.tp,
)

# sol.planes[3].ub = sol.planes[1].lb
sol.planes[3] = bad_plane
@test is_feasable(sol) == false
if lg3()
    println(to_s_long(sol))
end

ln1("fait.")

######################################################
# PARTIE 2
# Test sur l'instance alp_07 car très peu de solutions faisables.
# PRINCIPE
# - on teste l'instance trié par target croissant 
#   => doit être faisable
# - on teste l'instance trié par target décroissant 
#   => doit être infaisable
# - on teste NB_SHUFFLE mélanges (e.g 10_000 mélanges)
# - on doit trouver emte 97% et 98% de solutions infaisable !
# 

# inst = instance_build_mini10()
# inst = Instance("data/07.alpx")

lg1("faisabilité d'une solution PARTIE 2 : test sur alp07... ")

DATA_DIR = dirname(@__FILE__) * "/data/"
NB_SHUFFLES = 10_000

inst = Instance("$DATA_DIR/alp_07.alp") # instance SEQATA ou ALAP
# inst = Instance("$DATA_DIR/alp_07.alp") # instance ALAP

@test inst.nb_planes == 44

# Création de la solution vide
sol = Solution(inst)

initial_sort!(sol, presort = :target)
@test is_feasable(sol) == true

initial_sort!(sol, presort = :rtarget)
@test is_feasable(sol) == false


nb_infeasable = 0
if lg3()
    println("\nAvant $(NB_SHUFFLES) mélanges ")
end

for i = 1:NB_SHUFFLES
    shuffle!(sol)
    if !is_feasable(sol)
        global nb_infeasable += 1
        if lg4()
            println("shuffle! NOT FEASABLE : ", sol)
        end
    end
end
if lg2()
    println(
        "Après $(NB_SHUFFLES) mélanges => nb_infeasable=",
        nb_infeasable,
        " soit ",
        nb_infeasable / NB_SHUFFLES,
        " (test in [0.97, 0.982] ?",
    )
end
@test 0.97 <= nb_infeasable / NB_SHUFFLES <= 0.982
ln1("fait.")
