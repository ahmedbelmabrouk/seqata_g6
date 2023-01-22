# ===========
# Construction de l'instance mini
inst = instance_build_mini10()
@test inst.nb_planes == 10

# ===========
lg1("1. Création d'un MipSolver pour Cbc... ")
Args.set("external_lp_solver", :cbc)
solver = MipSolver(inst)

# Le nom du solver Cbc n'est pas standard :
#   @show solver_name(solver.model)
#   => "COIN Branch-and-Cut (Cbc)"
# Donc ceci ne fonctionne pas :
#   @test Symbol(lowercase(solver_name(solver.model))) == :cbc
# Je teste alors la présence de la sous-chaine cbc
@test occursin("cbc", lowercase(solver_name(solver.model)))
ln1(" fait.")

# ===========
lg1("2. Résolution par MipSolver pour Cbc (cost=700.0 ?)... ")
solve!(solver)
bestsol = solver.bestsol
@test bestsol.cost == 700.0
ln1(" (cost=$(solver.bestsol.cost)) fait.")
