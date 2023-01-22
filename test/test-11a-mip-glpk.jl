# ===========
# Construction de l'instance mini
inst = instance_build_mini10()
@test inst.nb_planes == 10

# ===========
lg1("1. Création d'un MipSolver pour GLPK... ")
Args.set("external_lp_solver", :glpk)
solver = MipSolver(inst)

@test Symbol(lowercase(solver_name(solver.model))) == :glpk

ln1(" fait.")

# ===========
lg1("2. Résolution par MipSolver pour GLPK (cost=700.0 ?)... ")
solve!(solver)
bestsol = solver.bestsol
@test bestsol.cost == 700.0
ln1(" (cost=$(solver.bestsol.cost)) fait.")
