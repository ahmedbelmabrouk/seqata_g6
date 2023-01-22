# Test de ExploreSolver

# ===========
# Construction de l'instance mini
inst = instance_build_mini10()
@test inst.nb_planes == 10

# ===========
lg1("1. Création d'un ExploreSolver... ")

Log.pushlevel!(0) # on passe temporairement à 0
solver = ExploreSolver(inst)
Log.poplevel!()

ln1("fait.")

# ===========
lg1("2. Mélange de la solution initiale... ")
initial_sort!(solver.cursol, presort = :shuffle)
copy!(solver.bestsol, solver.cursol)
ln1(" => (initcost=$(solver.cursol.cost)) fait")

# ===========
itermax = 10_000
lg1("3. Résolution avec itermax=$itermax) itérations ")
lg1("(cost=700 <= 5500.0 ?)... ")

Log.pushlevel!(0) # on passe temporairement à 0
solve!(solver, itermax)
Log.poplevel!()

bestsol = solver.bestsol
@test bestsol.cost <= 5500.0
ln1(" => ok ($(bestsol.cost))")

# ===========
itermax = 100_000
lg1("4. Résolution avec itermax=$itermax) itérations supplémentaires")
lg1("(cost=700 <= 5500.0 ?)... ")

Log.pushlevel!(0) # on passe temporairement à 0
solve!(solver, itermax)
Log.poplevel!()

bestsol = solver.bestsol
@test bestsol.cost <= 5500.0
ln1(" => ok ($(bestsol.cost))")

#./
