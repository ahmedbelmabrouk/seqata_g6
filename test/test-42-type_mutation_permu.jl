
m = Mutation([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])

@test typeof(m.indices1) == Vector{Int}
@test m.indices1 == Int[1, 2, 4, 5]
@test typeof(m.indices2) == Vector{Int}
@test m.indices2 == Int[5, 4, 2, 1]
@test m.class == :permu

@test repr(m) == "P[1,2,4,5]->[5,4,2,1]"
@test repr(reverse(m)) == "P[5,4,2,1]->[1,2,4,5]"


##############################################################
m = Mutation([1, 2, 5, 9, 10], [2, 1, 5, 10, 9])
@test repr(reverse(m)) == "P[2,1,10,9]->[1,2,9,10]"

# Test de l'application de cette mutation à un solution
sol_ori = Solution(instance_build_mini10())
@test occursin("[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]", repr(sol_ori))

sol = Solution(sol_ori)
mutate!(sol, m)
@test occursin("[p2,p1,p3,p4,p5,p6,p7,p8,p10,p9]", repr(sol))

# La mutation inverse doit redonner la solution originale :
mutate!(sol, reverse(m))
@test repr(sol) == repr(sol_ori)


# Test d'une mutation nulle
m = Mutation([1, 2, 10], [1, 2, 10])
@test repr(m) == "P[]->[]"

# L'application de la mutation nulle ne doit pas change la solution initiale
sol = Solution(sol_ori)
mutate!(sol, m)
@test repr(sol) == repr(sol_ori)

#
# Trois mutations identiques, test de la normalisation :
#

m1 = Mutation([4, 5, 10], [10, 4, 5])
m2 = Mutation([10, 1, 2, 4, 5], [5, 1, 2, 10, 4])
m3 = Mutation([10, 1, 2, 4, 5], [5, 1, 2, 10, 4], do_normalize=false)

# Test le la normalisation des permutations
@test width(m1) == 7 # expected 7 == 10-4+1
@test length(m1) == 3 # expected 3 car trois éléments permutent
@test repr(m1) == "P[4,5,10]->[10,4,5]"
@test repr(m2) == "P[4,5,10]->[10,4,5]"
@test repr(m1) == repr(m2)
@test repr(m3) == "P[10,1,2,4,5]->[5,1,2,10,4]" # car non normalisée


############################################################################
# Test de contruction de mutations avec des indices

# 1 et 2 swaps
m = Mutation(:swap , 9, 10)
@test repr(m) == "S[9,10]->[10,9]"

m = Mutation(:swap , 9, 10, 1, 2)
@test repr(m) == "S[1,2,9,10]->[2,1,10,9]"

# 2 swaps superposé pour faire une rotation de trois éléments
m = Mutation(:swap , 5, 1, 10, 5)
@test repr(m) == "S[1,5,10]->[5,10,1]"

sol_ori = Solution(instance_build_mini10())
sol = Solution(sol_ori)
@test occursin("[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]", repr(sol))
mutate!(sol, m)
@test occursin("[p5,p2,p3,p4,p10,p6,p7,p8,p9,p1]", repr(sol))


# 2 swaps inverses qui s'annulent
m = Mutation(:swap , 5, 1, 5, 1)
@test repr(m) == "S[]->[]"

sol_ori = Solution(instance_build_mini10())
sol = Solution(sol_ori)
@test occursin("[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]", repr(sol))
mutate!(sol, m)
@test occursin("[p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]", repr(sol))


# 1 et 2 shift
m = Mutation(:shift , 10, 8)
@test repr(m) == "T[8,9,10]->[10,8,9]"

m = Mutation(:shift , 10, 8, 1, 3)
# [1,2,3,8,9,10] -> [1,2,3,10,8,9] -> [2,3,1,10,8,9]
@test repr(m) == "T[1,2,3,8,9,10]->[2,3,1,10,8,9]"

# 2 shift inverses qui s'annulent
m = Mutation(:shift , 10, 8, 8, 10)
@test repr(m) == "T[]->[]"

############################################################################
# Test du merge d'une mutations avec des indices

m = Mutation([1, 2, 3, 4, 5], [2, 3, 4, 5, 1])
@test repr(m) == "P[1,2,3,4,5]->[2,3,4,5,1]"

m2 = merge_mutation(m, :swap, 1, 3)
@test repr(m2) == "X[1,2,3,4,5]->[4,3,2,5,1]"

m2 = merge_mutation(m, :swap, 8, 10)
@test repr(m2) == "X[1,2,3,4,5,8,10]->[2,3,4,5,1,10,8]"

m2 = merge_mutation(m, :shift, 1, 3)
@test repr(m2) == "X[1,2,3,4,5]->[3,4,2,5,1]"

m2 = merge_mutation(m, :shift, 8, 10)
@test repr(m2) == "X[1,2,3,4,5,8,9,10]->[2,3,4,5,1,9,10,8]"

############################################################################
# Test du merge deux mutations

m = Mutation([1, 2, 3, 4, 5], [2, 3, 4, 5, 1])
@test repr(m) == "P[1,2,3,4,5]->[2,3,4,5,1]"

# m2 = merge_mutation(m, :swap, 1, 3)
m2 = merge_mutation(m, Mutation(:swap, 1, 3))
@test repr(m2) == "X[1,2,3,4,5]->[4,3,2,5,1]"

# m2 = merge_mutation(m, :swap, 8, 10)
m2 = merge_mutation(m, Mutation(:swap, 8, 10))
@test repr(m2) == "X[1,2,3,4,5,8,10]->[2,3,4,5,1,10,8]"

# m2 = merge_mutation(m, :shift, 1, 3)
m2 = merge_mutation(m, Mutation(:shift, 1, 3))
@test repr(m2) == "X[1,2,3,4,5]->[3,4,2,5,1]"

# m2 = merge_mutation(m, :shift, 8, 10)
m2 = merge_mutation(m, Mutation(:shift, 8, 10))
@test repr(m2) == "X[1,2,3,4,5,8,9,10]->[2,3,4,5,1,9,10,8]"
