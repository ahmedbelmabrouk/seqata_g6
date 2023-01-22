do_big_test = false
# do_big_test = true # fonctionne au 25/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
# Test de quelques exceptions ou assertionc

lg2() && println("L2 test des exceptions levés par generate_mutations_shift2")

@test_throws AssertionError generate_mutations_shift2(n, idx_first = 0)
@test_throws AssertionError generate_mutations_shift2(n, idx_first = 2, idx_last = 2)
@test_throws AssertionError generate_mutations_shift2(n, shift_max = 0)

############################################################################
ln2("L2 test shift2 sur alp10 avec shift_max=3")
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 3)
# On veut :

mutations = generate_mutations_shift2(n; params...)

# Pour ce premier test, on trie des mutation pour simplifier la vérif
Base.sort!(mutations)

@test typeof(mutations[1]) == Mutation
@test repr(mutations[1]) == "T[1,2]->[2,1]"
@test repr(mutations[end]) == "T[9,10]->[10,9]"

ln2("L2 Mutations shift2 pour shift_max=3  : ", length(mutations))
ln2("L2 avec params suivants ", params)
ln5.("L5   ", mutations)

@test isa(mutations, Vector{T} where {T<:Mutation})
@test length(mutations) == 124


############################################################################
ln2("L2 test shift2 sur alp10 avec shift_max=-1")
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => -1)

mutations = generate_mutations_shift2(n; params...)

# Pour ce premier test, on trie ses mutations pour simplifier la vérif
Base.sort!(mutations)

@test typeof(mutations[1]) == Mutation
@test repr(mutations[1]) == "T[1,2]->[2,1]"
@test repr(mutations[end]) == "T[9,10]->[10,9]"

ln2("L2 Mutations shift2 pour shift_max=4 : ", length(mutations))
ln2("L2 avec params suivants ", params)
ln5.("L5   ", mutations)

@test isa(mutations, Vector{T} where {T<:Mutation})
@test length(mutations) == 2602

############################################################################
lg2() && println("L2 test shift2 sur alp10 avec shift_max=1")
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 1)
# On veut :
#   T[1,2]->[2,1]
#   T[2,3]->[3,2]
#   T[3,4]->[4,3]
#   T[4,5]->[5,4]
#   T[5,6]->[6,5]
#   T[6,7]->[7,6]
#   T[7,8]->[8,7]
#   T[8,9]->[9,8]
#   T[9,10]->[10,9]
# 

mutations = generate_mutations_shift2(n; params...)
ln2("L2 test shift2 sur alp10 avec shift_max=1 : ", length(mutations))
ln2("L2 avec params suivants ", params)
ln5.("L5   ", mutations)

@test isa(mutations, Vector{T} where {T<:Mutation})
@test length(mutations) == 9

############################################################################
lg2() && println("L2 test shift2 sur alp10 avec shift_max=2 (len=208->33 via unique!)")
params = Dict(
    :idx_first => 1,
    :idx_last => -1,
    :shift_max => 2, # donc largeur de win = 3
)
# On veut :
#   T[1,2]->[2,1]
#   T[1,2,3]->[2,3,1]
#   T[1,2,3]->[3,1,2]
#   T[1,3]->[3,1]     # idem permution [1,2,3]->[3,2,1]
#   T[2,3]->[3,2]
#   T[2,3,4]->[3,4,2]
#   T[2,3,4]->[4,2,3]
#   T[2,4]->[4,2]
#   T[3,4]->[4,3]
#   ...
#   T[8,10]->[10,8]
#   T[9,10]->[10,9]
# 

mutations = generate_mutations_shift2(n; params...)
Base.sort!(mutations)
ln2("L2 test shift2 sur alp10 avec shift_max=2 : ", length(mutations))
ln2("L2 avec params suivants ", params)
ln5.("L5   ", mutations)
@test isa(mutations, Vector{T} where {T<:Mutation})
@test length(mutations) == 33

############################################################################
lg2() && println("L2 test shift2 sur alp10 pour [5:8] et shift_max=-1")
# toutes les shift2 de taille illimitée entre les indices 5 à 8 inclus
params = Dict(:idx_first => 5, :idx_last => 8, :shift_max => -1)

# On veut :
#   T[5,6]->[6,5]
#   T[5,6,7]->[6,7,5]
#   T[5,6,7]->[7,5,6]
#   T[5,6,7,8]->[6,5,8,7]
#   T[5,6,7,8]->[6,7,8,5]
#   T[5,6,7,8]->[6,8,5,7]
#   T[5,6,7,8]->[7,5,8,6]
#   T[5,6,7,8]->[7,8,5,6]
#   T[5,6,7,8]->[7,8,6,5]
#   T[5,6,7,8]->[8,5,6,7]
#   T[5,6,7,8]->[8,7,5,6]
#   T[5,6,8]->[6,8,5]
#   T[5,6,8]->[8,5,6]
#   T[5,7]->[7,5]
#   T[5,7,8]->[7,8,5]
#   T[5,7,8]->[8,5,7]
#   T[5,8]->[8,5]
#   T[6,7]->[7,6]
#   T[6,7,8]->[7,8,6]
#   T[6,7,8]->[8,6,7]
#   T[6,8]->[8,6]
#   T[7,8]->[8,7]
#

mutations = generate_mutations_shift2(n; params...)
if lg2()
    println("L2 Mutations shift2 pour dist[12]_max=-1 in 5..8 : ", length(mutations))
    if lg3()
        println("L3 avec params suivants ", params)
        lg5() && println.("L5   ", sort(mutations))
    end
end
@test isa(mutations, Vector{T} where {T<:Mutation})
@test length(mutations) == 22

############################################################################
if do_big_test
    ln2("=========================================================")
    SIZE = 500     # taille de l'instance
    SHIFT_MAX = 13  # xxx

    RESULT = Dict(
        # shift_max => nb_permu
        1  => 499,    #  
        1  => 499,    #       
        2  => 1993,   #      
        3  => 8454,   #  
        4  => 27302,  #  
        5  => 63932,  #  
        6  => 124200, #  
        7  => 213926, #  
        8  => 338894, #   8.6s 
        9  => 504852, #  (504852 avant unique!)     14.6s
        10 => 717512, #  (717512 avant unique!)    29.2s
        11 => 982550, #  (982550 avant unique sizehint 8520336)  42.74s
        12 => 1305606, # (1305606 avant unique, sizzehint=11875968)    62.8s (paxille)
        13 => 1692284, # (1692284 avant unique!)    153s sur paxille
    )
    ln2("L2 test BIG shift2 pour SIZE=$SIZE et SHIFT_MAX=$SHIFT_MAX")

    # ce local évite un warning d'avertissement du masquage de la global précédente
    local params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => SHIFT_MAX)
    t=@elapsed   local muts = generate_mutations_shift2(SIZE; params...)
    ln2("L2 $(round(t,digits=2))s pour générer $(length(muts)) shift2 (ms=$(ms())s)")
    ln2("L2   avec params: ", params)
    ln2("L2   rappel: length(muts)=",length(muts), " RESULT[SHIFT_MAX]=",RESULT[SHIFT_MAX])
    @test length(muts) == RESULT[SHIFT_MAX]
end
