using DataStructures: OrderedDict

do_big_test = false
# do_big_test = true # fonctionne au xxxx/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
# Test de quelques exceptions ou assertionc

@test_throws AssertionError generate_mutations_shift2g(n, idx_first = 0)
@test_throws AssertionError generate_mutations_shift2g(n, idx_first = 2, idx_last = 2)
@test_throws AssertionError generate_mutations_shift2g(n, shift_max = 0)
@test_throws AssertionError generate_mutations_shift2g(n, gap_max = 0)

############################################################################
# Test des type retournés

params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 3, :gap_max => 3)
muts = generate_mutations_shift2g(n; params...)

@test typeof(muts[1]) == Mutation
# On trie pour garantir les assertions suivantes
muts = sort(muts)
@test repr(muts[1]) == "T[1,2,3,4]->[2,1,4,3]"
@test repr(muts[end]) == "T[7,8,9,10]->[8,7,10,9]"
@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test isa(muts, Vector{T} where {T<:Mutation})

ln2("L2 Mutations shift2g n=", length(muts))
ln3("L2 avec params suivants ", params)
ln4.("L4   ", muts)
@test length(muts) == 270

############################################################################
# On veut : les paire de petits shift=swap consécutifs
#    T(1,2,3,4)
#    T(2,3,4,5)
#    T(3,4,5,6)
#    T(4,5,6,7)
#    T(5,6,7,8)
#    T(6,7,8,9)
#    T(7,8,9,10)
# 
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 1, :gap_max => 1)
muts = generate_mutations_shift2g(n; params...)
ln2("L2 Mutations shift2g n=", length(muts))
ln3("L2 avec params suivants ", params)
ln3.("L3   ", muts)
@test length(muts) == 7

############################################################################
# toutes les shift2g de taille illimitée entre les indices 5 à 8 inclus
params = Dict(:idx_first => 3, :idx_last => 7, :shift_max => -1, :gap_max => 1)
muts = generate_mutations_shift2g(n; params...)

ln2("L2 Mutations shift2g n=", length(muts))
ln3("L2 avec params suivants ", params)
ln3.("L3   ", muts)
@test length(muts) == 6

############################################################################
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => -1, :gap_max => -1)
# On veut tous les shifts disjoints quelque soit leut taille

muts = generate_mutations_shift2g(n; params...)
ln2("L2 Mutations shift2g n=", length(muts))
ln3("L2 avec params suivants ", params)
ln4.("L4   ", muts)
@test length(muts) == 532


############################################################################
if do_big_test
    ln2("=========================================================")
    SIZE = 500     # taille de l'instance
    GAP_MAX = 1    # 1 si les deux shists se touchent

    specif = OrderedDict(
        # gap_max => nb_permu
        1  => 4461     , # 
        2  => 8913     , # 
        3  => 13356    , # 
        4  => 17790    , # 
        5  => 22215    , # 
        6  => 26631    , # 
        7  => 31038    , # 
        8  => 35436    , # 
        9  => 39825    , # 
        10 => 44205    , # 
        -1 => 1107817  , #  environ 8s sur paxille 3968048
         )
    # for GAP_MAX in keys(specif)
    for (GAP_MAX,n) in specif
        # n = specif[GAP_MAX]
        # ln2("L2 test BIG shift2g pour SIZE=$SIZE et GAP_MAX=$GAP_MAX")

        # ce local évite un warning d'avertissement du masquage de la global précédente
        local params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 2, :gap_max => GAP_MAX)

        local t=@elapsed   local muts = generate_mutations_shift2g(SIZE; params...)
        # ln1("   $GAP_MAX => $n  # n=",length(muts))
        ln2("L2 $(round(t,digits=2))s pour générer $(length(muts)) shift2 (ms=$(ms())s)")
        ln2("L2   avec params: ", params)
        # ln2("L2   rappel: length(muts)=",length(muts), " RESULT[GAP_MAX]=",RESULT[GAP_MAX])
        # @test length(muts) == RESULT[GAP_MAX]
        ln2("L2   rappel: length(muts)=",length(muts), " n=",n)
        @test length(muts) == n
    end
end
