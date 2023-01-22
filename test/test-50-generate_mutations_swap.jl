do_big_test = false
# do_big_test = true # fonctionne au 24/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
# Test de quelques exceptions ou assertions

@test_throws AssertionError generate_mutations_swap(n, idx_first = 0)
@test_throws AssertionError generate_mutations_swap(n, idx_first = 2, idx_last = 2)
@test_throws AssertionError generate_mutations_swap(n, shift_min = 0)
@test_throws AssertionError generate_mutations_swap(n, shift_max = -2) # 0 serait accepté

############################################################################
params = Dict(:idx_first => 1, :idx_last => -1, :shift_min => 1, :shift_max => 5)
muts = generate_mutations_swap(n; params...) ; Base.sort!(muts)

@test typeof(muts[1]) == Mutation
@test repr(muts[1]) == "S[1,2]->[2,1]"

ln2("L2 Mutations permu mini10 pour shift_max=5 : len=", length(muts))
ln2("L3 avec params suivante ", params)
ln3.("L3   ", muts)

@test isa(muts, Vector{T} where {T<:Mutation})         # La classe réelle
@test isa(muts, Vector{T} where {T<:AbstractMutation}) # La classe mère
@test length(muts) == 35


############################################################################
params = Dict(:idx_first => 3, :idx_last => 7, :shift_min => 1, :shift_max => -1)
muts = generate_mutations_swap(n; params...) ; Base.sort!(muts)

@test typeof(muts[1]) == Mutation
@test repr(muts[1]) == "S[3,4]->[4,3]"

ln2("L2 Mutations swap mini10 entre 3 et 7 : len=", length(muts))
ln2("L3 avec params suivante ", params)
ln3.("L3   ", muts)
@test length(muts) == 10

############################################################################
# TEST DE swap avec shift_min

params = Dict(:idx_first => 1, :idx_last => -1, :shift_min => 7, :shift_max => -1)
muts = generate_mutations_swap(n; params...) ; Base.sort!(muts)

@test typeof(muts[1]) == Mutation
@test repr(muts[1]) == "S[1,8]->[8,1]"

ln2("L2 Mutations swap mini10 avec shift_min=7 len=", length(muts))
ln2("L3 avec params suivante ", params)
ln4.("L4   ", muts)

@test length(muts) == 6

############################################################################
############################################################################
############################################################################
### TEST DE PERFORMANCE POUR LES MUTATIONS DE CLASSE swap

############################################################################
# TEST GROSSE INSTANCE  (création de inst et de sol initiale)

INST="../data/09.alp" ; INST_NB_MUTS=4950    # shift
INST="../data/10.alp" ; INST_NB_MUTS=11175   # shift
INST="../data/11.alp" ; INST_NB_MUTS=19900   # shift
INST="../data/12.alp" ; INST_NB_MUTS=31125   # shift
INST="../data/13.alp" ; INST_NB_MUTS=124750; # shift


ln2("Lecture $(INST), et création sol ($(ms())ms) ...")
sol_ori = Solution(Instance(INST)) #  premier appel => précompilation
sol = Solution(sol_ori)
n = length(sol.planes)
ln2("Lecture $(INST), et création sol_ori et copie sol ($(ms())ms) FAIT.")

@test repr(sol) == repr(sol_ori)
# Plus tard dans le code pour restaurer sol on fera : copy!(sol, sol_ori)

############################################################################
# TEST GROSSE INSTANCE 1
if do_big_test
    ln2("L2 BIG Mutations swap mini10 pour permu_size=10 ")
    GC.gc()
    ln2("L2 BIG Mutations shift pour shift_max=10 ; START car do_big_test=true")
    params = Dict(:idx_first => 1, :idx_last => -1, :shift_min => 1, :shift_max => -1)
    t=@elapsed   muts = generate_mutations_swap(n; params...)
    ln2("L2 $(round(t,digits=2))sec pour générer $(length(muts)) swaps (ms=$(ms())s)")
    ln2("L2     avec params suivants : ", params)
    @test length(muts) == INST_NB_MUTS
end
