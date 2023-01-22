do_big_test = false
# do_big_test = true # fonctionne au 24/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
# Test de quelques exceptions ou assertionc

@test_throws AssertionError generate_mutations_shift(n, idx_first = 0)
@test_throws AssertionError generate_mutations_shift(n, idx_first = 2, idx_last = 2)
@test_throws AssertionError generate_mutations_shift(n, shift_max = 0)

############################################################################
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 4)

muts = generate_mutations_shift(n; params...)

@test typeof(muts[1]) == Mutation
Base.sort!(muts)
@test repr(muts[1]) == "T[1,2]->[2,1]"

ln2("L2 Mutations shift mini10 pour shift_max=4 : ", length(muts))
ln2("L3 avec params suivante ", params)
ln4.("L4   ", muts)

@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test length(muts) == 51

############################################################################
params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 1)


muts = generate_mutations_shift(n; params...)
ln2("L2 Mutations shift pour shift_max=1 : ", length(muts))
ln2("L3 avec params suivante ", params)
ln4.("L4   ", muts)

@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test length(muts) == 9

############################################################################
# toutes les shift de taille 3 entre les indices 5 à 8 inclus
params = Dict(:idx_first => 5, :idx_last => 8, :shift_max => -1)


muts = generate_mutations_shift(n; params...)
ln2("L2 Mutations shift pour shift_max=3 in 5..8 : ", length(muts))
ln2("L3 avec params suivante ", params)
ln4.("L4   ", muts)

@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test length(muts) == 9

############################################################################
# toutes les shift d'écart minimum de 5
params = Dict(:idx_first => 1, :idx_last => 10, :shift_min => 6, :shift_max => -1)


muts = generate_mutations_shift(n; params...)
ln2("L2 Mutations shift pour shift_min=6 in 1..10 : ", length(muts))
ln2("L3 avec params suivante ", params)
ln4.("L4   ", muts)

@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test length(muts) == 20

############################################################################
if do_big_test
    ln2("L2 BIG Mutations shift pour shift_max=10")
    params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => 10)

    t=@elapsed   muts = generate_mutations_shift(n; params...)
    ln2("L2 $(round(t,digits=2)) sec pour générer n=$(length(muts)) shift mini10 avec ")
    lg2("L2   ") && @show (t, bytes, gctime, memallocs)
    ln2("L2   avec params suivants ", params)

    @test length(muts) == 81   # (n-1)*(n-1)
    ln4.("L4   ", muts)
end
