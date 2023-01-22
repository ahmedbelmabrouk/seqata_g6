do_big_test = false
# do_big_test = true # fonctionne au 25/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
params = Dict(
    :idx_first => 1,
    :idx_last => 5,
    :add_shift => true,
    :add_shift2 => false,
    :add_shift2g => false,
    :add_swap => true,
    :add_permu => false,
    :shift_min => 1,
    :shift_max => -1,
    # :gap_max => -1,
    # :permu_size => 4,
)

muts = generate_neighborhood(n; params...)

ln2("L2 Mutations mini10 pour add_shift et add_swap : ", length(muts))
ln2("L2 avec params suivants ", params)
ln3.("L3   ", muts)

@test length(muts) == 22
@test length(muts) == length(Set(muts))

############################################################################
params = Dict(
    :idx_first => 1,
    :idx_last => 5,
    :add_shift => true,
    :add_shift2 => false,
    :add_shift2g => true,
    :add_swap => true,
    :add_permu => true,
    :shift_min => 1,
    :shift_max => 3,
    :gap_max => -1,
    :permu_size => 3,
)

muts = generate_neighborhood(n; params...)

ln2("L2 Mutations toutes classes sauf shift2 entre 1 et 5 :", length(muts))
ln3("L2 avec params suivants ", params)
ln3.("L3   ", muts)
ln3("L3 Même muts mais triées ")
ln3.("L3   ", sort(muts))

# # @test length(muts) == 41 # OLD
# @test length(muts) == 39
# unique!(muts)
@test length(muts) == 26

# @test length(muts) == length(Set(muts))

############################################################################
# Log.level(3)
params = Dict(
    :idx_first => 1,
    :idx_last => 5,
    :add_shift => false,
    :add_shift2 => true,
    :add_shift2g => true,
    :add_swap => false,
    :add_permu => false,
    :shift_min => 1,
    :shift_max => 3,
    :gap_max => -1,
    :permu_size => 3,
)

muts = generate_neighborhood(n; params...)

ln2("L2 Mutations classes shift2 et shift2g entre 1 et 5 : ", length(muts))
ln3("L2 avec params suivants ", params)
ln3.("L3   ", muts)
ln3("L3 Même muts mais triées ")
ln3.("L3   ", sort(muts))

@test length(muts) == 44

############################################################################
# Log.level(1)
params = Dict(
    :idx_first => 1,
    :idx_last => 6,
    :add_shift => true,
    :add_shift2 => true,
    :add_shift2g => true,
    :add_swap => true,
    :add_permu => true,
    :shift_min => 1,
    :shift_max => 4,
    :gap_max => -1,
    :permu_size => 4,
)

muts = generate_neighborhood(n; params...)
Base.sort!(muts)

ln2("L2 Mutations TOUTES CLASSES entre 1 et 6 et shift_max=4 : ", length(muts))
ln3("L2 avec params suivants ", params)
ln3.("L3   ", muts)
ln3("L3 Même muts mais triées (était unique)")
# ln3.("L3   ", sort(muts))
ln3.("L3   ", muts)
@test length(muts) == 148 # unique pour [1..6] et shift_max=> 4

############################################################################
# Log.level(3)
if do_big_test
    ln2("L2  BIG generate_neighborhood")
    GC.gc()
    params = Dict(
        :idx_first => 1,
        :idx_last => -1,
        :add_shift => true,
        :add_shift2 => true,
        :add_shift2g => true,
        :add_swap => true,
        :add_permu => true,
        :shift_min => 1,
        :shift_max => -1,
        :gap_max => -1,
        :permu_size => -1,
    )

    ln2("L2 params pour BIG MUTS ", params)
    ln2("L2 calcul de BIG MUTS pour toutes classes de taille illimitée... ")

    t=@elapsed   muts = generate_neighborhood(n; params...) # e.g 84s ; 82.9s
    ln2("=> generate_neighborhood FAIT : len=$(length(muts))  DUREE=$(round(t,digits=2))s")

    # # @test length(muts) == 3629448 # sans add_shift2
    # @test length(muts) == 3632050 # avec add_shift2
    @test length(muts) == 3628799 # avec add_shift2

    lg2("L2 collect(values(Set(muts))) du big muts... ")
    # collect(values(Set(muts))) est plus rapide que unique!(muts) : 14.3s vs 22.5s
    # mais ne résultat n'est pas trié.
    # unique!(muts)
    t=@elapsed   muts = collect(values(Set(muts)))
    ln2("FAIT : DUREE=$(round(t,digits=2))s")

    # @test length(muts) == 3628799 # sans add_shift2
    @test length(muts) == 3628799 # avec add_shift2

    lg2("L2 shuffle du big muts... ")
    t=@elapsed   muts = shuffle(muts)
    ln2("FAIT : DUREE=$(round(t,digits=2))s")

    lg2("L2 sort du big muts... ")
    t=@elapsed   muts = sort(muts)
    ln2("FAIT : DUREE=$(round(t,digits=2))s")

    lg2("L2 Conversion de big muts en Set... ")
    t=@elapsed   mutset = Set(muts)
    ln2("FAIT : DUREE=$(round(t,digits=2))s")

    lg2("L2 Reconversion du mutset en Vector... ")
    t=@elapsed   set2 = values(mutset)
    ln2("FAIT : DUREE=$(round(t,digits=2))s")
end
