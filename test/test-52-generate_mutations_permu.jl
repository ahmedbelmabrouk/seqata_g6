do_big_test = false
# do_big_test = true # fonctionne au 24/02/2022
# Log.level(3)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)

############################################################################
# Test de quelques exceptions ou assertionc

@test_throws AssertionError generate_mutations_permu(n, idx_first = 0)
@test_throws AssertionError generate_mutations_permu(n, permu_size = 1)
@test_throws AssertionError generate_mutations_permu(n, idx_first = 2, idx_last = 2)

############################################################################
params = Dict(:permu_size => 3, :idx_first => 1, :idx_last => -1)

muts = generate_mutations_permu(n; params...)

@test typeof(muts[1]) == Mutation
# Attention muts n'est pas forcément trié *selon Set ou Vector)
# @test repr(muts[1]) == "P[2,3]->[3,2]" # car permu normalisée

if lg2()
    println("L2 Mutations permu mini10 pour permu_size=3 : ", length(muts))
    if lg3()
        println("L3 avec params suivante ", params)
        lg4() && println.("L3   ", sort(muts))
    end
end
@test isa(muts, Vector{T} where {T<:AbstractMutation})
# @test length(muts) == 40 # OLD
@test length(muts) == 33

############################################################################
params = Dict(:permu_size => 2, :idx_first => 1, :idx_last => -1)

muts = generate_mutations_permu(n; params...)
if lg2()
    println("L2 Mutations permu mini10 pour permu_size=2 : ", length(muts))
    if lg3()
        println("L3 avec params suivante ", params)
        lg4() && println.("L3   ", muts)
    end
end
@test isa(muts, Vector{T} where {T<:Mutation})
@test isa(muts, Vector{T} where {T<:AbstractMutation})
@test length(muts) == 9

############################################################################
# toutes les permu de taille 3 entre les indice 5 à 8 inclus
params = Dict(:permu_size => 3, :idx_first => 5, :idx_last => 8)

muts = generate_mutations_permu(n; params...)
if lg2()
    println("L2 Mutations permu mini10 pour permu_size=3 in 5..8 : ", length(muts))
    if lg3()
        println("L3 avec params suivante ", params)
        lg4() && println.("L3   ", muts)
    end
end
@test isa(muts, Vector{T} where {T<:Mutation})
@test isa(muts, Vector{T} where {T<:AbstractMutation})
# @test length(muts) == 10   # de 5 + 5 car P[567] et P[678]
@test length(muts) == 9   # de 5 (P[567]) + 5  (P[678])- 1 (par unique!)

############################################################################
if do_big_test
    ln2("=========================================================")
    ln2("L2 BIG Mutations permu mini10 pour permu_size=9  START (ms=$(ms())s)")
    params = Dict(
        # :permu_size => 10, # 3628799 dédoublonnage **inutile** en 60s
        :permu_size => 9, # 725758 dédoublonnage **utile** en 8s (=> 685439)
        :idx_first => 1,
        :idx_last => -1,
    )

    t=@elapsed   muts = generate_mutations_permu(n; params...)
    ln2("L2 $(round(t,digits=2)) sec pour générer n=$(length(muts)) permus mini10 avec ")
    ln2("L2   avec params suivants ", params)

    @test length(muts) == 685439    # 19/06/2019 avec normalisation intégrée  

    ln4("L4 Extrait des mutations générées [1:10; 100_001:100_010; end-9:end] :")
    ln4.("L4   ", muts[[1:10; 100_001:100_010; end-9:end]])


    # GC.gc()
    t=@elapsed   muts2 = reverse.(muts)
    ln2("L2 $(round(t,digits=2)) sec pour reverse.(muts) en (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = normalize.(muts)
    ln2("L2 $(round(t,digits=2)) sec pour normalize.(muts) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   bool = allunique(muts)
    ln2("L2 $(round(t,digits=2))s pour allunique(muts)=$bool (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   unique!(muts2)
    ln2("L2 $(round(t,digits=2))s pour unique(muts) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = collect(Set(muts))
    ln2("L2 $(round(t,digits=2))s pour collect(Set(muts)) (ms=$(ms())s)")
    ln2("L2      collect(values(Set(muts)) de taille : ", length(muts2))

    # GC.gc()
    t=@elapsed   muts2 = Random.shuffle(muts)
    ln2("L2 $(round(t,digits=2))s pour Random.shuffle(muts))) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = sort(muts2)
    ln2("L2 $(round(t,digits=2))s pour sort(muts2)) (ms=$(ms())s)")    
end

############################################################################
# Log.level(4)

if do_big_test
    ln2("=========================================================")
    SIZE = 500     # taille de l'instance
    PERMU_SIZE = 8 # largeur les permus à générer

    # Attention, 
    # De plus :
    #    unique!(muts) 
    # est 2.5 fois plus long que
    #    muts = values(Set(muts))
    # L'opération @timed ou @elapsed ne ralentit pas l'exécution
    RESULT = Dict(
        # PERMU_SIZE => nb_permu
        2 => 499, # idem sans unique
        3 => 1993, # 2490 sans unique
        4 => 8951, # 11431 sans unique
        5 => 47639, # 59024 sans unique (rapide)
        6 => 297119, # 355905 sans unique ; gen:0.5s 
        7 => 2134799, # 2489266 sans unique ; gen:7.1s paxille 
        8 => 17398079, # 19877267   ; gen:106s paxille 
        9 => 158739839, # 178536468 hint:181077120; TROP LONG ! gen:xx normalize:xx unique:xx
    )
    # performance pour SIZE=8 (paxille): DUR DUR ET IRRÉGULIER !
    #     n=17398079
    #     140.18s pour génération normalisée (avec unique interne)
    #     106.07s..1147.25s pour génération normalisée (avec collect(Set(muts)))
    #     73.42s pour unique(muts)
    #     30.96s pour collect(Set(muts))
    #     4.65s pour Random.shuffle(muts)))
    #     58.6s pour sort(muts2))
    #     n=17398079
    #     51.6s pour générer 17398079 permus normalisée par collect(Set(muts))
    #     0.23s pour reverse 
    #     52.92s pour normalize.(muts)
    # performance pour SIZE=8 (marange beaucoup de RAM)): 
    #     n=17398079
    #     generate_mutations_permu sizehint=20119680 llen1=19877267 len2=17398079
    #     52.12s pour générer 17398079 permus 
    #     0.24 sec pour reverse.(muts) 
    #     53.75 sec pour normalize.(muts)
    #     8.79s pour allunique(muts)=true
    #     8.98s pour collect(Set(muts))
    #     1.06s pour Random.shuffle(muts)))
    #     39.31s pour sort(muts2))
    # performance pour SIZE=9 (marange beaucoup de RAM)): 
    #     n=158739839
    #     1143.28s génération normalisée (avec collect(Set(muts))) 

    ln2("L2 BIG PERMU sur SIZE=$SIZE pour permu_size=$PERMU_SIZE")

    params = Dict(:idx_first => 1, :idx_last => -1, :permu_size => PERMU_SIZE)

    # GC.gc()
    t=@elapsed   muts = generate_mutations_permu(SIZE; params...)
    ln2("L2 $(round(t,digits=2))s pour générer $(length(muts)) permus (ms=$(ms())s)")
    ln2("L2   avec params: ", params)
    ln2("L2   rappel: length(muts)=",length(muts), " RESULT[PERMU_SIZE]=",RESULT[PERMU_SIZE])
    # @test length(muts) == RESULT[PERMU_SIZE]

    # GC.gc()
    t=@elapsed   muts2 = reverse.(muts)
    ln2("L2 $(round(t,digits=2)) sec pour reverse.(muts) en (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = normalize.(muts)
    ln2("L2 $(round(t,digits=2)) sec pour normalize.(muts) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   bool = allunique(muts)
    ln2("L2 $(round(t,digits=2))s pour allunique(muts)=$bool (ms=$(ms())s)")

    # unique!(muts2) ou unique(muts2) semble trop lent pour grosse instances !
    # GC.gc()
    # t=@elapsed   unique!(muts2)
    # ln2("L2 $(round(t,digits=2))s pour unique(muts) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = collect(Set(muts))
    ln2("L2 $(round(t,digits=2))s pour collect(Set(muts)) (ms=$(ms())s)")
    ln2("L2      collect(values(Set(muts)) de taille : ", length(muts2))

    # GC.gc()
    t=@elapsed   muts2 = Random.shuffle(muts)
    ln2("L2 $(round(t,digits=2))s pour Random.shuffle(muts))) (ms=$(ms())s)")

    # GC.gc()
    t=@elapsed   muts2 = sort(muts2)
    ln2("L2 $(round(t,digits=2))s pour sort(muts2)) (ms=$(ms())s)")

end
