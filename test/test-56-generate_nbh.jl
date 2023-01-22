do_big_test = false
do_big_test = true # fonctionne au 25/02/2022
# Log.level(3)
# Log.level(2)

sol = Solution(instance_build_mini10())
@test sol.cost == 25650.0
n = length(sol.planes)


############################################################################
ln2("=== test s1,t1,S1,T1 AssertionError")
@test_throws AssertionError generate_nbh(n, "s1")
@test_throws AssertionError generate_nbh(n, "t1")
@test_throws AssertionError generate_nbh(n, "S1")
@test_throws AssertionError generate_nbh(n, "T1")


############################################################################
ln2("=== test s2 et t2 (mutations identiques)")
muts_s2,label_s2 = generate_nbh(n, "s2") ; Base.sort!(muts_s2)
muts_t2,label_t2 = generate_nbh(n, "t2") ; Base.sort!(muts_t2)
@test label_s2 == "s2"
@test label_t2 == "t2"
@test isa(muts_s2, Vector{T} where {T<:Mutation})

ln2("    len $label_s2: ",length(muts_s2))
@test length(muts_s2) == 9

ln2("    len $label_t2: ",length(muts_t2))
@test length(muts_t2) == 9

# On doit trier car l'ordre n'est pas garanti
@test repr(muts_s2[1]) == "S[1,2]->[2,1]"
@test repr(muts_t2[1]) == "T[1,2]->[2,1]"

muts_st2,label_st2 = generate_nbh(n, "s2+t2") ; Base.sort!(muts_st2)

ln2("    len $label_st2: ",length(muts_st2))
@test length(muts_st2) == 9

ln3("Valeurs de muts_s2 :")
ln3.("   ", muts_s2)
ln3("Valeurs de muts_t2 :")
ln3.("   ", muts_t2)
ln3("Valeurs de muts_st2 :")
ln3.("   ", muts_st2)


############################################################################
ln2("=== test s6")
muts,label = generate_nbh(n, "s6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 5
ln3.("   ", muts)

############################################################################
ln2("=== test S6")
muts,label = generate_nbh(n, "S6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 35
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test t6")
muts,label = generate_nbh(n, "t6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 10
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test T6")
muts,label = generate_nbh(n, "T6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 61
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test P3")
muts,label = generate_nbh(n, "P3") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 33
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test d6")
muts,label = generate_nbh(n, "d6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 15
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test D6")
muts,label = generate_nbh(n, "D6") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 87
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test s2+s3+s4+s5+s6+t2+t3+t4+t5+t6 (+P3) est équivalent à D6")
muts,label = generate_nbh(n, "s2+s3+s4+s5+s6+t2+t3+t4+t5+t6+P3") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 87
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)

############################################################################
ln2("=== test 2T4")
muts,label = generate_nbh(n, "2T4") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 124
# ln3.("   ", muts[[1:3;end-3:end]])
ln2.("   ", muts)


############################################################################
ln2("=== test 2T3g2")
muts,label = generate_nbh(n, "2T3g2") ; Base.sort!(muts)
ln2("    len $label: ",length(muts))
@test length(muts) == 93
# ln3.("   ", muts[[1:3;end-3:end]])
ln2.("   ", muts)

############################################################################
ln2("=== test p4 **À FAIRE** (pour n=6)")
# muts,label = generate_nbh(6, "p4") ; Base.sort!(muts)
# ln2("    len $label: ",length(muts))
# @test length(muts) == 93
# # ln3.("   ", muts[[1:3;end-3:end]])
# ln2.("   ", muts)

############################################################################
ln2("=== test p5.4 **À FAIRE** (pour n=6)")
# muts,label = generate_nbh(6, "p5.4") ; Base.sort!(muts)
# ln2("    len $label: ",length(muts))
# @test length(muts) == 93
# # ln3.("   ", muts[[1:3;end-3:end]])
# ln2.("   ", muts)


############################################################################
ln2("=== test D2+2T2+2T2g1+P3")
muts,label = generate_nbh(n, "D2+2T2+2T2g1+P3"); Base.sort!(muts)
ln2("D2+2T2+2T2g1+p3 length=", length(muts))
ln2("    len $label: ",length(muts))
# @test length(muts) == 40
@show typeof(muts), label
# ln3.("   ", muts[[1:3;end-3:end]])
ln3.("   ", muts)




############################################################################
if do_big_test
    ln2("=========================================================")
    # SIZE = 500     # taille de l'instance
    # SHIFT_MAX = 13  # xxx

    # RESULT = Dict(
    #     # shift_max => nb_permu
    #     1  => 499,    #  
    #     1  => 499,    #       
    #     2  => 1993,   #      
    #     3  => 8454,   #  
    #     4  => 27302,  #  
    #     5  => 63932,  #  
    #     6  => 124200, #  
    #     7  => 213926, #  
    #     8  => 338894, #   8.6s 
    #     9  => 504852, #  (504852 avant unique!)     14.6s
    #     10 => 717512, #  (717512 avant unique!)    29.2s
    #     11 => 982550, #  (982550 avant unique sizehint 8520336)  42.74s
    #     12 => 1305606, # (1305606 avant unique, sizzehint=11875968)    62.8s (paxille)
    #     13 => 1692284, # (1692284 avant unique!)    153s sur paxille
    # )
    # ln2("L2 test BIG shift2 pour SIZE=$SIZE et SHIFT_MAX=$SHIFT_MAX")

    # # ce local évite un warning d'avertissement du masquage de la global précédente
    # local params = Dict(:idx_first => 1, :idx_last => -1, :shift_max => SHIFT_MAX)
    # t=@elapsed   local muts = generate_mutations_shift2(SIZE; params...)
    # ln2("L2 $(round(t,digits=2))s pour générer $(length(muts)) shift2 (ms=$(ms())s)")
    # ln2("L2   avec params: ", params)
    # ln2("L2   rappel: length(muts)=",length(muts), " RESULT[SHIFT_MAX]=",RESULT[SHIFT_MAX])
    # @test length(muts) == RESULT[SHIFT_MAX]
end
