# Fichier de test de manipulation des forme d'instanve ORLIB, AMPL, ALP et ALPX
# La version pour SEQATA est simplifiée pour ne gérer que le seul format "alp"


# =========== LECTURE FORMAT alp
format = "alp"
relfile = "test/data/alp_01.alp"
exp_inst_name = "alp_01_p10"

lg1("Lecture du fichier au format $format $relfile...")
inst = Instance("$APPDIR/$relfile")
@test inst.name == exp_inst_name
@test length(inst.planes) == 10
@test repr(inst.sep_mat) == "[3 15; 15 8]"
lg1(" fait\n")


# =========== EXPORT AU FORMAT alpx, alp, ampl
relfile = "test/data/alp_01.alp"
inst = Instance("$APPDIR/$relfile")


# export alp
lg1("Export au formats alp... ")
alp_str = to_s_long(inst, format = "alp")
ln3("alp_str:\n$alp_str")
# @test length(alp_str) == 753
# println(alp_str)
# @test length(alp_str) == 772
@test length(alp_str) == 750
@test startswith(alp_str[29:43], "name alp_01_p10")
# Format alp => les pénalités sont simple (ep tp)
@test occursin(r" ep +tp ", alp_str)
lg1(" fait\n")







# ===========
relfile = "test/data/alp_01.alp"
lg1("Lecture du fichier au format alp $relfile...")
inst = Instance("$APPDIR/$relfile")
@test inst.name == "alp_01_p10"
alp_alpx_str = to_s_long(inst)
# println(alp_alpx_str)
# @test length(alp_alpx_str) == 1164  # FORMAT ALPX (pour alap)
@test length(alp_alpx_str) == 750     # FORMAT ALP (pour seqata)
@test alp_alpx_str[29:43] == "name alp_01_p10"
lg1(" fait\n")

# =========== Statistiques sur l'instance
lg1("Statistiques sur l'instance $(inst.name)... ")
@test inst.nb_planes == length(inst.planes)
@test inst.nb_planes == 10
@test inst.nb_kinds == 2
@test inst.name == "alp_01_p10"
lg1(" fait\n")

# ===========
lg1("Lecture et test du premier avion (p1)... ")
p1 = inst.planes[1]
p2 = inst.planes[2]
p3 = inst.planes[3]
@test p1.name == "p1"
@test to_s_alp(p1)[1:9] == "plane  p1"
@test p1.id == 1
@test p1.kind == 1
@test p1.at == 55
@test p1.lb == 130
@test p1.target == 156
@test p1.ub == 560
@test p1.ep == 10.0
@test p1.tp == 10.0

@test repr(inst.sep_mat) == "[3 15; 15 8]"
@test sum(inst.sep_mat) == 41
@test p2.kind == 1
@test p3.kind == 2
@test get_sep(inst, p2, p3) == 15

lg1(" fait\n")

#./
