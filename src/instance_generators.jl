
export instance_build_mini10

# Génère un certain nombre d'instances de test
#
# TODO : MODULE LE CONSTRUCTEUR DE TYPE Instance POUR LIRE UNE CHAINE
# DANS UN DES FORMATS ACCEPTÉ
#
# La description d'un avion est de la forme :
#
#   #    name kind  at   E    T   L    ep    tp
#   plane  p1    1  55  130  156 560  10.0  10.0
#
# ou bien
#   #    name kind  at  E     T   L    dt1  c1   dt2 c2    dt3 c3
#   plane  p1    1  55  130  156 560  -26 260.0   0 0.0    404 4040.0
#
function instance_build_mini10()
    inst = Instance()
    inst.name = "mini10"
    inst.nb_kinds = 2
    inst.nb_planes = 10
    inst.freeze_time = 10

    inst.planes = Vector{Plane}()

    # inst.sep_mat = Matrix{Int}(undef, 2, 2) # ok
    inst.sep_mat = zeros(Int, 2, 2) # ok
    inst.sep_mat[1, 1] = 3
    inst.sep_mat[1, 2] = 15
    inst.sep_mat[2, 1] = 15
    inst.sep_mat[2, 2] = 8

    add_plane(inst, "  p1    1    55   130   156   560    10.0 10.0")
    add_plane(inst, "  p2    1   121   196   259   745    10.0 10.0")
    add_plane(inst, "  p3    2    15    90    99   511    30.0 30.0")
    add_plane(inst, "  p4    2    22    97   107   522    30.0 30.0")
    add_plane(inst, "  p5    2    36   111   124   556    30.0 30.0")
    add_plane(inst, "  p6    2    46   121   136   577    30.0 30.0")
    add_plane(inst, "  p7    2    50   125   139   578    30.0 30.0")
    add_plane(inst, "  p8    2    52   127   141   574    30.0 30.0")
    add_plane(inst, "  p9    2    61   136   151   592    30.0 30.0")
    add_plane(inst, " p10    2    86   161   181   658    30.0 30.0")


    return inst
end
