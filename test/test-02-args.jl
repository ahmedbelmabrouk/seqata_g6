# ===========
lg1("Test de l'analyse des arguments passées en ligne de commande... ")

# Si le premier paramètre est inconnu (action "xxx"), on doit lever une erreur
@test_throws ArgumentError args = Args.parse_commandline("xxx")

# Si le fichier d'instance n'existe pas, on doit lever une erreur
@test_throws ArgumentError args = Args.parse_commandline("test --infile fichier_bidon.alp")

# le Log.level va être cassé par le test suivant
oldlevel = Log.level() # NE PAS utiliser Log.pushlevel(xx) ici !
# Test de l'analyse de la chaine passée en arguments
cli = "
    explore
    --infile $APPDIR/test/data/alp_01.alpx
    --loglevel 2
    --presort shuffle
    -L0
    -n100
    -d/tmp
"
args = Args.parse_commandline(cli)
Log.level(oldlevel) # restauration pour les tests suivants

# On vérifier que to_s_dict(args) retourne bien le dict des arguments
# Args.show_dict(args)
# println(Args.to_s_dict(args))
@test occursin(r"\spresort\s*=>\s*shuffle\s", Args.to_s_dict(args))

# Test d'accès aux arguments via le dict
@test args[:infile] == "$APPDIR/test/data/alp_01.alpx"
@test args[:loglevel] == 0

# Test d'accès aux arguments via l'accesseur
@test Args.get(:loglevel) == args[:loglevel]
@test Args.get("loglevel") == args[:loglevel]

# Test de la modification d'un arguments
Args.set(:loglevel, 3)
@test args[:loglevel] == 3
@test Args.get(:loglevel) == args[:loglevel]
@test Args.get("loglevel") == args[:loglevel]

lg1("fait.\n")
