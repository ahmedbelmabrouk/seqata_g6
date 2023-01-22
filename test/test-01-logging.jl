# ===========
# ATTENTION : CE TEST SUPPOSE QU'IL EST LUI-MEME FONCTIONNEL  (lg1(), ...)
lg1("Test des logging : pile des levels... ")
@testset "Test des logging : pile des levels" begin
    Log.level(1)
    @test Log.level() == 1

    Log.level(2)
    @test Log.level() == 2

    @test_throws ErrorException Log.poplevel!()

    @test Log.pushlevel!(3) == 3
    @test Log.level() == 3
    @test Log.pushlevel!(4) == 4
    @test Log.level() == 4

    @test Log.poplevel!() == 3
    @test Log.level() == 3

    @test Log.poplevel!() == 2
    @test Log.level() == 2

    @test_throws ErrorException Log.poplevel!()
    @test Log.level() == 2

end
ln1("fait.")

##########################
# CAPTURE DE STDOUT POUR RÉCUPÉRATION ULTÉRIEURE DANS UNE VARIABLE 

Log.pushlevel!(2)
original_stdout = stdout
(read_pipe, write_pipe) = redirect_stdout();

# --------------------------------------------
lg2("A") # Affiche A sans saut de ligne
lg2("B") # Affiche B sans saut de ligne
ln2("")  # Affiche un saut de ligne

# --------------------------------------------
ln2("C") # Affiche C avec saut de ligne

# --------------------------------------------
# Affiche des valeura de type varié sans saut de ligne malgré le ln3()
# car suffix est positionné pour annuler le "\n" par défaut.
n = 5
f = 3.1
r = 22 // 7
ln2("n=", n, " f=", f, " r=", r, prefix = "UN_PREFIX==", suffix = "==UN_SUFFIX")

# --------------------------------------------
lg2("-1-") # On ajoute à la fin de la ligne précédente
ln2() # n'affiche rien, même pas un saut de ligne
lg2("-2-") # On ajoute à la fin de la ligne précédente
ln2("") # affiche un saut de ligne

# --------------------------------------------
ln2("FIN") # affiche "FIN" puis un saut de ligne

##########################
# RESTAURATION DE STDOUT
redirect_stdout(original_stdout)
close(write_pipe)
Log.poplevel!()

##########################
# VÉRIFICATION DE L'ECRITURE CAPTURÉE 
lg1("Test des logging : capture de stdout... ")
@testset "Test des logging via capture de stdout par redirect_stdout" begin

    # --------------------------------------------
    # txt = read(read_pipe, String)
    line = readline(read_pipe, keep = true)
    # @show line
    @test line == "AB\n"

    # --------------------------------------------
    line = readline(read_pipe, keep = true)
    # @show line
    @test line == "C\n"

    # --------------------------------------------
    line = readline(read_pipe, keep = true)
    # @show line
    # @test line == "UN_PREFIX==n=5 f=3.1 r=22//7==UN_SUFFIX" # NON car suffix écrasé !
    @test line == "UN_PREFIX==n=5 f=3.1 r=22//7\n"

    # --------------------------------------------
    #
    line = readline(read_pipe, keep = true)
    # @show line
    @test line == "-1--2-\n"

    # --------------------------------------------
    line = readline(read_pipe, keep = true)
    # @show line
    @test line == "FIN\n"

end
# println("\nVALEUR ACTUELLE DE Log.level()=", Log.level())
ln1("fait.")
