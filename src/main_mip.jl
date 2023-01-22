@ms include("mip_solver.jl")

export main_mip

function main_mip()
    println("="^70)
    println("Début de l'action mip")
    inst = Instance(Args.get(:infile))
    sv = MipSolver(inst)
    solve!(sv)

    bestsol = sv.bestsol
    print_sol(bestsol)
    print("Création du fichier \"$(guess_solname(bestsol))\"... ")
    write(bestsol)
    println("FAIT !")

    println("Fin de l'action mip")
end
