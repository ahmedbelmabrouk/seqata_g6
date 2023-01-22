@ms include("mip_discret_solver.jl")

export main_dmip

function main_dmip()
    println("="^70)
    println("Début de l'action dmip")
    inst = Instance(Args.get(:infile))
    sv = MipDiscretSolver(inst)
    solve!(sv)

    bestsol = sv.bestsol
    print_sol(bestsol)
    print("Création du fichier \"$(guess_solname(bestsol))\"... ")
    write(bestsol)
    println("FAIT !")

    println("Fin de l'action dmip")
end
