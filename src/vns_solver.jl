export VnsSolver

"""
VnsSolver

Résoud le problème global par méthode à voisinage variable.

"""
mutable struct VnsSolver
    inst::Instance
    VnsSolver() = new() # Constructeur par défaut (sans paramètre)
end

function VnsSolver(inst::Instance; startsol::Union{Nothing,Solution} = nothing)
    ln3("Début constructeur de VnsSolver")

    this = VnsSolver()
    this.inst = inst

    println("VnsSolver non implanté...")
    println("   AU BOULOT ! ")
    exit(1)

    return this
end

# END TYPE VnsSolver
