# cd(@__DIR__)
# using Pkg
# Pkg.activate(".")
using BoundaryValueDiffEq, OrdinaryDiffEq, Printf#, BenchmarkTools
using Parameters
# import BoundaryValueDiffEq: BVProblem

include("welcome.jl")
include("constants.jl")
include("modules/CheckMod.jl")
include("modules/WMod.jl")
include("modules/mPBEMod.jl")
include("modules/pGMMod.jl")

using .Check_Mod

@with_kw struct Para
    beta::Float64
    elc::Float64 = 1.6021765e-19    # Elementary charge [C]
    epsilon_o::Float64 = 8.85418782e-12 # Vacuum's permitivity [F/m]
    epsilon_w::Float64 = 78.3       # Water's permitivity [F/m]
    Avog::Float64 = 6.0221415e23    # Avogadro's number [mol-1]
    h::Float64 = 10e-9  
    kappa::Float64
    # n::Int64
    STconst::Float64
    ionS::Float64
    ion_conc::Vector{Float64} = Vector{Float64}[]
    ion_charges::Vector{Int64}  = Vector{Float64}[]
    ion_hyd_radii::Vector{Float64}  = Vector{Float64}[]
    ion_Wcal::Vector{Float64}  = Vector{Float64}[]
    ion_types::Array{String}  = Vector{String}[]
end

const kB = 1.3806503e-23         # Boltzmann constant [J/K]
const elc = 1.6021765e-19        # Elementary charge  [C]
const Avog = 6.0221415e23        # Avogadro constant  [1/mol]
const epsilon_o = 8.85418782e-12 # Vacuum's permittivity [F/m]
const epsilon_w = 78.3           # Water's permittivity [F/m]

function main(list=false; 
    print_flag::Bool=true, T::Float64=297.15, h::Float64=10e-9, ionS::Float64=0.0, ion_tot::Tuple=(),
    abstol::Float64=1e-9, reltol::Float64=1e-9, fun=mPBE!, bc=bc!)
    """
    Solve the modified Poisson–Boltzmann equation (mPBE) and compute the Gibbs-Marangoni pressure for a given electrolyte solution.
    
    # Args:
        list (array or false): 
            Array of arrays, each specifying an ion as 
            [concentration [M], charge, hydrated radius [m], type (optional)]. 
            If type is not provided, alpha ion is assumed.
            If false, a default list is used.
        print_flag (Bool): 
            Statement to (not) print the checks and results.
        T (Float64): 
            Temperature [K].
        h (Float64): 
            Separation distance [m].
        ionS (Float64): 
            Initial ionic strength [M].
        ion_tot (Tuple): 
            Additional ion parameters (unused by default).
        abstol (Float64): 
            Absolute tolerance for the BVP solver.
        reltol (Float64): 
            Relative tolerance for the BVP solver.
        fun (Function): 
            Function for the mPBE ODE system (default: mPBE!).
        bc (Function): 
            Boundary condition function (default: bc!).
    
    # Returns:
        pGM (Float64): 
            Gibbs-Marangoni pressure [Pa].
        surface_tension_gradient (Float64): 
            Calculated surface tension gradient [mN/m.M].
        retcode: 
            Return code from the BVP solver indicating success or failure.
    
    # Example:
        pGM, surface_tension, retcode = main()
    """
    

    if list == false
        # contains molar concentration, charge, hydrated radius, type {'alpha', 'beta', 'proton'}
        # list = [[0.1, +1, 2.5e-10], [0.1, -1, 2.0e-10]]
        
        # Extra tough concentration of NH4Cl from other work
        # list = [[0.0053, +1, 2.5e-10], [0.0053, -1, 2.0e-10]]
        
        # Extra tough concentration of NaAc from other work
        list = [[0.95, +1, 2.5e-10, "alpha"], [0.95, -1, 3.22e-10, "beta"]]
    end
    # n = length(list)/2 # Number of salts [-]

    beta::Float64 = 1 / (kB * T)            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
    STconst = -1e-6 * 1e10 * beta / Avog    # Surface tension constant [mol/N.m] -> [M.m.Å/mN]

    for ion in list
        ionS += ion[1] * ion[2]^2
    end
    ionS = 0.5*ionS # Ionic strength [M]

    kappa::Float64 = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)) # Inverse Debye-Hückel lenght [1/m]
    if kappa == 0.0
        return [0, 0, 0, 0]
    end
    boundary = 10e10/kappa         # 10x Debye-Hückel lenght [Å]

    params = Para(
        beta = beta,
        kappa = kappa,
        STconst = STconst,
        ionS = ionS
    )

    n_present = 0
    for (i, ion) in enumerate(list)
        if ion[1] > 0.0
            n_present += 1
        end
    end
    ion_conc = zeros(Float64, n_present)
    ion_charges = zeros(Int64, n_present)
    ion_hyd_radii = zeros(Float64, n_present)
    ion_Wcal = zeros(Float64, n_present)
    ion_types = Array{String}(undef, n_present)

    i = 1
    for (_, ion) in enumerate(list)
        if ion[1] == 0.0
            continue
        end
        Wcal = W(ion[2], ion[3], params)
        ion_conc[i] = ion[1]
        ion_charges[i] = ion[2]
        ion_hyd_radii[i] = ion[3]
        ion_Wcal[i] = Wcal
        if last(ion) isa String
            ion_types[i] = ion[4]
        else
            ion_types[i] = "alpha"
        end
        i += 1
    end
    # println(ion_types)
    
    # --- Electrostatic potential ---
    # mPBE_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa)
    # mPBE_param = (mPBE_constants, ion_tot, n)
    mPBE_param = Para(params,
        ion_conc = ion_conc,
        ion_charges = ion_charges,
        ion_hyd_radii = ion_hyd_radii,
        ion_Wcal = ion_Wcal,
        ion_types = ion_types)
    # println(eltype(mPBE_param))

    paramized_fun(du, u, p, t) = fun(du, u, mPBE_param, t)
    
    tspan = (boundary, 0.0)
    u0  = [0.0, 0.0]

    bvp_problem = BVProblem(paramized_fun, bc, u0, tspan)
    bvp_sol     = solve(
        bvp_problem,
        # Rodas4P worked for the 'tough' NH4Cl concentration at tolerances of 1e-9
        MultipleShooting(;ode_alg = Rodas4P(autodiff = false), nshoots = 20),   
        # MultipleShooting(;ode_alg = TRBDF2(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = Rodas5P(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = ROS3(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = DFBDF(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = FBDF(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = QNDF(autodiff = false), nshoots = 20),
        # Shooting(;ode_alg = AutoVern7(Rodas4(autodiff = false))),
        abstol = abstol,
        reltol = reltol)

    # println(bvp_sol[1], bvp_sol[-1])
    # --- Gibbs-Marangoni pressure ---
    time_range = range(0.0, boundary, 100001)
    # pGM_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa, STconst, ionS, h)
    # pGM_param = (pGM_constants, ion_tot, n)
    pGM_sol = pGM!(time_range, bvp_sol, mPBE_param)

    # println("\nChecks:")
    Check_Mod.check_input(list, print_flag=print_flag)  
    Check_Mod.check_mPBE(bvp_sol, print_flag=print_flag)
    Check_Mod.check_pGM(pGM_sol[2], print_flag=print_flag)
    
    if print_flag == true
        println("\nResults:")
        println("   • Surface excess × q × c   = ", round.(pGM_sol[2]; digits=3), " = ", round.(sum(pGM_sol[2]); digits=3), " Å")
        @printf("   • Surface tension gradient = %.3e mN/m.M", pGM_sol[3])
        @printf("\n   • Electrostatic potential  = %.3f mV", 1000*last(bvp_sol.u)[1])
        @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa\n", pGM_sol[4])
    end

    return pGM_sol[4], pGM_sol[3], SciMLBase.successful_retcode(bvp_sol)
end
