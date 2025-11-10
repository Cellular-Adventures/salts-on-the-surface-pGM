module Check_Mod
using OrdinaryDiffEq

function check_input(list; charge_tot::Float64=0.0, print_flag::Bool=false)
    """Check if sum of charges = 0 and number of salts is correct.

    Args:
        list (tuple):
            1. q (float): Valency of the ion [-]
            2. ah (float): Hydrated radius of the ion [m]
            3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
        print_flag (boolean): Statement to (not) print the checks [-]

    Returns:
        nothing
    """
    for ion in list
        charge_tot += ion[1]*ion[2] # Sum of charges [-]
    end

    if charge_tot != 0
        println("   ✖ WARNING: Sum of charges is NOT 0!")
    else
        if print_flag
            println("   ✔ Sum of charges is 0.")
        end
    end
    nothing
end

function check_mPBE(sol; print_flag::Bool=false)
    if SciMLBase.successful_retcode(sol)
        if print_flag
            println("   ✔ BVP succesfully solved.")
        end
    else
        println("   ✖ WARNING:\tBVP NOT succesfully solved!")
        println("           \t\tReturn code: $(sol.retcode)")
    end
    nothing
end

function check_pGM(sol; print_flag::Bool=false)
    if sum(sol) > 0.05
        println("   ✖ WARNING: Electroneutrality condition is NOT satisfied!")
    else
        if print_flag
            println("   ✔ Electroneutrality condition is satisfied.")
        end
    end
    nothing
end

end # end of module
