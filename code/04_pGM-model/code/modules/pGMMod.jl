using Trapz
using Parameters

function pGM!(time, potential, p; pGM::Float64=0.0)
    """
    Calculate the surface excesses, total surface tension & total Gibbs-Marangoni pressure at separation 'h'

    Args:
        time (array):
            1. t (float): Distance from G/L-interface [m]
        potential (array):
            1. Electrostatic potential at distance t [V]
        p (tuple):
                1. constants (tuple):
                    1.1 STconst (float): Surface tension constant [M.m.Å/mN]
                    1.2 Avog (float): Avogadro constant [1/mol]
                    1.3 ionS (float): Ionic strength [M]
                    1.4 beta (float): Thermodynamic beta [1/J]
                    1.5 h (float): Separation distance [nm]
                2. ion_tot (tuple):
                    2.1. q (float): Valency of the ion [-]
                    2.2. ah (float): Hydrated radius of the ion [m]
                    2.3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
                3. n_salts (integer): Number of salts present in the system [-]
        print_flag (boolean): Statement to (not) print the checks [-]

    Returns:
        conc_matrix (array):
            1. conc (float): Ion's local concentration fraction at distance t from G/L-interface [-]
        gammacon_qc (array):
            1. gammacon_qc (float): Ion's surface excess × charge × concentration [Å.M]
        tension (float): Total surface tension [mN/m.M]
        pGM (float): Total Gibbs-Marangoni pressure [Pa]
    """
    @unpack beta, elc, epsilon_o, epsilon_w, Avog, kappa, STconst, ionS, h = p
    @unpack ion_conc, ion_charges, ion_hyd_radii, ion_Wcal, ion_types = p
    
    conc_matrix = zeros(Float64, size(time)[1], length(ion_conc))
    for (i, t) in enumerate(time)
        for (j, (c_i, ch_i, hr_i, W_i, type_i)) in enumerate(zip(ion_conc, ion_charges, ion_hyd_radii, ion_Wcal, ion_types))

            if type_i == "beta"
                if t < 1e10*hr_i
                    U = W_i*1e10*hr_i/t * exp(-2kappa * (1e-10t - hr_i)) - 2.1
                else
                    U = W_i*1e10*hr_i/t * exp(-2kappa * (1e-10t - hr_i))
                end
            elseif type_i == "proton"
                if t < 1e10*hr_i
                    U = beta/(4pi*1e-10t) * (elc^2)/(4epsilon_o*epsilon_w) * exp(-2kappa*1e-10t) - 3.05
                else
                    U = beta/(4pi*1e-10t) * (elc^2)/(4epsilon_o*epsilon_w) * exp(-2kappa*1e-10t)
                end
            else
                if t < 1e10*hr_i
                    U = 1000
                else
                    U = W_i*1e10*hr_i/t * exp(-2kappa * (1e-10t - hr_i))
                end
            end
            conc_matrix[i, j] = exp(-U - ch_i*potential(t)[1]) - 1 # Fractionate ionic concentration profile [M/M] = [-], Eq.6
        end
    end

    gammacon = zeros(1, length(ion_conc))
    gammacon_qc = zeros(1, length(ion_conc))
    
    for (i, (c_i, ch_i)) in enumerate(zip(ion_conc, ion_charges))
        # Surface excess over cocentration [M.Å/M] = [Å], Eq.6
        gammacon[i] = trapz(time, conc_matrix[:, i]) 
        
        # Surface excess × charge × concentration [Å.M]
        gammacon_qc[i] = gammacon[i] * c_i * ch_i
    end
    
    # Surface tension [Å] / [M.m.Å/mN] -> [mN/m.M]
    tension = sum(gammacon) / STconst
    
    for (index, c_i) in enumerate(ion_conc)
        # Surface excess squared [Å^2]
        gammacon_squared = gammacon[index]^2
        
        # Gibbs-Marangoni pressure [m4.mol.kg/m5.mol.s2] -> [Pa], Eq.5
        pGM += 4Avog * c_i * 1000 * gammacon_squared * 1e-20 / (beta * h^2)
    end

    return (conc_matrix, gammacon_qc, tension, pGM)
end
