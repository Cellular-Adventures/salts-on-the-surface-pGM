# Gibbs-Marangoni pressure
Numerical solver of the **modified Poisson-Boltzmann equation** from Duignan (2021) in Julia. This equation is modified as it not only includes an electrostatic potential but also an adsorption energy function of the present ions [1].

<img src="https://github.com/stijnrommens/GM_pressure/assets/113170925/cc816c5d-426a-4f06-9b4b-0b75e1747da0" width="400">

The model needs the following input:
- ionic concentrations;
- hydrated radius and charge of each ion;

As output, the electrostatic potential of salts at the interface is obtained. This can be used to calculate the surface excess of ions, which results in the **Gibbs-Marangoni pressure** (at 10 nm separation).

Sources:\
[1] Duignan, T. T. (2021). The surface potential explains ion-specific bubble coalescence inhibition. _Journal of Colloid and Interface Science, 600_, 338–343. https://doi.org/10.1016/j.jcis.2021.04.144

# Requirements
Julia (v1.9) needs to be installed. The attached _.tomls_ in this repo should automatically install the required packages.

# How to use
Open a Julia REPL in the base folder (`~\GM_pressure>julia`), and activate the environment (`julia>] activate .`).

In the Julia REPL, include the file `main.jl` (`julia>include("code/main.jl")`). This gives you access to the `main()` function which is used to calculate the GM pressure for an electrolyte mixture.

The first and most important input to `main()` is the list of ions. In this list, each input is expected to be a list describing an individual ion as `[concentration, charge, hydrated radius, type]`.

The function returns three values: a `float` for the GM pressure in Pa, a `float` for the surface tension gradient ($\delta\gamma/\delta c$) in $\mathrm{\left[mN* m^{-1}* M^{-1}\right]}$ and a `bool` success flag indicating if the system solved successfully.

## Example 1:
Calculating the GM pressure for an 0.5 M NaCl solution:

```julia
main([
  [0.5, +1, 2.5e-10, "alpha"],
  [0.5, -1, 2.0e-10, "alpha"]
])
```

which will return

```julia
(8759.22..., 1.471..., true)
```

## Example 2:
Calculating the GM pressure for an 0.95 M NaCH3COO solution:

```julia
main([
  [0.95, +1, 2.5e-10, "alpha"],
  [0.95, -1, 3.22e-10, "beta"]
])
```

which will return

```julia
(6048.09..., 0.886..., true)
```

In this case, the solver will throw some warnings during the solution finding. As long as the final return flag is true, these warnings were overcome later during the solving of the system of equations and can safely be ignored.

# Debugging & Tweaking
Certain scenarios (e.g. high ionic strengths) may result in warnings/errors from the BVP-solver. This will most likely return values which go to infinity.
`main()` accepts several keyword arguments that can help increase stability (at the cost of speed) or increase speed (at the cost of stability):

- absolute tolerance (**abstol**); currently at 1e-9 (default is 1e-6)
- relative tolerance (**reltol**); currently at 1e-9 (default is 1e-3)

Additionally, the algorithm used can be changed around line 166. 

Accuracy/validity are tested by checking if the surface excess of all ions fulfils the electroneutrality condition by checking that $\sum^n \Gamma_i * q_i * c_i = 0$.
