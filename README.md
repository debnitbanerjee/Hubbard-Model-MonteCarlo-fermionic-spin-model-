# Hubbard-Model-MonteCarlo-fermionic-spin-model-
Julia implementation of a finite-temperature Monte Carlo simulation for fermions coupled to auxiliary classical magnetic fields on a 2D square lattice. To study antiferromagnetic phase transitions, k- structure factors, and critical temperature (Tc) behaviors.


1.Constructs the single-particle Hamiltonian matrix for a given configuration of classical spins and exactly diagonalizes it to compute the fermionic free energy.
2. Uses the Metropolis Hastings algorithm to update the classical magnetic fields, simulating the system's thermal equilibration at a given temperature.
3. Computes the magnetic structure factor S(q) to detect spatial ordering, specifically tracking (π,π) odering.
4. Sweeps through temperatures to find Tc by tracking the inflection points of the structure factor
