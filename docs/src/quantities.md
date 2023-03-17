# Quantity definitions
A more thorough discussion of pseudopotentials, from which many of the equations below have been adapted can be found [here](https://www.tcm.phy.cam.ac.uk/~jry20/gipaw/tutorial_pp.pdf).

## Norm-conserving pseudopotentials
Following the procedure of Vanderbilt 1991, separable norm-conserving pseudopotentials can be directly produced using the following procedure:

- solve the all-electron Kohn-Sham equations for the isolated atom, yielding the all-electron potential ``V(r)`` and atomic wave-functions ``|\phi_{l,n}\rangle`` with energies ``\varepsilon_{l,n}``
- generate a _local potential_ ``V_{\mathrm{loc}}(r)`` such that ``V_{\mathrm{loc}}(r) = V(r)`` for ``r > r_L`` (``r_L`` is the inner pseudization radius for the local potential); ``V_\mathrm{loc}(r)`` for ``r < r_L`` can be any smooth regular function
- generate pseudo-atomic wavefunctions ``|\tilde{\phi}_{l,n}\rangle`` such that ``\tilde{\phi}_{l,n}(r) = \phi_{l,n}(r)`` for ``r > r_{c,l,n}`` (``r_{c,l,n}`` is the inner cutoff radius for the n-th pseudo-atomic wavefunction at angular momentum ``l``); ``\tilde{\phi}_{l,n}(r)`` for ``r < r_{c,l,n}`` can be any smooth regular function
- generate the corresponding functions ``|\chi_{l,n}\rangle`` (vanishing for ``r > r_{c,l,n}``):
```math
|\chi_{l,n}\rangle \equiv (\varepsilon_{l,n} - T - V_{\mathrm{loc}})|\tilde{\phi}_{l,n}\rangle
```
- generate the KB projectors ``|\beta_{l,m}\rangle``:
```math
|\beta_{l,m}\rangle \equiv \sum_{m} (B^{-1})_{l,nm} |\chi_{l,m}\rangle
```
where ``B_{l,nm} = \langle \tilde{\phi_{l,n}} | \chi_{l,n} \rangle`` and ``|\beta_{l,m}\rangle`` satisfy ``\langle \beta_{l,n} | \tilde{\phi}_{l,m} \rangle = \delta_{nm}``

## Form of the pseudopotentials
PseudoPotentialIO assumes a the above separable (Kleinman-Bylander) form for all the pseudopotentials.
Therefore, the total pseudopotential is defined as:
```math
\hat{V}^{\mathrm{PsP}} \rightarrow \hat{V}_{\mathrm{KB}} = \hat{V}^\prime_{\mathrm{loc}} + \hat{V}_{\mathrm{NL}}
```

The non-local part of the potential ``\hat{V}_{l,\mathrm{NL}}`` at angular momentum ``l`` is defined as
```math
\hat{V}_{l,\mathrm{NL}} \equiv \sum_{nm} | \beta_{l,n} \rangle D_{l,nm} \langle \beta_{l,m} |
```
where ``\beta_{l,n}`` is the n-th non-local projector at angular momentum ``l`` in Kleinman-Bylander (KB) form, and ``D_{l,nm}`` are the KB energies or projector coupling coefficients at angular momentum ``l``.

## Storing the quantities
- The local part of the potential ``\hat{V}^\prime_{\mathrm{loc}}`` is stored in numerical pseudopotentials as the vector `Vloc`, without any prefactor, i.e. the stored quantity is
```math
\hat{V}_{\mathrm{loc}}(r)
```
- The KB projectors ``\beta_{l,n}`` are stored in numerical pseudopotentials as `β[l][n]` with a prefactor of ``r^2``, i.e. the stored quantity is
```math
r^2 \beta_{l,n}(r)
```
- The KB energies / projector coupling coefficients are stored in all pseudopotentials as `D[l][n,m]`.
- If available, the ``\chi_{l,n}`` functions are stored in numerical pseudopotentials as `χ[l][n]` with a prefactor of ``r^2``, i.e. the stored quantity is
```math
r^2 \chi_{l,n}(r)
```
- If available, the pseudo-atomic valence charge density ``\rho_{\mathrm{val}}(r) = \sum_{l=0}^{l_\mathrm{max}} \sum_{m=-l}^{l} \sum_{n} |\tilde{\phi}_{l,n}(r)|^2`` is stored with a prefactor of ``r^2``, i.e. the stored quantity is
```math
r^2 \rho_{\mathrm{val}}(r)
```
- If available, the core charge density (non-linear core correction) ``\rho_{\mathrm{core}}(r)`` is stored with a prefactor of ``r^2``, i.e. the stored quantity is
```math
r^2 \rho_{\mathrm{core}}(r)
```
