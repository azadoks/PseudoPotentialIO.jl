# UPF.jl

Support reading UPF (Unified Pseudopotential Format) versions 1 and 2.

At the moment, it is almost a direct translation of Simon Pintarelli's Python library [`upf_to_json`](https://github.com/simonpintarelli/upf_to_json).

## `upf_to_json` data structure
```
pseudo_potential
    header
        element::String
        pseudo_type::String
        core_correction::Bool
        spin_orbit::Bool
        paw_core_energy::Float64
        z_valence::Float64
        l_max::Int
        cutoff_radius_index::Int
        mesh_size::Int
        number_of_wfc::Int
        number_of_proj::Int
    radial_grid::Vector{Float64}
    core_charge_density::Vector{Float64}
    local_potential::Vector{Float64}
    beta_projectors::Vector{
        label::String
        angular_momentum::Int
        radial_function::Vector{Float64}
    }
    D_ion::Vector{Float64}
    augmentation::Vector{
        i::Int
        j::Int
        angular_momentum::Int
        radial_function::Vector{Float64}
    }
    atomic_wave_functions::Vector{
        label::String,
        angular_momentum::Int,
        occupation::Float64,
        radial_function::Vector{Float64}
    }
    paw_data
        aug_integrals::Vector{Float64}
        aug_multipoles::Vector{Float64}
        ae_wfc::Vector{
            radial_function::Vector{Float64}
            angular_momentum:Int
        }
        ps_wfc::Vector{
            radial_function::Vector{Float64}
            angular_momentum::Int
        }
        occupations::Vector{Float64}
        ae_core_charge_density::Vector{Float64}
        ae_local_potential::Vector{Float64}
    total_charge_density::Vector{Float64}
```

## `libpspio` data structure (norm-conserving only)
- [x] Norm-conserving (local+non-local)
- [ ] Norm-conserving (semilocal)
- [x] Multiple projectors per angular momentum
- [x] Atomic wavefunctions
- [x] Spin-orbit coupling
- [ ] Ultrasoft
- [ ] PAW

To support ultrasoft, add a block for augmentation.
To support paw, add agumentation, paw, and gipaw blocks.
```
struct PSPData
    # General data
    psp_info::PSPInfo,   # Metadata
        author::String,       # Author
        code_name::String,    # Generator name
        code_version::String, # Generator version
        code_input::String,   # Generator input
        description::String,  # Description
        time::Datetime        # Generation date
    format::Int,         # File format
    symbol::String,      # Atomic symbol
    z::Float64,          # Atomic number
    z_valence::Float64,  # Pseudo-ion charge
    ne_valence::Float64, # Number of electrons (usually == zvalence, except for ions)
    l_max::Int,          # Maximum angular momentum
    wave_eq::Int,        # Type of equation which was solved: Dirac | Scalar relativistic | Schroedinger
    total_energy::Float, # Total energy of pseudo-ion

    # Mesh
    mesh::PSPMesh, # Radial mesh on which _all_ functions are discretized
        type::Int,           # Mesh type
        a::Float64,          # Mesh parameter
        b::Float64,          # Mesh parameter
        np::Int,             # Number of points in mesh
        r::Vector{Float64},  # Mesh points
        rab::Vector{Float64} # Factor required for discrete integration: rab[i] = (dr(x)/dx)_{x=i}

    # States
    qn_to_istate::Matrix{Int}, # Lookup table for quantum number -> state index
    n_states::Int,             # Number of states / wavefunctions
    states::Vector{PSPState},  # States / wavefunctions
        qn::QuantumNumber,  # Quantum numbers (n, l, j) for this state / wavefunction
            n::Int,    # Principal quantum number
            l::Int,    # Angular momentum
            j::Float64 # Total angular momentum
        occ::Float64,       # Occupation of the electronic state
        eigenval::Float64,  # Eigenvalue of the electronic state
        label::String,      # Label describing the state (e.g. "2S" or "4D1.5")
        rc::Float64,        # Cutoff radius used for pseudo-state generation
        wf::PSPMeshFunction # Wavefunction
            mesh::PSPMesh      # In `libpspio`: a pointer to the mesh struct
            f::Vector{Float64} # Function values and derivatives on the mesh

    # Pseudo-potentials
    scheme::Int,                      # Pseudo-potential generation scheme
    n_potentials::Int,                # Number of pseudo-potentials
    potentials::Vector{PSPPotential}, # Pseudo-potentials
        qn::QuantumNumber, # Quantum numbers for the potential
        v::PSPMeshFunction # Pseudo-potential on a radial mesh

    # Non-local projectors
    n_projectors::Int,                   # Number of non-local projectors
    projector_energies::Vector{Float64}, # Dij terms for interactions between projectors
    projectors::Vector{PSPProjector},    # Non-local projectors
        qn::QuantumNumber,    # Quantum numbers for the non-local projector
        energy::Float64,      # Non-local projector energy
        proj::PSPMeshFunction # Non-local projector on a radial mesh
    projectors_l_max::Int,               # Maximum angular momentum of projectors
    l_local::Int,                        # Angular momentum channel used to generate the local potential
    v_local::Vector{PSPPotential},       # Local potential for non-local pseudo-potential form

    # Exchange-correlation data (incl. non-linear core corrections)
    xc::PSPExchangeCorrelation, # Exchange-correlation
        exchange::Int,                # Exchange functional in libxc convention
        correlation::Int,             # Correlation functional in libxc convention
        nlcc_scheme::Int,             # Scheme used for obtaining non-linear core correction
        nlcc_pf_scale::Float64,       # Prefactor for model scale (multiplies core-valence crossover radius)
        nlcc_pf_value::Float64,       # Prefactor for model amplitude (multiplies core-valence crossover value)
        nlcc_density::PSPMeshFunction # Core charge density and its derivatives on a radial mesh

    # Valence charge density
    rho_valence::PSPMeshFunction  # Valence charge density
end
```

## `quantumESPRESSO` data structure
```
struct PseudoUPF
    # General data
    generated::String  # Generator name
    author::String     # Author
    date::String       # Generation date
    comment::String    # Description
    psd::String        # Atomic symbol
    type::String       # Pseudo type (NC | SL | US / USPP | PAW | 1/r)
    rel::String        # Type of relativistic treatment (no | scalar | full)
    tvanp::Bool        # Is ultrasoft
    tcoulombp::Bool    # Is Coulomb
    nlcc::Bool         # Has non-linear core corrections
    is_gth::Bool       # Is Goedecker-Teter-Hutter
    is_multiproj::Bool # Has multiple non-local projectors per angular momentum (NC only, assumed true with US & PAW)
    dft::String        # Exchange-correlation type
    zp::Float64        # Pseudo-ion charge
    etotps::Float64    # Total energy of pseudo-ion
    ecutwfc::Float64   # Suggested wavefunction energy cutoff
    ecutrho::Float64   # Suggested charge density energy cutoff

    # Non-local projectors and atomic wavefunctions
    nv::String                  # UPF version string (e.g. 2.0.1)
    lmax::Int                   # Maximum angular momentum in non-local projectors
    lmax_rho::Int               # Maximum angular momentum in charge density (should be 2 * lmax)
    vnl::Matrix{Float64,3}      # Semi-local potential for single-projector norm-conserving pseudo-potentials
    nwfc::Int                   # Number of atomic wavefunctions
    nbeta::Int                  # Number of non-local projectors
    kbeta::Vector{Int}          # Cutoff radius index for each non-local projector
    kkbeta::Int                 # Maximum cutoff radius index for all non-local projectors
    lll::Vector{Int}            # Angular momentum of each non-local projector
    beta::Matrix{Float64}       # Non-local projectors
    els::String                 # Label of each atomic wavefunction
    els_beta::String            # Label of each non-local projector
    nchi::Vector{Int}           # Principal quantum number of each atomic wavefunciton
    lchi::Vector{Int}           # Angular momentum of each atomic wavefunction
    oc::Vector{Float64}         # Occupation of each atomic wavefunction
    epseu::Vector{Float64}      # One-particle pseudo-energy for each atomic wavefunction
    rcut_chi::Vector{Float64}   # Inner cutoff radius for each atomic wavefunction
    rcutus_chi::Vector{Float64} # Outer cutoff radius for each atomic wavefunction (ultrasoft)
    chi::Matrix{Float64}        # Atomic wavefunctions (used for initial wavefunction)
    rho_at::Vector{Float64}     # Atomic charge-density (used for initial charge density)
    mesh::Int                   # Number of points in mesh
    xmin::Float64               # Minimum `x` of the linear mesh
    rmax::Float64               # Maximum radius of the mesh
    zmesh::Float64              # Nuclear charge used for the mesh
    dx::Float64                 # Δx used for the linear mesh
    r::Vector{Float64}          # Radial grid
    rab::Vector{Float64}        # Factor required for discrete integration: rab[i] = (dr(x)/dx)_{x=i}
    rho_atc::Vector{Float64}    # Atomic core charge

    # Local potential
    lloc::Int             # Angular momentum channel used to generate the local potential
    rcloc::Float64        # Radius above which the local potential is equal to the all-electron potential
    vloc::Vector{Float64} # Local atomic potential
    dion::Matrix{Float64} # Dij terms for interactions between projectors

    # Augmentation
    q_with_l::Bool           # `qfunc` is pseudized separately for each angular momentum
    nqf::Int                 # Number of Q coefficients
    nqlc::Int                # Number of angular momenta in Q
    qqq_eps::Float64         # `qfunc` is null if its norm is < `qqq_eps`
    rinner::Vector{Float64}  # Inner cutoff radius `r_L`
    qqq::Matrix{Float64}     # Qij
    qfunc::Matrix{Float64}   # `Qij(|r|)` for `|r| > r_L` without angular momentum dependence
    qfuncl::Array{Float64,3} # `Qij(|r|)` for `|r| > r_L` with angular momentum dependence (optional, except for PAW)
    qfcoef::Array{Float64,4} # Coefficients for `Q` for `|r| < r_L`

    # Atomic wavefunctions
    has_wfc::Bool          # Has all-electron and pseudo- wavefunctions (different from chi; one for each projector)
    aewfc::Matrix{Float64} # All-electron wavefunctions for each non-local projector
    pswfc::Matrix{Float64} # Pseudo-wavefunctions for each non-local projector

    # Spin-orbit coupling
    has_so::Bool            # Has spin-orbit coupling
    nn::Vector{Float64}     # Principal quantum number for each atomic wavefunction
    rcut::Vector{Float64}   # Inner? cutoff radius for each non-local projector
    rcutus::Vector{Float64} # Outer? cutoff radius for each non-local projector (ultrasoft)
    jchi::Vector{Float64}   # Total angular momentum for each atomic wavefunction (j=l+1/2 or j=l-1/2)
    jjj::Vector{Float64}    # Total angular momentum for each non-local projector (j=l+1/2 or j=l-1/2)

    # PAW
    paw_data_format::Int # Version of the PAW data format
    tpawp::Bool          # Is PAW
    paw::PAWData
        ae_rho_atc::Vector{Float64} # All-electron core charge
        pfunc::Array{Float64,3}     # Psi_i(r) * Psi_j(r)
        pfunc_rel::Array{Float64,3} # Psi_i(r) * Psi_j(r) small component
        ptfunc::Array{Float64,3}    # "as `pfunc`, but for pseudo
        aewfc_resl::Matrix{Float64} # "as `pfunc_rel`, but for pseudo"
        ae_vloc::Vector{Float64}    # All-electron local potential
        oc::Vector{Float64}         # Starting occupation used in initializing becsum for each non-local projector (US is for each wfc)
        augmom::Array{Float64,3}    # Multipole all-electron "pseudo"
        raug::Float64               # Cutoff radius for "augfunction"
        iraug::Int                  # Index on radial grid close to and greater than `raug`
        lmax_aug::Int               # Maximum angular momentum of augmentation functions
        core_energy::Float64        # Constant shift to recover all-electron energy
        augshape::String            # Shape of augmentation charge
    
    # GIPAW
    has_gipaw::Bool                       # Has GIPAW data
    paw_as_gipaw::Bool                    # "EMINE"
    gipaw_data_format::Int                # Version of the GIPAW data format
    gipaw_ncore_orbitals::Int             # ?
    gipaw_core_orbital_n::Vector{Float64} # ?
    gipaw_core_orbital_l::Vector{Float64} # ?
    gipaw_wfs_nchannels::Int              # ?
    gipaw_wfs_el::Vector{String}          # ?
    gipaw_wfs_ll::Vector{Int}             # ?
    gipaw_wfs_ae::Matrix{Float64}         # ?
    gipaw_wfs_rcut::Vector{Float64}       # ?
    gipaw_wfs_rcutus::Vector{Float64}     # ?
    gipaw_wfs_ps::Matrix{Float64}         # ?
end
```