function random_versor()
    r, i, j, k = randn(4)
    quat_inorm = inv(sqrt(r * r + i * i + j * j + k * k))
    return quat(r * quat_inorm, i * quat_inorm, j * quat_inorm, k * quat_inorm)
end

function rotate_vector(q::Quaternion, u::AbstractVector)
    q_u = Quaternion(0, u[1], u[2], u[3])
    q_v = q * q_u * conj(q)
    return [imag_part(q_v)...]
end

function rescale_vector(R, rmin, rmax)
    isinf(rmin) && return R
    isinf(rmax) && return R
    return R / norm(R) * (rand() * (rmax - rmin) + rmin)
end

QUANTITIES = [ValenceChargeDensity(), CoreChargeDensity(), BetaProjector(),
              ChiProjector(), BetaCoupling(), LocalPotential(),
              AugmentationFunction()]

@testset "AbstractPsP" begin
    @testset "$(splitpath(filepath)[end])" for filepath in TEST_FILEPATHS
        file = load_psp_file(filepath)
        if (is_norm_conserving(file) | is_ultrasoft(file)) & !has_spin_orbit(file)
            psp = load_psp(file)
            @test isa(identifier(file), AbstractString)
            @test isa(element(psp), PeriodicTable.Element)
            @test -1 <= max_angular_momentum(psp) <= 5
            @test 0 <= n_radials(BetaProjector(), psp)
            @test 0 <= n_radials(ChiProjector(), psp)
            @test 0 <= valence_charge(file)
            @test isa(is_norm_conserving(psp), Bool)
            @test isa(is_ultrasoft(psp), Bool)
            @test isa(is_paw(psp), Bool)
            @test isa(has_spin_orbit(psp), Bool)

            for quantity in QUANTITIES
                @test isa(has_quantity(quantity, psp), Bool)
            end

            for l in angular_momenta(psp)
                @test size(get_quantity(BetaCoupling(), psp, l)) ==
                      (n_radials(BetaProjector(), psp, l),
                       n_radials(BetaProjector(), psp, l))
                @test eltype(get_quantity(BetaCoupling(), psp, l)) <: Real
                for n in 1:n_radials(BetaProjector(), psp, l)
                    @test get_quantity(BetaCoupling(), psp, l, n, n) ==
                          get_quantity(BetaCoupling(), psp, l, n)
                    for m in (n + 1):n_radials(BetaProjector(), psp, l)
                        @test isa(get_quantity(BetaCoupling(), psp, l, n, m), Real)
                        @test get_quantity(BetaCoupling(), psp, l, n, m) ==
                              get_quantity(BetaCoupling(), psp, l, m, n)
                    end
                end
            end

            K = rand(3) .+ eps(Float64)
            K_rot = rotate_vector(random_versor(), K)

            for quantity in [BetaProjector(), ChiProjector()]
                if has_quantity(quantity, psp)
                    for l in angular_momenta(psp), n in 1:n_radials(quantity, psp, l)
                        rmin = isa(psp, NumericPsP) ? first(psp.r) : 0
                        rmax = cutoff_radius(quantity, psp, l, n)
                        if !isnothing(rmax)
                            R = rescale_vector(rand(3), rmin, rmax)
                            R_rot = rotate_vector(random_versor(), R)

                            @test psp_quantity_evaluator(RealSpace(), quantity, psp, l, n)(R) ≈
                                psp_quantity_evaluator(RealSpace(), quantity, psp, l, n)(R_rot)
                            @test psp_quantity_evaluator(FourierSpace(), quantity, psp, l, n)(K) ≈
                                psp_quantity_evaluator(FourierSpace(), quantity, psp, l, n)(K_rot)
                        end
                    end
                end
            end

            for quantity in [LocalPotential(), ValenceChargeDensity(), CoreChargeDensity()]
                if has_quantity(quantity, psp)
                    rmin = isa(psp, NumericPsP) ? first(psp.r) : 0
                    rmax = cutoff_radius(quantity, psp)
                    if !isnothing(rmax)
                        R = rescale_vector(rand(3), rmin, rmax)
                        R_rot = rotate_vector(random_versor(), R)
                        @test psp_quantity_evaluator(RealSpace(), quantity, psp)(R) ≈
                            psp_quantity_evaluator(RealSpace(), quantity, psp)(R_rot)
                        @test psp_quantity_evaluator(FourierSpace(), quantity, psp)(K) ≈
                            psp_quantity_evaluator(FourierSpace(), quantity, psp)(K_rot)
                    end
                end
            end

            @test relativistic_treatment(psp) in (:scalar, :full)
            @test formalism(psp) in
                  (NormConservingPsP, UltrasoftPsP, ProjectorAugmentedWavePsP)
        end
    end
end
