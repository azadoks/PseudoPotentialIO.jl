# Write and read UPF files from several tasks at once, with other tasks allocating (so that
# garbage collections start on other threads mid-call), and compare every output with a
# serial one. EzXML's `prettyprint` freed its buffer while copying the text out under these
# conditions: 4-8% of the files written were garbage (PseudoPotentialIO #7). Meaningful
# only with several threads (`julia -t 4`; `Pkg.test` passes the thread count on).

# Number of results of `f(input)` that differ from `reference` (`same`), from `ntasks`
# tasks each calling `f` on every input `rounds` times, while `npressure` tasks allocate
function count_concurrent_mismatches(f, inputs, references; same=isequal, rounds=3,
                                     ntasks=4, npressure=2)
    nbad = Threads.Atomic{Int}(0)
    stop = Threads.Atomic{Bool}(false)
    pressure = map(1:npressure) do _
        Threads.@spawn while !stop[]
            garbage = [rand(100) for _ in 1:1000]
            GC.safepoint()
        end
    end
    workers = map(1:ntasks) do _
        Threads.@spawn for _ in 1:rounds, (input, reference) in zip(inputs, references)
            same(f(input), reference) || Threads.atomic_add!(nbad, 1)
        end
    end
    foreach(wait, workers)
    stop[] = true
    foreach(wait, pressure)
    return nbad[]
end

@testset "UPF from several threads" begin
    filepaths = [UPF2_CASE_FILEPATHS["Dy.GGA-PBE-paw-v1.0.UPF"],
                 UPF2_CASE_FILEPATHS["Si.pbe-n-rrkjus_psl.1.0.0.UPF"],
                 UPF2_CASE_FILEPATHS["Fe.upf"],
                 UPF1_CASE_FILEPATHS["ag_lda_v1.4.uspp.F.UPF"]]
    psps = load_psp_file.(filepaths)
    upf_text(psp) = sprint(io -> PseudoPotentialIO.save_psp(io, psp))
    texts = upf_text.(psps)

    @testset "save_psp_file" begin
        outdir = mktempdir(; cleanup=true)
        counter = Threads.Atomic{Int}(0)
        function write_file(psp)
            path = joinpath(outdir, "$(Threads.atomic_add!(counter, 1)).upf")
            save_psp_file(path, psp)
            return read(path, String)
        end
        @test count_concurrent_mismatches(write_file, psps, texts) == 0
    end

    @testset "load_psp_file" begin
        # read back the files written (texts), compared through their text written again:
        # every field, to the last digit written
        read_text(text) = upf_text(UpfFile(IOBuffer(text)))
        references = read_text.(texts)
        @test count_concurrent_mismatches(read_text, texts, references) == 0
    end
end
