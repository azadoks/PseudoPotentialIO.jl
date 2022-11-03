# @testset "Load all PSP8" begin
# @testset "Load all PSP8" begin
#     for (root, dirs, files) in walkdir("./psp8/"), file in files
#         psp = Psp8PsP(joinpath(root, file))
#         @test isa(psp, Psp8PsP)
#         @test format(psp) == "PSP v8"
#     end
# end
