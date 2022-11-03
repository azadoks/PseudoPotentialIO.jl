# @testset "Load all HGH" begin
#     for (root, dirs, files) in walkdir("./hgh/"), file in files
#         psp = HghPsP(joinpath(root, file))
#         @test isa(psp, Psp8HGH)
#         @test format(psp) == "HGH"
#     end
# end
