@testset "Utilities" begin

    @testset "Homogenization" begin
        @polyvar x y z
        @test Utilities.ishomogenous([x^2+y^2+x*y, x^5])
        @test Utilities.ishomogenous([x^2+y^2+x*y, x^4+1]) == false

        #test weighted degree
        @test Utilities.ishomogenous(x^3+x*y, [(x, 2), (y, 4)])
        @test Utilities.homogenize(x+x*y, [(x, 2), (y, 4)], z) == x*z^4+x*y

        @test Utilities.ishomogenous(Utilities.homogenize([x^2+y^2+x*y, x^4+1]))
        @test Utilities.ishomogenous([x^2+z^2 + y, x^4+z^4], [x,z]) == false
        @test Utilities.ishomogenous(Utilities.homogenize([x^2+z^2*y, x^4+z^4*y], [x,z]), [x,z]) == true


        @polyvar x y a b

        g = [x^2+y^2, x + y + 1, y]
        f = [2*a, a^2-2*a*b]
        @test Utilities.compose(g, f) isa Utilities.Composition
        @test g ∘ f isa Utilities.Composition
        @test Utilities.ishomogenous(g ∘ f) == false
        @test Utilities.ishomogenous(Utilities.homogenize(g ∘ f))
    end

    @testset "Misc" begin
        A = rand(Complex{Float64}, 12, 12)
        b = rand(Complex{Float64}, 12)
        C, d = copy(A), copy(b)
        @test norm(Utilities.solve!(C, d) - A \ b) < 1e-10

        A = rand(Complex{Float64}, 15, 12)
        b = rand(Complex{Float64}, 15)
        C, d = copy(A), copy(b)
        Utilities.solve!(C, d)
        @test norm(d[1:12] - A \ b) < 1e-10
    end
end
