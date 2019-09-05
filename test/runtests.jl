using JuliaFFT
using Test: @testset, @test

@testset "Known-output FFT" begin
    xs = [1.0, -1.0, 1.0, -1.0]
    ys = fft(xs)
    @test all(ys .== [0.0, 0.0, 4.0, 0.0])

    xs = zeros(8)
    xs[1] = 1.0
    ys = fft(xs)
    @test all(ys .== ones(8))

    xs = zeros(4)
    xs[2] = 1.0
    ys = fft(xs)
    ys_true = [1.0+0im, 0.0-1im, -1.0+0im, 0.0+1.0im]
    @test all(abs.(ys .- ys_true) .< 1e-8)
end
