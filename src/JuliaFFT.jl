module JuliaFFT

export fft

function fft(xs)
    """
        fft(xs)

    Computes the discrete Fourier transform of the array `xs`.

    ``y_k = \\sum_{j=1}^N \\exp\\left[ - 2 \\pi i \\frac{(j-1)(k-1)}{N} \\right] x_j``

    The implementation is fully generic; as long as the input types have
    addition and multiplication (by complex numbers) defined, the algorithm will
    work.  Not particularly fast (about 2ms for FFT of 1024 element array on my
    modern-ish MacBook).

    Currently only works for 1D arrays.

    """
    n = size(xs,1)

    if n == 1
        return xs
    else
        xso = xs[1:2:end]
        xse = xs[2:2:end]

        ffo = fft(xso)
        ffe = fft(xse)

        runity = exp.(-2im*pi*range(0, stop=n/2-1)/n)

        return vcat(ffo .+ runity .* ffe, ffo .- runity.*ffe)
    end
end

end # module
