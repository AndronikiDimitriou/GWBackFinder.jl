


function Omega_powerlaw(A_1,f,gamma1)
    """
    Calculate the power-law spectrum for a given amplitude and exponent.

    Parameters:
    - A_1: Amplitude parameter.
    - f: Frequency vector.
    - gamma1: Exponent parameter.

    Returns:
    - Spectrum values calculated based on the power-law formula.
    """
    # Calculate power-law spectrum using the given parameters
    return 10^A_1.*(f./0.001).^gamma1
end

function Omega_peak(A_1,f,ftt,Delta1)
    """
    Calculate a spectrum with a peak for a given amplitude, peak frequency, and width.

    Parameters:
    - A_1: Amplitude parameter.
    - f: Frequency vector.
    - ftt: Peak frequency parameter.
    - Delta1: Width parameter.

    Returns:
    - Spectrum values calculated based on the peak formula.
    """
    # Calculate the spectrum with a peak using the given parameters

    return 10^A_1.*exp.(-(log10.(f./ftt)).^2/Delta1.^2)
end

function Omega_wiggly(A_1,f,Delta,fw)
    """
    Calculate a wiggly spectrum for a given amplitude, frequency, and modulation parameters.

    Parameters:
    - A_1: Amplitude parameter.
    - f: Frequency vector.
    - Delta: Modulation amplitude parameter.
    - fw: Modulation frequency parameter.

    Returns:
    - Spectrum values calculated based on the wiggly formula.
    """
    # Calculate the wiggly spectrum using the given parameters
    return 10^(A_1) * 10 .^ (sin.(Delta * log10.(f ./ fw)))
end

function Omega_double_peak(A_1,A_2,f,f1,f2,Delta1,Delta2)
    """
    Calculate a spectrum with two Gaussian peaks for given parameters.

    Parameters:
    - A_1, A_2: Amplitude parameters for the two peaks.
    - f: Frequency vector.
    - f1, f2: Peak frequencies for the two peaks.
    - Delta1, Delta2: Modulation amplitude parameters for the two peaks.

    Returns:
    - Spectrum values calculated based on the double peak formula.
    """
    # Calculate the double peak spectrum using the given parameters
    return 10^(A_1)*exp.(-(log10.(f./f1)).^2/Delta1^2)+10^(A_2)*exp.(-(log10.(f./f2)).^2/Delta2^2)
end

function Omega_broken_powerlaw(A_002, f11, gamma1, gamma2, ftt)
    """
    Calculate a broken power-law spectrum for given parameters.

    Parameters:
    - A_002: Amplitude parameter.
    - f11: Break frequency.
    - gamma1, gamma2: Power-law exponents for the two segments.
    - ftt: Transition frequency between the two power-law segments.

    Returns:
    - Spectrum values calculated based on the broken power-law formula.
    """
    # Calculate the broken power-law spectrum using the given parameters
    
    term = 10^A_002 * ifelse.(ftt .- f11 .> 0, (f11 ./ ftt) .^ gamma1, (f11 ./ ftt) .^ gamma2)
    return term
end

function Omega_three_peaks(A_1,A_2,A_3,f,f1,f2,f3,Delta1,Delta2,Delta3)
    """
    Calculate a spectrum with three Gaussian peaks for given parameters.

    Parameters:
    - A_1, A_2, A_3: Amplitude parameters for the three peaks.
    - f: Frequency vector.
    - f1, f2, f3: Peak frequencies for the three peaks.
    - Delta1, Delta2, Delta3: Modulation amplitude parameters for the three peaks.

    Returns:
    - Spectrum values calculated based on the three peak formula.
    """
    # Calculate the three peak spectrum using the given parameters
    
    return 10^(A_1) * exp.(-(log10.(f ./ f1)).^2 / Delta1^2) .+
    10^(A_2) * exp.(-(log10.(f ./ f2)).^2 / Delta2^2) .+
    10^(A_3) * exp.(-(log10.(f ./ f3)).^2 / Delta3^2)
end

struct Zeta{T} 
    z_1::Vector{T}
    z_2::Vector{T}
    z_3::Vector{T}
    z_4::Vector{T}
    z_5::Vector{T}
    z_6::Vector{T}
end

function Zeta(z_1::AbstractVector{T}, z_2::AbstractVector{T}, z_3::AbstractVector{T}, z_4::AbstractVector{T}, z_5::AbstractVector{T}, z_6::AbstractVector{T}) where {T<:Real}
    return Zeta{T}(z_1, z_2, z_3, z_4, z_5, z_6)
end

abstract type zType end

struct zPowerLaw{T<:Real} <: zType
    val::Vector{T}
end

function zPowerLaw(val::AbstractVector{T}) where {T<:Real}
    return zPowerLaw{T}(val)
end

struct zPeak{T<:Real} <: zType
    val::Vector{T}
end

function zPeak(val::AbstractVector{T}) where {T<:Real}
    return zPeak{T}(val)
end

struct zWiggly{T<:Real} <: zType
    val::Vector{T}
end

function zWiggly(val::AbstractVector{T}) where {T<:Real}
    return zWiggly{T}(val)
end

struct zBroken_powerlaw{T<:Real} <: zType
    val::Vector{T}
end

function zBroken_powerlaw(val::AbstractVector{T}) where {T<:Real}
    return zBroken_powerlaw{T}(val)
end

struct zThree_peaks{T<:Real} <: zType
    val::Vector{T}
end

function zThree_peaks(val::AbstractVector{T}) where {T<:Real}
    return zThree_peaks{T}(val)
end

struct zDouble_peaks{T<:Real} <: zType
    val::Vector{T}
end

function zDouble_peaks(val::AbstractVector{T}) where {T<:Real}
    return zDouble_peaks{T}(val)
end

struct zTest{T<:Real} <: zType
    val::Vector{T}
end

function zTest(val::AbstractVector{T}) where {T<:Real}
    return zTest{T}(val)
end

struct zNoise{T<:Real} <: zType
    val::Vector{T}
end

function zNoise(val::AbstractVector{T}) where {T<:Real}
    return zNoise{T}(val)
end

function different_signals(z::zPowerLaw,z_N::zNoise, f,f_filtered, logbins:: AbstractVector, idx::AbstractVector)
    """
    A function that computes signals based on a combination of power-law and noise components.

    Parameters:
        - `z::zPowerLaw`: A struct containing parameters for the power-law component.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

    Returns:
    - Logarithm of the total signal in the AA channel.
    - Logarithm of the total signal in the TT channel .
    - Logarithm of the frequencies corresponding to the AA channel.
    - Logarithm of the frequencies corresponding to the TT channel.

    The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    """

    CUDA.device_reset!
    P,A = z_N.val[1],z_N.val[2]
    stds_omega = sqrt.(Omega_powerlaw(z.val[1],f,z.val[2])) 
    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # Signal part. Following Data generation described in https://arxiv.org/pdf/2009.11845.pdf

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2

    # Noise part AA channel. Following Data generation described in https://arxiv.org/pdf/2009.11845.pdf

    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    # Add signal and noise 

    CUDA.@sync c = c1 + c2

    # Calculate mean over chunks

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    # Bin the data and frequencies after the 2971 element in 1000 bins and calculate a value for each of them using nverse variance weighting. (see https://arxiv.org/pdf/2009.11845.pdf)

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    # Concatenate  unbinned data before the 2971 element with the binned data. 

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    # Noise part TT channel. Following Data generation described in https://arxiv.org/pdf/2009.11845.pdf

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    # Add signal and noise 

    CUDA.@sync c = c1 + c2


    # Calculate mean over chunks

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    # Bin the data and frequencies after the 2971 element in 1000 bins and calculate a value for each of them using nverse variance weighting. (see https://arxiv.org/pdf/2009.11845.pdf)

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    # Concatenate  unbinned data before the 2971 element with the binned data. 

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    # Concatenate data with the value of the noise parameter P 

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)
end


function different_signals(z::zPeak,z_N::zNoise, f,f_filtered, logbins::AbstractVector, idx::AbstractVector)
    """
    A function that computes signals based on a combination of gaussian peak and noise components.

    Parameters:
        - `z::zPeak`: A struct containing parameters for the guassian peak.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

        Returns:
        - Logarithm of the total signal in the AA channel.
        - Logarithm of the total signal in the TT channel .
        - Logarithm of the frequencies corresponding to the AA channel.
        - Logarithm of the frequencies corresponding to the TT channel.

    The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    
    """
    CUDA.device_reset!

    P,A = z_N.val[1],z_N.val[2]
    stds_omega = sqrt.(Omega_peak(z.val[1],f,z.val[2],z.val[3])) 
    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2




    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2


    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)
end

function different_signals(z::zWiggly,z_N::zNoise, f,f_filtered, logbins::AbstractVector,idx::AbstractVector)
    """
    A function that computes signals based on a combination of wiggly signal and noise components.

    Parameters:
        - `z::zWiggly`: A struct containing parameters for the wiggly signal.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

        Returns:
        - Logarithm of the total signal in the AA channel.
        - Logarithm of the total signal in the TT channel .
        - Logarithm of the frequencies corresponding to the AA channel.
        - Logarithm of the frequencies corresponding to the TT channel.
        The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    
    """
    CUDA.device_reset!

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_wiggly(z.val[1],f,z.val[2],z.val[3])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2




    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2


    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)
end

function different_signals(z::zDouble_peaks,z_N::zNoise, f,f_filtered,logbins::AbstractVector, idx::AbstractVector)
    """
    A function that computes signals based on a combination of two gaussian peaks and noise components.

    Parameters:
        - `z::zDouble_peaks`: A struct containing parameters for the two gaussian peak.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

        Returns:
        - Logarithm of the total signal in the AA channel.
        - Logarithm of the total signal in the TT channel .
        - Logarithm of the frequencies corresponding to the AA channel.
        - Logarithm of the frequencies corresponding to the TT channel.
        The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    
    """
    CUDA.device_reset!

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_double_peak(z.val[1],z.val[2],f,z.val[3],z.val[4],z.val[5],z.val[6])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2




    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2


    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)

end

function different_signals(z::zBroken_powerlaw,z_N::zNoise, f,f_filtered, logbins::AbstractVector, idx::AbstractVector)
    """
    A function that computes signals based on a combination of a broken power-law and noise components.

    Parameters:
        - `z::zBroken_powerlaws`: A struct containing parameters for the broken power-law.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

        Returns:
        - Logarithm of the total signal in the AA channel.
        - Logarithm of the total signal in the TT channel .
        - Logarithm of the frequencies corresponding to the AA channel.
        - Logarithm of the frequencies corresponding to the TT channel.
        The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    
    """
    CUDA.device_reset!

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_broken_powerlaw(z.val[1],f,z.val[2],z.val[3],z.val[4])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2




    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2


    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)
end

function different_signals(z::zThree_peaks,z_N::zNoise, f,f_filtered,logbins::AbstractVector, idx::AbstractVector)
    """
    A function that computes signals based on a combination of three gaussian peaks and noise components.

    Parameters:
        - `z::zThree_peaks`: A struct containing parameters for the three gaussian peaks.
        - `z_N::zNoise`: A struct containing parameters for the noise component.
        - `f`: Frequency vector.
        - `f_filtered`: Filtered frequency vector.
        - `logbins::AbstractVector`: Vector of the edeges of 1000 logarithmically spaced bins.
        - `idx::AbstractVector`: Vector of indices showing in which bin the frequency belongs to.

        Returns:
        - Logarithm of the total signal in the AA channel.
        - Logarithm of the total signal in the TT channel .
        - Logarithm of the frequencies corresponding to the AA channel.
        - Logarithm of the frequencies corresponding to the TT channel.
        The function uses CUDA for parallel processing and generates random samples to compute the power-law and noise components separately.
    
    """
    CUDA.device_reset!

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_three_peaks(z.val[1],z.val[2],z.val[3],f,z.val[4],z.val[5],z.val[6],z.val[7],z.val[8],z.val[9])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2




    stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2


    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_AA = vcat(Data[1:2970], weighted_data)
    f_total_AA = vcat(f[1:2970], weighted_f)

#"""

    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    CUDA.@sync c = c1 + c2

    Data = Array(view(mean(c, dims=2), :, 1))  #mean over chunks

    weight_norm = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1))

        weight_norm[i] = res
    end

    weighted_data = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        data_filtered = Data[2971:end][idx.==i]
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

        weighted_data[i] = sum(data_filtered .* w)
    end

#"""
    weighted_f = zeros(length(logbins))
    @threads for i in eachindex(logbins)
        num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, P, A)) .^ (-1)
        w = num ./ weight_norm[i]

    weighted_f[i] = sum(f_filtered[idx.==i] .* w)
    end

    Data_total_TT = vcat(Data[1:2970], weighted_data)
    f_total_TT = vcat(f[1:2970], weighted_f)
    #"""

    return  vcat(log10.(Data_total_AA),P), log10.(Data_total_TT), log10.(f_total_AA), log10.(f_total_TT)
end

function different_signals(z::zType)
    @error "Not Implemented"
end


