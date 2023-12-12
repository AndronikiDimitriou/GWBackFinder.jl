function Omega_powerlaw(A_1,f,gamma1)
    return 10^A_1.*(f./0.001).^gamma1
end

function Omega_peak(A_1,f,ftt,Delta1)
    return 10^A_1.*exp.(-(log10.(f./ftt)).^2/Delta1.^2)
end

function Omega_wiggly(A_1,f,Delta,fw)
    return 10^(A_1) * 10 .^ (sin.(Delta * log10.(f ./ fw)))
end

function Omega_double_peak(A_1,A_2,f,f1,f2,Delta1,Delta2)
    return 10^(A_1)*exp.(-(log10.(f./f1)).^2/Delta1^2)+10^(A_2)*exp.(-(log10.(f./f2)).^2/Delta2^2)
end

function Omega_broken_powerlaw(A_002, f11, gamma1, gamma2, ftt)
    term = 10^A_002 * ifelse.(ftt .- f11 .> 0, (f11 ./ ftt) .^ gamma1, (f11 ./ ftt) .^ gamma2)
    return term
end

function Omega_three_peaks(A_1,A_2,A_3,f,f1,f2,f3,Delta1,Delta2,Delta3)
    return 10^(A_1) * exp.(-(log10.(f ./ f1)).^2 / Delta1^2) .+
    10^(A_2) * exp.(-(log10.(f ./ f2)).^2 / Delta2^2) .+
    10^(A_3) * exp.(-(log10.(f ./ f3)).^2 / Delta3^2)
end

struct Zeta{T} # this is a structure, modify it with the parameters you need
    z_1::AbstractVector{T}
    z_2::AbstractVector{T}
    z_3::AbstractVector{T}
    z_4::AbstractVector{T}
    z_5::AbstractVector{T}
    z_6::AbstractVector{T}
end

abstract type zType end

struct zPowerLaw{T<:Real} <: zType
    val::PyList{Any}{T}
end

struct zPeak{T<:Real} <: zType
    val::AbstractVector{T}
end

struct zWiggly{T<:Real} <: zType
    val::AbstractVector{T}
end

struct zBroken_powerlaw{T<:Real} <: zType
    val::AbstractVector{T}
end

struct zThree_peaks{T<:Real} <: zType
    val::AbstractVector{T}
end

struct zDouble_peaks{T<:Real} <: zType
    val::AbstractVector{T}
end

struct zTest{T<:Real} <: zType
    val::AbstractVector{T}
end


struct zNoise{T<:Real} <: zType
    val::AbstractVector{T}
end


function different_signals(z::zPowerLaw,z_N::zNoise, channel::String, f,f_filtered, logbins:: AbstractVector, idx::AbstractVector)
    P,A = z_N.val[1],z_N.val[2]
    stds_omega = sqrt.(Omega_powerlaw(z.val[1],f,z.val[2])) 
    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end


function different_signals(z::zPeak,z_N::zNoise,channel::String, f,f_filtered, logbins::AbstractVector, idx::AbstractVector)
    
    P,A = z_N.val[1],z_N.val[2]
    stds_omega = sqrt.(Omega_peak(z.val[1],f,z.val[2],z.val[3])) 
    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"

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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end

function different_signals(z::zWiggly,z_N::zNoise,channel::String, f,f_filtered, logbins::AbstractVector,idx::AbstractVector)

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_wiggly(z.val[1],f,z.val[2],z.val[3])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"

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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end

function different_signals(z::zDouble_peaks,z_N::zNoise,channel::String, f,f_filtered,logbins::AbstractVector, idx::AbstractVector)

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_double_peak(z.val[1],z.val[2],f,z.val[3],z.val[4],z.val[5],z.val[6])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"

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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end

function different_signals(z::zBroken_powerlaw,z_N::zNoise,channel::String, f,f_filtered, logbins::AbstractVector, idx::AbstractVector)
    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_broken_powerlaw(z.val[1],f,z.val[2],z.val[3],z.val[4])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"

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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end

function different_signals(z::zThree_peaks,z_N::zNoise,channel::String, f,f_filtered,logbins::AbstractVector, idx::AbstractVector)

    P,A = z_N.val[1],z_N.val[2]

    stds_omega = sqrt.(Omega_three_peaks(z.val[1],z.val[2],z.val[3],f,z.val[4],z.val[5],z.val[6],z.val[7],z.val[8],z.val[9])) 

    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"

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
            data_filtered = Data[971:end][idx.==i]
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
    #"""

    elseif channel=="TT"
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
            data_filtered = Data[971:end][idx.==i]
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
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)


    return  log10.(Data_total), log10.(f_total)
end

function different_signals(z::zType)
    @error "Not Implemented"
end


