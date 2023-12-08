function model(P,f)
    
    A  = 3+randn()*3*.2
    ## n part
    
    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A))  #compute sqrt(h^2Ωnoise(fi))
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    Data = Array(view(mean(c1, dims=2), :, 1))  #mean over chunks
    
    Data_noise = Data[1:970]
    f_noise = f[1:970]
    return  log10.(Data_noise), log10.(f_noise) #, log10.(f_total), log10.(stds_omega.^2)
end

