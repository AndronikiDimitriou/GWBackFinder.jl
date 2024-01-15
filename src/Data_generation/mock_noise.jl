function model_noise_train_data(f,P)
    """
    Model function that generates data based on given parameters. We assume that TT-channel is dominated by noise in the frequencies below 0.003 Hz. 

    Parameters:
    - P: Noise paramerer.
    - f: Frequency vector.

    Returns:
    - log10(Data_noise): Logarithm of the generated data.
    - log10(f_noise): Logarithm of the corresponding frequency vector.
    """
    
    # Generate random noise parameter A from a gaussian distribution with mean 3 and std 0.6.
    A  = 3+randn()*3*.2

    # Noise generation as described in https://arxiv.org/pdf/2009.11845.pdf
    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A)) 
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    Data = Array(view(mean(c1, dims=2), :, 1))  #mean over chunks
    
    Data_noise = Data[1:2970]

    return  log10.(Data_noise) 
end


function model_noise_plot(f,P)
    """
    Model function that generates data based on given parameters. We assume that TT-channel is dominated by noise in the frequencies below 0.003 Hz. 

    Parameters:
    - P: Noise paramerer.
    - f: Frequency vector.

    Returns:
    - log10(Data_noise): Logarithm of the generated data.
    - log10(f_noise): Logarithm of the corresponding frequency vector.
    """
    
    # Generate random noise parameter A from a gaussian distribution with mean 3 and std 0.6.
    A  = 3+randn()*3*.2

    # Noise generation as described in https://arxiv.org/pdf/2009.11845.pdf
    stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, P, A)) 
    stds_n_cuda = CuArray(Float32.(stds_n))

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise3 = stds_n_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_n), 94)
    CUDA.@sync samples_noise4 = stds_n_cuda .* std_eps


    CUDA.@sync c1 = (samples_noise3 .^ 2 + samples_noise4 .^ 2) / 2

    Data = Array(view(mean(c1, dims=2), :, 1))  #mean over chunks
    
    Data_noise = Data[1:2970]

    f_noise = f[1:2970]
    return  log10.(Data_noise) , log10.(f_noise)
end