
function f26(S, S0, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, n21, n22, n23, n24, n25, n26, As,Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3)

    if S <= Sb1_new3

        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9) *
                 (Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (Sb5_new3 / S0)^(n5 - n4) *
                 (Sb4_new3 / S0)^(n4 - n3) * (Sb3_new3 / S0)^(n3 - n2) * (Sb2_new3 / S0)^(n2 - n1) * (Sb1_new3 / S0)^(n1 - n0) * (S / S0)^(n0)

    elseif Sb1_new3 < S <= Sb2_new3

        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (Sb5_new3 / S0)^(n5 - n4)*(Sb4_new3 / S0)^(n4 - n3) * (Sb3_new3 / S0)^(n3 - n2) * (Sb2_new3 / S0)^(n2 - n1) * (S / S0)^(n1)

    elseif Sb2_new3 < S <= Sb3_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (Sb5_new3 / S0)^(n5 - n4)*(Sb4_new3 / S0)^(n4 - n3) * (Sb3_new3 / S0)^(n3 - n2) * (S / S0)^(n2)

    elseif Sb3_new3 < S <= Sb4_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (Sb5_new3 / S0)^(n5 - n4)*(Sb4_new3 / S0)^(n4 - n3) * (S / S0)^(n3)

    elseif Sb4_new3 < S <= Sb5_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (Sb5_new3 / S0)^(n5 - n4)*(S / S0)^(n4)

    elseif Sb5_new3 < S <= Sb6_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (Sb6_new3 / S0)^(n6 - n5) * (S / S0)^(n5)

    elseif Sb6_new3 < S <= Sb7_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (Sb7_new3 / S0)^(n7 - n6) * (S / S0)^(n6)

    elseif Sb7_new3 < S <= Sb8_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (Sb8_new3 / S0)^(n8 - n7) * (S / S0)^(n7)

    elseif Sb8_new3 < S <= Sb9_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(Sb9_new3 / S0)^(n9 - n8) * (S / S0)^(n8)

    elseif Sb9_new3 < S <= Sb10_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (Sb10_new3 / S0)^(n10 - n9)*(S / S0)^(n9)

    elseif Sb10_new3 < S <= Sb11_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (Sb11_new3 / S0)^(n11 - n10) * (S / S0)^(n10)

    elseif Sb11_new3 < S <= Sb12_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (Sb12_new3 / S0)^(n12 - n11) * (S / S0)^(n11)

    elseif Sb12_new3 < S <= Sb13_new3
        result = As * (Sb13_new3 / S0)^(n13 - n12) * (S / S0)^(n12)

    elseif Sb13_new3 < S <= Sb14_new3
        result = As * (S / S0)^(n13)

    elseif Sb14_new3 < S <= Sb15_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (S / S0)^(n14)

    elseif Sb15_new3 < S <= Sb16_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (S / S0)^(n15)

    elseif Sb16_new3 < S <= Sb17_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (S / S0)^(n16)

    elseif Sb17_new3 < S <= Sb18_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17) * (S / S0)^(n17)

    elseif Sb18_new3 < S <= Sb19_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (S / S0)^(n18)

    elseif Sb19_new3 < S <= Sb20_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (S / S0)^(n19)

    elseif Sb20_new3 < S <= Sb21_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (S / S0)^(n20)

    elseif Sb21_new3 < S <= Sb22_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21) * (S / S0)^(n21)

    elseif Sb22_new3 < S <= Sb23_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21) * (Sb22_new3 / S0)^(n21 - n22) * (S / S0)^(n22)

    elseif Sb23_new3 < S <= Sb24_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21) *(Sb22_new3 / S0)^(n21 - n22) * (Sb23_new3 / S0)^(n22 - n23) * (S / S0)^(n23)

    elseif Sb24_new3 < S <= Sb25_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21) *(Sb22_new3 / S0)^(n21 - n22) * (Sb23_new3 / S0)^(n22 - n23) * (Sb24_new3 / S0)^(n23 - n24) * (S / S0)^(n24)

    elseif Sb25_new3 < S <= Sb26_new3
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21) *(Sb22_new3 / S0)^(n21 - n22) * (Sb23_new3 / S0)^(n22 - n23) * (Sb24_new3 / S0)^(n23 - n24) * (Sb25_new3 / S0)^(n24 - n25) * (S / S0)^(n25)

    else
        result = As * (Sb14_new3 / S0)^(n13 - n14) * (Sb15_new3 / S0)^(n14 - n15) * (Sb16_new3 / S0)^(n15 - n16) * (Sb17_new3 / S0)^(n16 - n17)*(Sb18_new3 / S0)^(n17 - n18) * (Sb19_new3 / S0)^(n18 - n19) * (Sb20_new3 / S0)^(n19 - n20) * (Sb21_new3 / S0)^(n20 - n21)*(Sb22_new3 / S0)^(n21 - n22) * (Sb23_new3 / S0)^(n22 - n23) * (Sb24_new3 / S0)^(n23 - n24) * (Sb25_new3 / S0)^(n24 - n25)*(Sb26_new3 / S0)^(n25 - n26) * (S / S0)^(n26)

    end
    return result
end


function model(z1,f,channel,idx,f_filtered,logbins_28,logbins,Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3)
    
    e  = 15+randn()*15*.2
    z  =-12 .+z1[1:27]*(12-(-12))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    z2 =-14+z1[27]*(-6-(-14)) 
    z3 = z1[28]
    print(z,z2,z3)
    ## n part
    stds_omega = [sqrt.(f26(fi, logbins_28[14], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11], z[12], z[13], z[14], z[15], z[16], z[17], z[18], z[19], z[20], z[21], z[22], z[23], z[24], z[25], z[26], z[27], 10^(z2), Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3)) for fi in f]
    stds_omega_cuda = CuArray(Float32.(stds_omega))

    # repeat for each chunk
    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise5 = stds_omega_cuda .* std_eps

    std_eps = CUDA.randn(Float32, length(stds_omega), 94)
    CUDA.@sync samples_noise6 = stds_omega_cuda .* std_eps

    CUDA.@sync c2 = (samples_noise5 .^ 2 + samples_noise6 .^ 2) ./ 2


    if channel=="AA"
        stds_n = sqrt.(Omega_noiseh2_AA.(f, 2.5 * 1e9, 3 * 1e8, e, z3))  #compute sqrt(h^2Ωnoise(fi))
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
            res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1))

            weight_norm[i] = res
        end

        weighted_data = zeros(length(logbins))
        @threads for i in eachindex(logbins)
            data_filtered = Data[971:end][idx.==i]
            num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1)
            w = num ./ weight_norm[i]

            weighted_data[i] = sum(data_filtered .* w)
        end

    #"""
        weighted_f = zeros(length(logbins))
        @threads for i in eachindex(logbins)
            num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1)
            w = num ./ weight_norm[i]

        weighted_f[i] = sum(f_filtered[idx.==i] .* w)
        end
    #"""

    elseif channel=="TT"
        stds_n = sqrt.(Omega_noiseh2_TT.(f, 2.5 * 1e9, 3 * 1e8, e, z3))  #compute sqrt(h^2Ωnoise(fi))
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
            res = sum((Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1))

            weight_norm[i] = res
        end

        weighted_data = zeros(length(logbins))
        @threads for i in eachindex(logbins)
            data_filtered = Data[971:end][idx.==i]
            num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1)
            w = num ./ weight_norm[i]

            weighted_data[i] = sum(data_filtered .* w)
        end

    #"""
        weighted_f = zeros(length(logbins))
        @threads for i in eachindex(logbins)
            num = (Omega_noiseh2_AA.(f_filtered[idx.==i], 2.5 * 1e9, 3 * 1e8, e, z3)) .^ (-1)
            w = num ./ weight_norm[i]

        weighted_f[i] = sum(f_filtered[idx.==i] .* w)
        end
    #"""
    end
    Data_total = vcat(Data[1:970], weighted_data)
    f_total = vcat(f[1:970], weighted_f)

    return  vcat((log10.(Data_total),e)), log10.(f_total)
end

