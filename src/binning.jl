
function binning(f)

    ## 1000 bins
    f_filtered = f[971:end]
    logbins = 10 .^ range(log10(0.001), stop=log10(0.5), length=1000)
    idx = [searchsortedlast(logbins, f_filtered[i]) for i in eachindex(f_filtered)]

    ## 27 bins
    logbins_27 = 10 .^ range(log10(3 * 1e-5), stop=log10(0.5), length=27)
    idx_27 = [searchsortedlast(logbins_27, f[i]) for i in eachindex(f)]
    idx_27[1]=1
    return idx,idx_27,logbins_27,logbins,f_filtered
end