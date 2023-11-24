
function binning()
    ## 1000 bins
    f = range(start=3 * 1e-5, stop=0.5, step=1e-6)
    f_filtered = f[971:end]
    logbins = 10 .^ range(log10(0.001), stop=log10(0.5), length=1000)
    idx = [searchsortedlast(logbins, f_filtered[i]) for i in eachindex(f_filtered)]

    ## 28 bins
    logbins_28 = 10 .^ range(log10(3 * 1e-5), stop=log10(0.5), length=28)
    idx_28 = [searchsortedlast(logbins_28, f[i]) for i in eachindex(f)]
    idx_28[1]=1
    return idx,idx_28,logbins_28
end