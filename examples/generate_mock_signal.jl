using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder.jl")
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter


#%%
f = range(start=3 * 1e-5, stop=0.5, step=1e-6)
idx,idx27,logbins_27,logbins,f_filtered=GWBackFinder.binning(f)
Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21, Sb22, Sb23, Sb24, Sb25, Sb26, Sb27  = logbins_27[1], logbins_27[2], logbins_27[3], logbins_27[4], logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12], logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18], logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], logbins_27[24], logbins_27[25], logbins_27[26], logbins_27[27]

#z= npzread("/home/zaldivar/Documents/Androniki/Github/GWBackFinder/src/GWBackFinder/data/thetas_for_julia.npy")

z=rand(1000000,29)

data,freq=GWBackFinder.model(z[100, :],f,idx,f_filtered,logbins_27,logbins, Sb1, Sb2  , Sb3  , Sb4  , Sb5  , Sb6  , Sb7  , Sb8  , Sb9  , Sb10  , Sb11  , Sb12  , Sb13  , Sb14  , Sb15  , Sb16  , Sb17  , Sb18  , Sb19  , Sb20  , Sb21  , Sb22  , Sb23  , Sb24  , Sb25  , Sb26  )
plt.clf()
plt.plot(freq, data[1][1])
xlabel("f")
ylabel("Ω_GW")

"""
plt.axvline(log10(Sb1))
plt.axvline(log10(Sb2))
plt.axvline(log10(Sb3))
plt.axvline(log10(Sb4))
plt.axvline(log10(Sb5))
plt.axvline(log10(Sb6))
plt.axvline(log10(Sb7))
plt.axvline(log10(Sb8))
plt.axvline(log10(Sb9))
plt.axvline(log10(Sb10))
plt.axvline(log10(Sb11))
plt.axvline(log10(Sb12))
plt.axvline(log10(Sb13))
plt.axvline(log10(Sb14))
plt.axvline(log10(Sb15))
plt.axvline(log10(Sb16))
plt.axvline(log10(Sb17))
plt.axvline(log10(Sb18))
plt.axvline(log10(Sb19))
plt.axvline(log10(Sb20))
plt.axvline(log10(Sb21))
plt.axvline(log10(Sb22))
plt.axvline(log10(Sb23))
plt.axvline(log10(Sb24))
plt.axvline(log10(Sb25))
plt.axvline(log10(Sb26))
plt.axvline(log10(Sb27))
"""

plt.show()
plt.gcf()

#%%

@showprogress for i in 1:size(z)[1]
    Data_total,freq = GWBackFinder.model(z[i, :],f,idx,f_filtered,logbins_27,logbins, Sb1, Sb2, Sb3,
    Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12,
    Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21,
    Sb22, Sb23, Sb24, Sb25, Sb26)
    GWBackFinder.write_sample(Data_total,"/data/users/Androniki/Dani2/$i.jld2")
end

