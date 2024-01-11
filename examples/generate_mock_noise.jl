using Pkg
using Revise
Pkg.activate("/home/zaldivar/Documents/Androniki/Github/GWBackFinder2.jl")
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter


#%%
### Define the frequency range
f = range(start=3 * 1e-5, stop=0.5, step=1e-6) 

### Split in 27 bins. idx27 shows in which bin the frequency belongs to and logbins_27 are the boundaries of the bins . After we coarse grain the frequencies from f = 10−3 Hz to the maximum frequency
### fmax = 0.5 Hz), i.e., we bin them in 1000 intervals of equal log-spacing. idx shows in which of the 1000 bins the frequency belongs to and logbins the boundaries of the 1000 bins.
idx,idx27,logbins_27,logbins,f_filtered=GWBackFinder.binning(f)

### Define breakpoints, i.e., read the boundaries of the 27 bins.
Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21, Sb22, Sb23, Sb24, Sb25, Sb26 = logbins_27[1], logbins_27[2], logbins_27[3], logbins_27[4], logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12], logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18], logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], logbins_27[24], logbins_27[25], logbins_27[26]
#%%
### Define values for parameters randomly. 10000000 different z (z -> 27 slopes, 1 amplitude, 1 noise parameter A) 
z= npzread("/data/users/Androniki/Dani/z_noise.npy")
#%%

### Look at an example and plot it 
data,freq=GWBackFinder.model_noise(f,z[1])
plt.clf()
plt.plot(freq, data[1:2970])
xlabel("f")
ylabel("Ω_GW")

plt.axvline(log10(Sb13))

#plt.axvline(log10(0.001))
plt.axvline(log10(0.003))
#plt.axhline(-14+z[100,26]*(-6-(-14)) )
"""
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

### Generate and save data
@showprogress for i in 1:1000000
    Data_total= GWBackFinder.model_noise(f,z[i])
    #print(length(Data_total))
    GWBackFinder.write_sample(Data_total,"/data/users/Androniki/Dani_new_noise/$(i-1).jld2")
end
