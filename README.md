# GWBackFinder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AndronikiDimitriou.github.io/GWBackFinder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AndronikiDimitriou.github.io/GWBackFinder.jl/dev/)
[![Build Status](https://github.com/AndronikiDimitriou/GWBackFinder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AndronikiDimitriou/GWBackFinder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AndronikiDimitriou/GWBackFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AndronikiDimitriou/GWBackFinder.jl)

# Installation


Clone/download this repository and then start Julia in the package folder, then type in the REPL:
```
using Pkg
Pkg.activate("./")
using GWBackFinder

```

# Examples

* ## Mock signal (see generate_mock_signal.jl)

Import other packages 
```
using PyPlot
using NPZ
using JLD2
using ProgressMeter
```

### Define the frequency range
```
f = range(start=3 * 1e-5, stop=0.5, step=1e-6) 
```
### Split in 26 bins. idx26 shows in which bin the frequency belongs to and logbins_26 are the boundaries of the bins . After we coarse grain the frequencies from $f = 3*10^{−3}$ Hz to the maximum frequency $f_{max} = 0.5$ Hz), i.e., we bin them in 1000 intervals of equal log-spacing. idx shows in which of the 1000 bins the frequency belongs to and logbins the boundaries of the 1000 bins.
```
idx,idx26,logbins_26,logbins,f_filtered=GWBackFinder.binning(f) 
```
### Define breakpoints
```
Sb1 , Sb2 , Sb3 , Sb4 , Sb5 , Sb6 , Sb7 , Sb8 , Sb9 , Sb10 , Sb11 , Sb12 ,
Sb13 , Sb14 , Sb15 , Sb16 , Sb17 , Sb18 , Sb19 , Sb20 , Sb21 , Sb22 , Sb23 ,
Sb24 , Sb25 , Sb26  = logbins_26[1],logbins_26[2], logbins_26[3], logbins_26[4],
logbins_26[5], logbins_26[6], logbins_26[7], logbins_26[8], logbins_26[9], logbins_26[10], logbins_26[11], logbins_26[12],
logbins_26[13], logbins_26[14], logbins_26[15], logbins_26[16], logbins_26[17], logbins_26[18],
logbins_26[19], logbins_26[20], logbins_26[21], logbins_26[22], logbins_26[23], 
logbins_26[24], logbins_26[25], logbins_26[26]
```

### Define values for parameters (z $\rightarrow$ 26 slopes, 1 amplitude, 1 noise parameter A)
```
z=rand(1000000,28) # for example uniform prior
```

### Generate and save data

Save as number.jld2, or use a different path

```
@showprogress for i in 1:size(z)[1]
    Data_total,freq = GWBackFinder.model(z[i, :],f,idx,f_filtered,logbins_26,logbins, Sb1, Sb2, Sb3,
    Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12,
    Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21,
    Sb22, Sb23, Sb24, Sb25)
    GWBackFinder.write_sample(Data_total,"$i.jld2") 
end
```

### Plot an example 

```
data,freq=GWBackFinder.model(z[100, :],f,idx,f_filtered,logbins_26,logbins, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9,
Sb10,Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20  , Sb21, Sb22, Sb23, Sb24, Sb25, Sb26) #to run it uncomment the frequency calculation in ./GWBackFinder.jl/src/Data_generation/mock_signal.jl
plt.clf()
plt.plot(freq, data[1:3970])
xlabel("f[Hz]")
ylabel("Ω_GW")
plt.show()
plt.title("Mock signal")
plt.gcf()
#%%
```
![plot](/examples/mock_signal.png)


* ## Mock noise (see generate_mock_signal_noise.jl)

Import other packages 
```
using PyPlot
using NPZ
using JLD2
using ProgressMeter
```

### Define the frequency range
```
f = range(start=3 * 1e-5, stop=0.5, step=1e-6) 
```
### Define binned frequencies
```
idx,idx26,logbins_26,logbins,f_filtered=GWBackFinder.binning(f) 
```

### Define values for noise parameter P
```
z=rand(1000000,1) # for example uniform prior
```

### Generate and save data

Save as number.jld2, or use a different path

```
@showprogress for i in 1:1000000
    Data_total= GWBackFinder.model_noise(f,z[i])
    #print(length(Data_total))
    GWBackFinder.write_sample(Data_total,"/data/users/Androniki/Dani_new_noise/$(i-1).jld2")
end
```

### Plot an example 

```
data,freq=GWBackFinder.model_noise(f,z[1]) #to run it uncomment the frequency calculation in ./GWBackFinder.jl/src/Data_generation/mock_noise.jl
plt.clf()
plt.plot(freq, data[1:2970])
xlabel("f[Hz]")
ylabel("Ω_GW")
plt.title("Mock noise")
plt.show()
plt.gcf()

#%%
```
![plot](/examples/mock_noise.png)


* ## Different signals (see choose_template_signals.jl)
Import other packages 

```
using GWBackFinder
using PyPlot
```

### Define the frequency range
```
f = range(start=3 * 1e-5, stop=0.5, step=1e-6)
```
### Split in 26 bins. idx26 shows in which bin the frequency belongs to and logbins_26 are the boundaries of the bins . After we coarse grain the frequencies from f = 10−3 Hz to the maximum frequency fmax = 0.5 Hz), i.e., we bin them in 1000 intervals of equal log-spacing. idx shows in which of the 1000 bins the frequency belongs to and logbins the boundaries of the 1000 bins.

```
idx,idx26,logbins_26,logbins,f_filtered = GWBackFinder.binning(f)
```

### Define different signals

```
z_power_law_test         =  GWBackFinder.zPowerLaw([-9,0.6])
z_peak_test              =  GWBackFinder.zPeak([-9.,0.002,0.2])
z_wiggly_test            =  GWBackFinder.zWiggly([-10,2,0.1])
z_Broken_powerlaw_test   =  GWBackFinder.zBroken_powerlaw([-11,-1,2/3,0.002])
z_Double_peaks_test      =  GWBackFinder.zDouble_peaks([-11.,-10.,0.001,0.01,0.1,0.1])
z_three_peaks_test       =  GWBackFinder.zThree_peaks([-10.,-10.,-10.,5*10^(-4),2*10^(-3),
                            8*10^(-3),0.1,0.1,0.1])
z_noise                  =  GWBackFinder.zNoise([15.,3.])


plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_power_law_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Powerlaw")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_peak_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Peak")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_wiggly_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Wiggly")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Broken_powerlaw_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.title("Broken powerlaw")
plt.legend()
plt.show()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ  = GWBackFinder.different_signals(z_Double_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.show()
plt.title("Double_peak")
plt.legend()
plt.gcf()

plt.clf()
Data_ΑΑ,Data_ΤΤ,f_total_ΑΑ,f_total_ΤΤ = GWBackFinder.different_signals(z_three_peaks_test,z_noise,f,f_filtered,logbins,idx)
plt.plot(f_total_ΑΑ, Data_ΑΑ[1:3970],label="AA channel")
plt.plot(f_total_ΤΤ, Data_ΤΤ,label="TT channel")
plt.title("Three peaks")
plt.show()
plt.legend()
plt.gcf()
```
![](/examples/peak.png)

