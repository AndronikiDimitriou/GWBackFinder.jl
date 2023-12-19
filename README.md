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
Pkg.instantiate()
```

Import other packages 
```
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter
```

# Examples

## Mock signal

```
f = range(start=3 * 1e-5, stop=0.5, step=1e-6) # define the frequency range

idx,idx27,logbins_27,logbins,f_filtered=GWBackFinder.binning(f) # define binned frequencies

Sb1 , Sb2 , Sb3 , Sb4 , Sb5 , Sb6 , Sb7 , Sb8 , Sb9 , Sb10 , Sb11 , Sb12 ,
Sb13 , Sb14 , Sb15 , Sb16 , Sb17 , Sb18 , Sb19 , Sb20 , Sb21 , Sb22 , Sb23 ,
Sb24 , Sb25 , Sb26 , Sb27  = logbins_27[1],logbins_27[2], logbins_27[3], logbins_27[4],
logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12],
logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18],
logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], 
logbins_27[24], logbins_27[25], logbins_27[26], logbins_27[27] # define breakpoints
```

### Prior
```
z=rand(1000000,29) # uniform prior
```

### Generate and save data

```
@showprogress for i in 1:size(z)[1]
    Data_total,freq = GWBackFinder.model(z[i, :],f,idx,f_filtered,logbins_27,logbins, Sb1, Sb2, Sb3,
    Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10, Sb11, Sb12,
    Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20, Sb21,
    Sb22, Sb23, Sb24, Sb25, Sb26)
    GWBackFinder.write_sample(Data_total,"$i.jld2")
end
```

### Plot an example 

```
data,freq=GWBackFinder.model(z[100, :],f,idx,f_filtered,logbins_27,logbins, Sb1, Sb2, Sb3, Sb4, Sb5, Sb6, Sb7, Sb8, Sb9, Sb10,\
Sb11, Sb12, Sb13, Sb14, Sb15, Sb16, Sb17, Sb18, Sb19, Sb20  , Sb21, Sb22, Sb23, Sb24, Sb25, Sb26)
plt.clf()
plt.plot(freq, data[1][1])
plt.show()
plt.gcf()
```
![alt text]([http://url/to/img.png](https://ibb.co/4RT16z5)https://ibb.co/4RT16z5)
