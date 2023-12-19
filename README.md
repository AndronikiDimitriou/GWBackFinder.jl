# GWBackFinder

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AndronikiDimitriou.github.io/GWBackFinder.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AndronikiDimitriou.github.io/GWBackFinder.jl/dev/)
[![Build Status](https://github.com/AndronikiDimitriou/GWBackFinder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AndronikiDimitriou/GWBackFinder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AndronikiDimitriou/GWBackFinder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AndronikiDimitriou/GWBackFinder.jl)

# Installation


Clone/download this repository and then start Julia in the package folder, then type in the REPL:

using Pkg
Pkg.activate("./")
Pkg.instantiate()
using GWBackFinder
using PyPlot
using NPZ
using JLD2
using ProgressMeter

# Examples

## Mock signal
f = range(start=3 * 1e-5, stop=0.5, step=1e-6)
idx,idx27,logbins_27,logbins,f_filtered=GWBackFinder.binning(f)
Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3 = logbins_27[1], logbins_27[2], logbins_27[3], logbins_27[4], logbins_27[5], logbins_27[6], logbins_27[7], logbins_27[8], logbins_27[9], logbins_27[10], logbins_27[11], logbins_27[12], logbins_27[13], logbins_27[14], logbins_27[15], logbins_27[16], logbins_27[17], logbins_27[18], logbins_27[19], logbins_27[20], logbins_27[21], logbins_27[22], logbins_27[23], logbins_27[24], logbins_27[25], logbins_27[26], logbins_27[27]

