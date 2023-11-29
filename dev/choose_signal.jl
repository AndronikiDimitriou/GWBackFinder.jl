using Pkg
using Revise
Pkg.activate(".")
using GWBackFinder
using PyPlot

#%%
f = range(start=3 * 1e-5, stop=0.5, step=1e-6)

idx,idx28,logbins_28,logbins,f_filtered=GWBackFinder.binning(f)
Sb1_new3, Sb2_new3, Sb3_new3, Sb4_new3, Sb5_new3, Sb6_new3, Sb7_new3, Sb8_new3, Sb9_new3, Sb10_new3, Sb11_new3, Sb12_new3, Sb13_new3, Sb14_new3, Sb15_new3, Sb16_new3, Sb17_new3, Sb18_new3, Sb19_new3, Sb20_new3, Sb21_new3, Sb22_new3, Sb23_new3, Sb24_new3, Sb25_new3, Sb26_new3, Sb27_new3 = logbins_28[1], logbins_28[2], logbins_28[3], logbins_28[4], logbins_28[5], logbins_28[6], logbins_28[7], logbins_28[8], logbins_28[9], logbins_28[10], logbins_28[11], logbins_28[12], logbins_28[13], logbins_28[14], logbins_28[15], logbins_28[16], logbins_28[17], logbins_28[18], logbins_28[19], logbins_28[20], logbins_28[21], logbins_28[22], logbins_28[23], logbins_28[24], logbins_28[25], logbins_28[26], logbins_28[27]

z_power_law_test = GWBackFinder.zPowerLaw([-11,0.6])
z_peak_test=GWBackFinder.zPeak([-9.,0.002,0.2])
z_wiggly_test=GWBackFinder.zWiggly([-10,2,0.1])
z_Broken_powerlaw_test=GWBackFinder.zBroken_powerlaw([-11,-4,2/3,0.002])
z_Double_peaks_test=GWBackFinder.zDouble_peaks([-11.,-10.,0.001,0.01,0.1,0.1])
z_three_peaks_test=GWBackFinder.zThree_peaks([-10.,-10.,-10.,5*10^(-4),2*10^(-3),8*10^(-3),0.1,0.1,0.1])
z_noise=GWBackFinder.zNoise([15.,3.])

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_power_law_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_peak_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_wiggly_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_Broken_powerlaw_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_Double_peaks_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

plt.clf()
Data,f_total=GWBackFinder.different_signals(z_three_peaks_test,z_noise,"TT",f,f_filtered,logbins,idx)
plt.plot(f_total, Data)
plt.show()
plt.gcf()

