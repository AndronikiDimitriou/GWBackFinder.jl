
function binning(f)
    """
    Bins the input array 'f' into two sets of indices based on logarithmic binning.

    Parameters:
    - f: Input array.

    Returns:
    - idx: Indices for 1000 logarithmic bins.
    - idx_27: Indices for 27 logarithmic bins.
    - logbins_27: Logarithmic bin edges for 27 bins.
    - logbins: Logarithmic bin edges for 1000 bins.
    - f_filtered: Filtered array 'f' excluding the first 2970 elements.
    """
    # Filter the array 'f' starting from the 2971st element
    f_filtered = f[2971:end]
    # Create 1000 logarithmic bins and get indices for each element in 'f_filtered'
    logbins = 10 .^ range(log10(0.003), stop=log10(0.5), length=1000)
    idx = [searchsortedlast(logbins, f_filtered[i]) for i in eachindex(f_filtered)]
    idx[1]=1
    # Create 27 logarithmic bins and get indices for each element in 'f'
    logbins_27 = 10 .^ range(log10(3 * 1e-5), stop=log10(0.5), length=27)
    idx_27 = [searchsortedlast(logbins_27, f[i]) for i in eachindex(f)]
    idx_27[1]=1
    return idx,idx_27,logbins_27,logbins,f_filtered
end


function write_sample(Data_total, path)
    """
    Writes the input 'Data_total' to a file using the JLD2 format.

    Parameters:
    - Data_total: Data to be saved.
    - path: Path to save the data.
    """
    JLD2.@save(path,data= Data_total)
end