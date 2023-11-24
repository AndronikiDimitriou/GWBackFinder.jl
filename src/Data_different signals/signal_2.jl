function Omega_Powerlaw(A_1,f,gamma1)
    return 10^A_1.*(f./0.001).^gamma1
end

function Omega_Peak(A_1,f,ftt,Delta1)
    return 10^A_1.*exp.(-(log10.(f./ftt)).^2/Delta1.^2)
end

function Omega_wiggly(A_1,f,Delta,fw)
    return 10^(A_1) * 10 .^ (sin.(Delta * log10.(f ./ fw)))
end

function Omega_double_peak(A_1,A_2,f,f1,f2,Delta1,Delta2)
    return 10^(A_1)*exp.(-(log10.(f./f1)).^2/Delta1^2)+10^(A_2)*exp.(-(log10.(f./f2)).^2/Delta2^2)
end

function Omega_broken_powerlaw(A_1,f,gamma1,gamma2,ftt)
    return 10^(A_1) .* (f .> ftt) .* ((f ./ ftt) .^ gamma1) .+ (f .<= ftt) .* ((f ./ ftt) .^ gamma2)
end

function Omega_three_peaks(A_1,A_2,A_3,f,f1,f2,f3,Delta1,Delta2,Delta3)
    return 10^(A_1) * exp.(-(log10.(f ./ f1)).^2 / Delta1^2) .+
    10^(A_2) * exp.(-(log10.(f ./ f2)).^2 / Delta2^2) .+
    10^(A_3) * exp.(-(log10.(f ./ f3)).^2 / Delta3^2)
end