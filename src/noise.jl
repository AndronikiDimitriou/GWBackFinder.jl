"""
Define functions to generate noise according to https://arxiv.org/pdf/2009.11845.pdf
"""
function noise_Pims(f,P,c)
    return P^2*(10^(-12))^2*(1+(2*10^(-3)/f)^4)*(2*pi*f/c)^2
end

function noise_Pacc(f,A,c)
  return A^2*(10^(-15))^2*(1+(0.4*10^(-3)/f)^2)*(1+(f*10^3/8)^4)*(1/(2*pi*f))^4*(2*pi*f/c)^2
end

function Omega_noiseh2_TT(f,L,c,P,A)
    return 4*pi^2*f^3*Sn_TT(f,L,c,P,A)/(3*(3.24*10^(-18))^2)
end

function Omega_noiseh2_AA(f,L,c,P,A)
    return 4*pi^2*f^3*Sn_AA(f,L,c,P,A)/(3*(3.24*10^(-18))^2)
end

function Sn_TT(f,L,c,P,A)
    return (N_TT(f,L,c,P,A)/R_TT(f,L,c))
end

function Sn_AA(f,L,c,P,A)
    return (N_AA(f,P,A,L,c)/R_AA(f,L,c))
end

function R_TT(f,L,c)
    return 16*(sin(2*pi*f*L/c))^2*(2*pi*f*L/c)^2*R_tilde_TT(f,L,c)
end

function R_AA(f,L,c)
    return 16*(sin(2*pi*f*L/c))^2*(2*pi*f*L/c)^2*R_tilde_AA(f,L,c)
end

function R_tilde_AA(f,L,c)
    return 9/20*1/(1+0.7*(2*pi*f*L/c)^2)
end

function R_tilde_TT(f,L,c)
    return 9/20*(2*pi*f*L/c)^6/(1.8*10^3+0.7*(2*pi*f*L/c)^8)
end

function N_TT(f,L,c,P,A)
    return 16*(sin(2*pi*f*L/c))^2*(2*(1-(cos(2*pi*f*L/c)))^2*noise_Pacc(f,A,c)+(1-cos(2*pi*f*L/c))*noise_Pims(f,P,c))
end

function N_AA(f,P,A,L,c)
    return 8*sin(2*pi*f*L/c)^2*(4*(1+cos(2*pi*f*L/c)+cos(2*pi*f*L/c)^2)*noise_Pacc(f,A,c)+(2+cos(2*pi*f*L/c))*noise_Pims(f,P,c))
end