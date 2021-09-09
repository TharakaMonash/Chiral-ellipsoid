%Chirality and Permitivitty
function [xi_c,e_c] = calcChiralParam(lambda)
    tau = 3.5248e-20;%3.2044e-20;;  %in J
    c_const = 299792458;
    hbar =1.0545718e-34;
    %parameters
    %beeta_c = 5e-22;
    beeta_c = 1e-22;
    %beeta_c = 2.5e-24;%Final one that matches comsol
    %beeta_c = 5e-22; % one that matches experimental results
    %beeta_c = 1e-23;   %one used in comsol
    gamma_c = 5e-25; %one used in comsol
    %gamma_c = 1e-17;

    lambda0 = 625e-9;%585e-9; %635e-9;
    w0 = hbar*2*pi*c_const/lambda0;
    f = c_const.*1e9./lambda;

    xi_c = beeta_c*((1./(hbar*2*pi*f-w0+1i*tau))+(1./(hbar*2*pi*f+w0+1i*tau)));
    e_c = 1.5 - gamma_c*((1./(hbar*2*pi*f-w0+1i*tau))-(1./(hbar*2*pi*f+w0+1i*tau)));
end