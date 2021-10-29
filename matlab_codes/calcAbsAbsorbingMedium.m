function absorption = calcAbsAbsorbingMedium(eMedium,lambda,eEllipsoid,longRadius,transverseRadius,chirality,Cx,Cy,K,orientation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the absorption of a single ellipsoid in an absorbing medium based on the
% following parameters.
% Parameters:
%   eMedium     : Permittivity of the background medium (Complex)
%   lambda      : wavelength as an array (in nm)
%   eEllipsoid  : Permittivity of the ellipsoid medium (Complex)
%   longRadius  : Longitudinal radius of the ellipsoid
%   transRadius : Transverse radius of the ellipsoid
%   chirality   : Chirality parameter of the ellipsoid
%   C           : Parameter that describe the electric field
%   K           : Parameter that describe the electric field
%   orientation : 1 for longitudinal ellipsoid, 0 for transverse ellipsoid.
% Returns:
%   abs         : Absorption of the ellipsoid as an array
% Author: Tharaka Perera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = transverseRadius/longRadius;
    E0=1;
    e = sqrt(1-R^2);
    L = ((1-e^2)/e^2) * ((1/(2*e))*log((1+e)/(1-e)) - 1);
    N = 1;
    SL= 299792458;
    V= 4/3*pi*longRadius*transverseRadius*transverseRadius*(10^(-9))^3;
    e1 = real(eEllipsoid);
    e2 = imag(eEllipsoid);
    
    xi1 = real(chirality);
    xi2 = imag(chirality);
    
    if orientation == 1%"long"
        P = [L, (1-L)/2, (1-L)/2];
    elseif orientation == 0%"trans"
        P = [(1-L)/2, (1-L)/2, L];
    end
    temp = 0;
    for i = 1:2
        if i==1
            Cs=Cx;
            Ct=Cy;
        elseif i==2
            Cs=Cy;
            Ct=Cx;
        end
        temp3 = ((-Cs*eEllipsoid*eMedium+Ct*sqrt(eMedium)*(2*P(i)*eEllipsoid+eMedium-P(i)*eMedium).*chirality)./(P(i)*conj(eEllipsoid) + conj(eMedium)-P(i)*conj(eMedium))) .*(conj(Cs*eMedium-Ct*P(i)*sqrt(eMedium).*chirality));
        temp2 =((Cs*sqrt(eMedium)-Ct*P(i)*chirality)*abs(eMedium)*Ct.*chirality + temp3)./(P(i)*(eEllipsoid-eMedium)+eMedium);
        temp = temp+ imag(temp2);
        disp(temp)
    end
    absorption = -V*E0^2*K^2*temp*SL./lambda;
end






