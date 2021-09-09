function [AbsL,AbsR] = calcAbsOriented(eMedium,lambda,eEllipsoid,longRadius,transRadius,chirality,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the absorption of a oriented ellipsoid based on the
% following parameters.
% Parameters:
%   eMedium     : Permittivity of the background medium
%   lambda      : wavelength as an array (in nm)
%   eEllipsoid  : Permittivity of the ellipsoid medium
%   longRadius  : Longitudinal radius of the ellipsoid
%   transRadius : Transverse radius of the ellipsoid
%   chirality   : Chirality parameter of the ellipsoid
%   theta       : angle between the electric field propagation direction
%   and the longitudinal axis of the ellipsoid
% Returns:
%   AbsL         : Absorption of an oriented ellipsoid as when excited with left polarized light
%   AbsR         : Absorption of an oriented ellipsoid as when excited with right polarized light
% Author: Tharaka Perera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Absorption of a longitudinal ellipsoid with left polarized light   
    AbsLlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,-1,1/sqrt(2),1);
    % Absorption of a longitudinal ellipsoid with right polarized light  
    AbsRlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,1, 1/sqrt(2),1);
    % Absorption of a transverse ellipsoid with left polarized light
    AbsLtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,-1,1/sqrt(2),0);
    % Absorption of a transverse ellipsoid with right polarized light
    AbsRtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,1, 1/sqrt(2),0);
    % Absorption of a theta oriented ellipsoid with left polarized light 
    AbsL = AbsLlong*(sind(theta))^2+AbsLtrans*(cosd(theta))^2;
    % Absorption of a theta oriented ellipsoid with right polarized light 
    AbsR = AbsRlong*(sind(theta))^2+AbsRtrans*(cosd(theta))^2;
end