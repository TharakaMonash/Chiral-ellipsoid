function [AbsL,AbsR] = calcAbsN2(eMedium,lambda,eEllipsoid,longRadius,transRadius,chirality,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the absorption of N longitudinal ellipsoids based on the
% following parameters.
% Parameters:
%   eMedium     : Permittivity of the background medium
%   lambda      : wavelength as an array (in nm)
%   eEllipsoid  : Permittivity of the ellipsoid medium
%   longRadius  : Longitudinal radius of the ellipsoid
%   transRadius : Transverse radius of the ellipsoid
%   chirality   : Chirality parameter of the ellipsoid
%   N           : Number of ellipsoids
% Returns:
%   AbsL         : Absorption of N ellipsoids as when excited with left
%   polarized light
%   AbsR         : Absorption of N ellipsoids as an array with right
%   polarized light
% Author: Tharaka Perera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Absorption of a longitudinal ellipsoid with left polarized light   
    AbsLlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,-1,1/sqrt(2),1);
    % Absorption of a longitudinal ellipsoid with right polarized light  
    AbsRlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,1, 1/sqrt(2),1);
    % Absorption of a transverse ellipsoid with left polarized light
    %AbsLtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,-1,1/sqrt(2),0);
    % Absorption of a transverse ellipsoid with right polarized light
    %AbsRtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transRadius, chirality,1, 1/sqrt(2),0);
    
    % Absorption of N ellipsoids with left polarized light   
    AbsL =N*(AbsLlong);
    % Absorption of N ellipsoids with right polarized light   
    AbsR =N*(AbsRlong);
end