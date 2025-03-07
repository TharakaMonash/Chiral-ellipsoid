function [epsRe, epsIm] = getEpsAuByLambda(lambda_nm, numInterpolationPts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the permittivity of gold based on refractive
% indices from [P. B. Johnson and R.-W. Christy, �Optical constants of the noble
% metals,� Phys. review B6,4370 (1972).]
% IMPORTANT: lambda_nm should be between 188nm and 1937nm
% following parameters.
% Parameters:
%   lambda_nm           : wavelength in nm as an array
%   numInterpolationPts : number of interpolation points
% Returns:
%   epsRe         : Real part of permittivity of silver
%   epsIm         : Imaginary part of permittivity of silver
% Author: Sudaraka Mallwaarachchi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % JC gold refractive indices
    refIndexMat = load('Au_JC1972_lambdaMum_n_k.mat');
    refIndexMat = cell2mat(struct2cell(refIndexMat));
    
    lambda_interp = linspace(0.1879, 1.937, numInterpolationPts);
    n_interp = spline(refIndexMat(:,1), refIndexMat(:,2), lambda_interp);
    k_interp = spline(refIndexMat(:,1), refIndexMat(:,3), lambda_interp);
    
    [~, idx] = min(abs(lambda_interp-lambda_nm/1000));
    n = n_interp(idx);
    k = k_interp(idx);
    
    epsRe = n^2-k^2;
    epsIm = 2*n*k;
end

