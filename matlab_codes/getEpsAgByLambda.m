%lambda_nm should be between 188nm and 1937nm
function [epsRe, epsIm] = getEpsAgByLambda(lambda_nm, numInterpolationPts)
    refIndexMat = load('Ag_JC1972_lambdaMum_n_k.mat');
    refIndexMat = cell2mat(struct2cell(refIndexMat));
    
    lambda_interp = linspace(0.1879, 1.937, numInterpolationPts); %lambda in um
    n_interp = spline(refIndexMat(:,1), refIndexMat(:,2), lambda_interp);
    k_interp = spline(refIndexMat(:,1), refIndexMat(:,3), lambda_interp);
    
    [~, idx] = min(abs(lambda_interp-lambda_nm/1000));
    n = n_interp(idx);
    k = k_interp(idx);
    
    epsRe = n^2-k^2;
    epsIm = 2*n*k;
end

