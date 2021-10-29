function [AbsL,AbsR] = calcAbsNabsorbing(eMedium,lambda,eEllipsoid,longRadius,transverseRadius,chirality,N)
    AbsLlong = calcAbsAbsorbingMedium(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1,-1,1/sqrt(2),1);
    AbsRlong = calcAbsAbsorbingMedium(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1,1, 1/sqrt(2),1);
    AbsLtrans = calcAbsAbsorbingMedium(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1,-1,1/sqrt(2),0);
    AbsRtrans = calcAbsAbsorbingMedium(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1,1, 1/sqrt(2),0);
    
    AbsL =N*(AbsLlong+AbsLlong+AbsLtrans)/3;
    AbsR =N*(AbsRlong+AbsRlong+AbsRtrans)/3;
end