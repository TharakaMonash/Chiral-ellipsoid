function [AbsL,AbsR] = calcAbsN(eMedium,lambda,eEllipsoid,longRadius,transverseRadius,chirality,N)
    AbsLlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,-1,1/sqrt(2),1);
    AbsRlong = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1, 1/sqrt(2),1);
    AbsLtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,-1,1/sqrt(2),0);
    AbsRtrans = calcAbsOneEllipsoid(eMedium, lambda, eEllipsoid, longRadius, transverseRadius, chirality,1, 1/sqrt(2),0);
    
    AbsL =N*(AbsLlong+AbsLlong+AbsLtrans)/3;
    AbsR =N*(AbsRlong+AbsLlong+AbsRtrans)/3;
end