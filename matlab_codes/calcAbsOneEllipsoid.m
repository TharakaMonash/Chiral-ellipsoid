function abs = calcAbsOneEllipsoid(eMedium,lambda,eAu,longRadius,transverseRadius,chirality,C,K,orientation)
    R = transverseRadius/longRadius;
    e = sqrt(1-R^2);
    L = ((1-e^2)/e^2) * ((1/(2*e))*log((1+e)/(1-e)) - 1);
    N = 1;
    SL= 299792458;
    V= 4/3*pi*longRadius*transverseRadius*transverseRadius*(10^(-9))^3;
    e1 = real(eAu);
    e2 = imag(eAu);
    
    xi1 = real(chirality);
    xi2 = imag(chirality);
    
    if orientation == 1%"long"
        P = [L, (1-L)/2, (1-L)/2];
    elseif orientation == 0%"trans"
        P = [(1-L)/2, (1-L)/2, L];
    end
    temp = 0;
    for i = 1:2
        %temp = temp + (sqrt(eMedium)*e2 - 2*C*((P(i)-1)*eMedium-P(i)*(e1+e2)).*xi1)./(((P(i)-1)*eMedium-P(i)*e1).^2+(P(i)*e2).^2);
        %temp = temp + (sqrt(eMedium)*e2 + 2*C*((1-P(i))*eMedium+P(i)*(e1+e2)).*xi1)./(((P(i)-1)*eMedium-P(i)*e1).^2+(P(i)*e2).^2);
        temp = temp + (sqrt(eMedium)*e2 - 2*C*((eMedium+P(i).*(e1-eMedium)).*xi2+P(i).*e2.*xi1))./(((P(i)-1)*eMedium-P(i)*e1).^2+(P(i)*e2).^2);
    end
    abs = N*V*K^2*SL*(eMedium^(3/2))*temp./lambda;
end