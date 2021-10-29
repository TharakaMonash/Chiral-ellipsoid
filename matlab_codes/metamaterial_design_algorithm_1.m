close all;
clear all;

lambda = 400:1:1200;
% Au permittivity
for j = 1:length(lambda)
    [eps1(j), eps2(j)] =  getEpsAuByLambda(lambda(j), 10e3);
end
e_Au = eps1 + 1i*eps2;

% Parameter initialisation
NoOdEllipsoids=1e12;
AR = 0.5;
Width = 10;
e_w = 1.75;
targetLambda = 800;
% Outer dimensions of chiral Nanorod
Length = Width/AR;

target=lorentz(lambda, targetLambda, 50, 1)
plot(lambda, target,'--k','DisplayName','Target');
hold on;
% Chirality parameter
[CL,e_c] = calcChiralParam(lambda);
%maxwell garnett to get effective medium parameters
frac =0.05;
e_eff = e_Au .*(2*frac*(-e_Au+e_c)+e_c+2*e_Au)./(e_c+2*e_Au-frac*(e_c-e_Au));
CL_modified = 3*frac*(CL.*e_Au./(e_c+2*e_Au-frac*(e_c-e_Au)));

linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};
prev=0;
it=1;
mserror=2;
jump = 0.5;
while abs(mserror) >1
    Length = Width/AR;
    [AbsL,AbsR]= calcAbsN2(e_w, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids); % longitudinal esponse only
    [maxval,maxlambda] = max(AbsL);
    mserror=lambda(maxlambda)-targetLambda;
    prev = maxlambda;
    if it==1
        plot(lambda, AbsL/max(AbsL),':b','DisplayName','Initial guess');
        hold on;
    end
    jump=jump/2;
    if mserror > 0
       AR=AR+jump; 
    else 
       AR=AR-jump; 
    end
    it=it+1;
end

plot(lambda, AbsL/max(AbsL),'r','DisplayName','Final');
legend()
leg = legend('Location','NorthEast','NumColumns',1);
leg.ItemTokenSize = [20,30];
xlabel('Wavelength (nm)');
ylabel('Absorption (arb. units)');