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
AR = 0.66;
Width = 10;
e_w1 = 2.5;
e_w2 = 4.5;
% Change these two parameters to get parameters for a specific range
targetLambda1 = 800;
targetLambda2 = 1000;
% Outer dimensions of chiral Nanorod
Length = Width/AR;

target1=[targetLambda1,targetLambda1];
target2=[targetLambda2,targetLambda2];
plot(target1, [0,1.5],'--r','DisplayName','Target lower');
hold on;
plot(target2, [0,1.5],'--b','DisplayName','Target higher');

% Chirality parameter
[CL,e_c] = calcChiralParam(lambda);
%maxwell garnett to get effective medium parameters
frac =0.05;
e_eff = e_Au .*(2*frac*(-e_Au+e_c)+e_c+2*e_Au)./(e_c+2*e_Au-frac*(e_c-e_Au));
CL_modified = 3*frac*(CL.*e_Au./(e_c+2*e_Au-frac*(e_c-e_Au)));

linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};

it=1;
mserror=2;
jump = 0.5;

t = 1000;
cooling_factor = 0.001;
pertubs_for_annealing=100;
while abs(t) > 0.00001 
    for i=1:1:pertubs_for_annealing
        % current value
        Length = Width/AR;
        [AbsL1,AbsR1]= calcAbsN2(e_w1, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids);
        [AbsL2,AbsR2]= calcAbsN2(e_w2, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids);
        [maxval,maxlambda1] = max(AbsL1);
        [maxva2,maxlambda2] = max(AbsL2);
        
        perturbed_AR = rand;
        Length = Width/perturbed_AR;
        [PAbsL1,PAbsR1]= calcAbsN2(e_w1, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids);
        [PAbsL2,PAbsR2]= calcAbsN2(e_w2, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids);
        [Pmaxval,Pmaxlambda1] = max(PAbsL1);
        [Pmaxva2,Pmaxlambda2] = max(PAbsL2);
        %myerror = (abs(lambda(maxlambda1)-targetLambda1)+abs(lambda(maxlambda2)-targetLambda2));
        mserror = abs(lambda(Pmaxlambda1)-targetLambda1)+abs(lambda(Pmaxlambda2)-targetLambda2)-(abs(lambda(maxlambda1)-targetLambda1)+abs(lambda(maxlambda2)-targetLambda2));
        %errors(it)=myerror; % plot this to see the errors and convergence
        %temp(it)=t;
        exp(-mserror/t)
        if mserror<0
            AR = perturbed_AR;
        elseif  rand < exp(-mserror/t)
            AR = perturbed_AR;
        end
        if it==1
            plot(lambda, AbsL1/max(AbsL1),':r','DisplayName','Initial');
            hold on;
            plot(lambda, AbsL2/max(AbsL2),':b','DisplayName','Initial');
            hold on;
        end
        it=it+1;
    end
    t = cooling_factor*t;
end

plot(lambda, AbsL1/max(AbsL1),'r','DisplayName','Final');
plot(lambda, AbsL2/max(AbsL2),'b','DisplayName','Final');
legend()
leg = legend('Location','NorthWest','NumColumns',3);
ylim([0,1.4]);
leg.ItemTokenSize = [20,30];
xlabel('Wavelength (nm)');
ylabel('Absorption (arb. units)');