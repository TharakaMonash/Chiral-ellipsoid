%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This template file can be used to calculate the absorption and cd spectra of chiral plasmonic ellipsoids based on the following parameters.
% Parameters:
%   eMedium     : Permittivity of the background medium
%   lambda      : wavelength as an array (in nm)
%   eEff        : Permittivity of the ellipsoid medium
%   Length      : length of the longitudinal axis of the ellipsoid
%   Width       : length of the Transverse axis of the ellipsoid
%   xiModified  : Chirality parameter of the ellipsoid
%   
% 
% Author: Tharaka Perera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

lambda = 400:1:1000;
%% Styling plotlines
lineStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};

%% Material property calculation
% Au permittivity
for j = 1:length(lambda)
    [epsAu1(j), epsAu2(j)] =  getEpsAuByLambda(lambda(j), 10e3);
end
eAu = epsAu1 + 1i*epsAu2;
% Ag permittivity
for j = 1:length(lambda)
    [epsAg1(j), epsAg2(j)] =  getEpsAgByLambda(lambda(j), 10e3);
end
eAg = epsAg1 + 1i*epsAg2;
% Chiral material
[xi,eC] = calcChiralParam(lambda);

%Maxwell garnett to get effective medium parameters of the ellipsoid
frac =0.05;
eEff = eAu.*(2*frac*(-eAu+eC)+eC+2*eAu)./(eC+2*eAu-frac*(eC-eAu));
xiModified = 3*frac*(xi.*eAu./(eC+2*eAu-frac*(eC-eAu)));

%% Dimensions of chiral plasmonic ellipsoid
Length = 30; 
Width = 10;
NoOfEllipsoids=1e12;

%% Results changing background permittivity 
%Background permittivity
eMedium =1.75:0.15:2.5;
% Array initialisation
AbsLAll=zeros(length(eMedium),length(lambda));
CDall=zeros(length(eMedium),length(lambda));

for count=1:1:length(eMedium)
    em = eMedium(count);
    [AbsL,AbsR]= calcAbsN(em, lambda, eEff, Length/2, Width/2, xiModified,NoOfEllipsoids);
    AbsLAll(count,:)=AbsL;
    CD = AbsR-AbsL;
    CDall(count,:)=CD;
end
%Calculate the scales
ABSscale=max(max(AbsLAll));
CDscale=max(max(CDall));
% Abs with medium permittivity plot
figure
set(gcf,'units','inches','position',[1,2,8,4]);
t = tiledlayout(1,2); % Requires R2019b or later
t.TileSpacing = 'none';
t.Padding = 'compact';
title(t,'Changes in chiroptical properties with background permittivity')

nexttile
for count=1:1:length(eMedium)
    plot(lambda, AbsLAll(count,:)/ABSscale,lineStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(eMedium(count))));
    hold on
end
legend('Location','northwest');
legend show
xlabel('Wavelength (nm)');
ylabel('Absorption (arb. units)');
xlim([400,800]);
% CD with medium permittivity plot
nexttile
for count=1:1:length(eMedium)
    plot(lambda, CDall(count,:)/CDscale,lineStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(eMedium(count))));
    hold on
end
legend('Location','northwest');
legend show
leg.ItemTokenSize = [16,10];
xlabel('Wavelength (nm)');
ylabel('CD (arb. units)');
xlim([400,800]);

%%  Results changing orientation
clear AbsLAll,CDall;
eMedium = 1.76;
anglestep =15;
h = zeros(1,2*length(0:anglestep:90));
% Array initialisation
AbsLAll = zeros(length(90:-anglestep:0),length(lambda));
CDall = zeros(length(90:-anglestep:0),length(lambda));

count =0;
for theta=90:-anglestep:0
    count = count +1;
    [AbsL,AbsR] = calcAbsOriented(eMedium, lambda, eEff, Length/2, Width/2, xiModified,theta);
    CD = AbsR-AbsL;
    AbsLAll(count,:)=AbsL;
    CDall(count,:)=CD;
end
Absscale=max(max(AbsLAll));
CDscale=max(max(abs(CDall)));

figure
set(gcf,'units','inches','position',[4,2,8,4]);
t = tiledlayout(1,2); % Requires R2019b or later
t.TileSpacing = 'none';
t.Padding = 'compact';
title(t,'Changes in chiroptical properties with orientation')
% Abs with orientation
nexttile
count =0;
for theta=90:-anglestep:0 
    count = count +1;
    plot(lambda, AbsLAll(count,:)/Absscale,lineStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat("\theta=",num2str(theta),char(176)));
    hold on
end 
legend('Location','northwest');
legend show
leg.ItemTokenSize = [16,10];
xlabel('Wavelength (nm)');
ylabel('Absorption (arb. units)');
xlim([400,800]);
% CD with orientation
nexttile
count =0;
for theta=90:-anglestep:0 
    count = count +1;
    plot(lambda, CDall(count,:)/CDscale,lineStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat("\theta=",num2str(theta),char(176)));
    hold on
end 
legend('Location','northwest');
legend show
leg.ItemTokenSize = [16,10];
xlabel('Wavelength (nm)');
ylabel('CD (arb. units)');
xlim([400,800]);
%% Results changing geometric factor

%Ellipsoid sizes
width_array = 5:1:19;
length_array = 30*(ones(1,length(width_array)));

long_peakL = zeros(1,length(length_array));
long_peakR = zeros(1,length(length_array));
long_peakCD = zeros(1,length(length_array));

for count=1:1:length(length_array)
    Width =width_array(count);
    Length=length_array(count);
    
    [AbsL,AbsR]=calcAbsOriented(eMedium, lambda, eEff, Length/2, Width/2, xiModified, NoOfEllipsoids);
    [maxval,index]=max(AbsL); 
    long_peakL(count)=lambda(index);
    [maxval,index]=max(AbsR); 
    long_peakR(count)=lambda(index);
    
    CD = AbsR-AbsL;
    [maxval,index]=max(abs(CD)); 
    long_peakCD(count)=lambda(index);

end

f = figure;
a1 = axes;

plot(width_array,long_peakL,'-ok','LineWidth',1,'DisplayName', 'Absorption' );
hold on
plot(width_array,long_peakCD,'-sr','LineWidth',1,'DisplayName', 'CD');
leg = legend('Location','northeast');
legend show
set(gcf,'units','inches','position',[6,2,5,4]);

set(gca,'units','inches','position',[0.75,0.5,4,3]);
ylim([500,1000]);

a1.Box="on";
set(gca,'Linewidth',1);
% Adding an extra axis
xt = get(a1,'XTick');
a2 = copyobj(a1,f);
set(a2,'Color','none');
set(a2,'Ytick',[]);
set(a2,'XAxisLocation','top');
set(a2,'Box','off');
delete(a2.Children);

% Convert the xTick values based on a different scale
R = xt/30;
e = sqrt(1-R.^2);
L = ((1-e.^2)/e.^2) .* ((1./(2*e)).*log((1+e)./(1-e)) - 1);
converted_values = round(L,2);
set(a2,'XTickLabel',converted_values)
xlabel(a2,'Geometric factor')
xlabel(a1,'Width (nm)');
ylabel(a1,'Max Wavelength (nm)');

