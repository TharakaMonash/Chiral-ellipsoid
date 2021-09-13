close all;
clear all;

update_files=true;
linStyles = {'-','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};
%% Analylitical Results
lambda = 200:1:800;
% Au permittivity
for j = 1:length(lambda)
    [eps1(j), eps2(j)] =  getEpsAuByLambda(lambda(j), 10e3);
end
e_Au = eps1 + 1i*eps2;

for j = 1:length(lambda)
    [eps1(j), eps2(j)] =  getEpsAgByLambda(lambda(j), 10e3);
end
e_Ag = eps1 + 1i*eps2;

%Background permittivity
e_w_array=[1.76,1.86,1.96,2.07,2.16];
e_w = 1.76;
% Outer dimensions of chiral Nanorod
PlotLegends = {
'Ag/Au=0,    R_l=29,    R_t=9.5'
'Ag/Au=0.5, R_l=29.5,  R_t=11'
'Ag/Au=1,    R_l=30,    R_t=13'
'Ag/Au=1.5, R_l=30,    R_t=15'
'Ag/Au=2,    R_l=30.5,  R_t=18'
    };
L=[58,59,60,60,61];
W=[19,22,26,30,36];
%L=[60,60,60,60,60];
%W=[26,28,30,32,34];
%Length = 62;
%Width = 29;
NoOdEllipsoids=1e12;
% Chirality parameter
[CL,e_c] = calcChiralParamCysteine(lambda);
%Cysteine volume fraction
frac =0.05;
% Ag volume fraction.
fracs =[0,0.333,0.5,0.6,0.677];
%fracs =[0.5,0.5,0.5,0.5,0.5];
for count=1:1:length(fracs)
    e_w=e_w_array(count);
    Length=L(count);
    Width=W(count);
    frac2=fracs(count);
    e_eff =e_Au.*(2*frac2*(-e_Au+e_Ag)+e_Ag+2*e_Au)./(e_Ag+2*e_Au-frac2*(e_Ag-e_Au));
    CL_modified = 3*frac*(CL.*e_eff./(e_c+2*e_eff-frac*(e_c-e_eff)));
    e_eff2 = e_eff .*(2*frac*(-e_eff+e_c)+e_c+2*e_eff)./(e_c+2*e_eff-frac*(e_c-e_eff));
    
    [AbsL,AbsR] = calcAbsN(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,NoOdEllipsoids); 
    %cdlong = calcCD(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,-1,1/sqrt(2),1);
    %cdtrans = calcCD(e_w, lambda, e_eff2, Length/2, Width/2,  CL_modified,-1,1/sqrt(2),0);

    %CD = (2*cdlong+cdtrans)/3;
    CD=AbsL-AbsR;
    CDall(count,:)=CD;
    Absall(count,:)=AbsL;
end

ABSscale = max(max(Absall));
CDscale = max(max(abs(CD)));


%For cysteine molecules
% e_eff =e_c;
% CL_modified = 3*frac*(CL.*e_eff./(e_c+2*e_eff-frac*(e_c-e_eff)));
% e_eff2 = e_eff .*(2*frac*(-e_eff+e_c)+e_c+2*e_eff)./(e_c+2*e_eff-frac*(e_c-e_eff));
% 
% Abs = (calcAbs(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,-1,1/sqrt(2),1)+calcAbs(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,-1,1/sqrt(2),0))/2;
% AbsR = (calcAbs(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,1,1/sqrt(2),1)+calcAbs(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,1,1/sqrt(2),0))/2;
% cdlong = calcCD(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,-1,1/sqrt(2),1);
% cdtrans = calcCD(e_w, lambda, e_eff2, Length/2, Width/2,  CL_modified,-1,1/sqrt(2),0);
% 
% 
% CD = (cdlong+cdtrans)/2;

linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};
linStyles2 = {':o',':+',':*',':<',':s',':d',':^',':v',':>',':p',':x',':h'};

fracs_2 = [0,0.5,1,1.5,2];

figure
for i=1:1:length(fracs)
    plot(lambda, Absall(i,:)/ABSscale,linStyles{i},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    %plot(lambda, Absall(i,:)/ABSscale,linStyles{i},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',PlotLegends{i});
    hold on
    xlim([300,800])
    ylim([0,1.4])
end
% plot(lambda, Abs/ABSscale ,'--k','LineWidth',1,'DisplayName',"Molecules")
set(gcf,'units','inches','position',[0.5,0.5,3.33,2])
%set(gca,'units','inches','position',[0.325,0.3,2.85,2.1])%size 8 font
set(gca,'units','inches','position',[0.35,0.35,2.85,1.6])%size 9 font
leg = legend('Location','northwest','NumColumns',3);
% legend('Location','northwest','NumColumns',2);
legend show
leg.ItemTokenSize = [16,10];

ylabel("Absorption (arb. units)")
xlabel("Wavelength (nm)")
box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

if (update_files)
    saveas(gcf,'Fig_6_abs_au_ag_new','epsc')
    %print('abs_au_ag','-depsc', '-opengl', '-r600');
end

figure
set(gcf,'units','inches','position',[0.5,0.5,3.33,2.75])
t = tiledlayout(5,1); % Requires R2019b or later
t.TileSpacing = 'none';
t.Padding = 'compact';

for i=1:1:length(fracs)
    nexttile
    plot(lambda, -CDall(i,:)/CDscale,'-b','LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    hold on
    plot(lambda, CDall(i,:)/CDscale,':r','LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    if i ~= 5
        set(gca,'xtick',[])
    end
    ylim([-1.2,1.2])
    xlim([300,800])
    if i==3
        ylabel("CD (arb. units)")
    end 
    %str = strcat('Ag/Au =  ',num2str(fracs_2(i)));
    str=PlotLegends{i};
    text(500,0.6,str)
    set(gca, 'FontName', 'Arial')
    set(gca,'FontSize',8)
    box on
    set(gca,'Linewidth',1)
end
% nexttile
% 
% plot(lambda, -CD/CDscale,'k','LineWidth',1,'DisplayName',"Molecules")
% hold on
% plot(lambda, CD/CDscale,':k','LineWidth',1,'DisplayName',"Molecules")
% ylim([-1.2,1.2])
% xlabel("Wavelength (nm)")
% text(650,0.5,"Molecules")
% xlim([300,800])
box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

if (update_files)
    saveas(gcf,'Fig_7_cd_au_ag_new','epsc')
    %print('cd_au_ag','-depsc', '-opengl', '-r600');
end
