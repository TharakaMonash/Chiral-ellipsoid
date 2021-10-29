close all;
clear all;
Folderpath = "..\experimental_data\";
Filenames = ["cd_0.csv","cd_0.5.csv","cd_1.csv","cd_1.5.csv","cd_2.csv"];

colors=[         
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


for i=1:length(Filenames)
    Filepath = strcat(Folderpath,Filenames(i));
    ExpDataTable = readtable(Filepath);
    if i==1
        ExpDataArray0 = table2array(ExpDataTable);
    elseif i==2
        ExpDataArray05 = table2array(ExpDataTable);
    elseif i==3
        ExpDataArray10 = table2array(ExpDataTable);
    elseif i==4
        ExpDataArray15 = table2array(ExpDataTable);
    else
        ExpDataArray2 = table2array(ExpDataTable);
    end
end
MAXEXP =60;


Folderpath = "..\experimental_data\";
Filenames = ["abs_0.csv","abs_0.5.csv","abs_1.0.csv","abs_1.5.csv","abs_2.csv"];

for i=1:length(Filenames)
    Filepath = strcat(Folderpath,Filenames(i));
    ExpAbsDataTable = readtable(Filepath);
    if i==1
        ExpAbsDataArray0 = table2array(ExpAbsDataTable);
        ExpAbsDataArray0(:,2) = ExpAbsDataArray0(:,2)/max(ExpAbsDataArray0(:,2));
    elseif i==2
        ExpAbsDataArray05 = table2array(ExpAbsDataTable);
        ExpAbsDataArray05(:,2) = ExpAbsDataArray05(:,2)/max(ExpAbsDataArray05(:,2));
    elseif i==3
        ExpAbsDataArray10 = table2array(ExpAbsDataTable);
        ExpAbsDataArray10(:,2) = ExpAbsDataArray10(:,2)/max(ExpAbsDataArray10(:,2));
    elseif i==4
        ExpAbsDataArray15 = table2array(ExpAbsDataTable);
        ExpAbsDataArray15(:,2) = ExpAbsDataArray15(:,2)/max(ExpAbsDataArray15(:,2));
    else
        ExpAbsDataArray2 = table2array(ExpAbsDataTable);
        ExpAbsDataArray2(:,2) = ExpAbsDataArray2(:,2)/max(ExpAbsDataArray2(:,2));
    end
end
MAXEXPABS =1;

expcol={':b',':r'}
expcolabs={':k',':r'}

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


e_w_array=[2.075,2.075,1.85,1.85,2];

PlotLegends = {
'Ag/Au=0,    R_l=29,    R_t=9.5'
'Ag/Au=0.5, R_l=29.5,  R_t=11'
'Ag/Au=1,    R_l=30,    R_t=13'
'Ag/Au=1.5, R_l=30,    R_t=15'
'Ag/Au=2,    R_l=30.5,  R_t=18'
    };
L=[58,59,60,60,61];
W=[19,22,26,30,36];

NoOdEllipsoids=1e12;
% Chirality parameter
[CL,e_c] = calcChiralParamCysteine(lambda);
%Cysteine volume fraction
frac =0.05;
% Ag volume fraction.
fracs =[0,0.333,0.5,0.6,0.677];

for count=1:1:length(fracs)
    e_w=e_w_array(count);
    Length=L(count);
    Width=W(count);
    frac2=fracs(count);
    e_eff =e_Au.*(2*frac2*(-e_Au+e_Ag)+e_Ag+2*e_Au)./(e_Ag+2*e_Au-frac2*(e_Ag-e_Au));
    CL_modified = 3*frac*(CL.*e_eff./(e_c+2*e_eff-frac*(e_c-e_eff)));

    e_eff2 = e_eff .*(2*frac*(-e_eff+e_c)+e_c+2*e_eff)./(e_c+2*e_eff-frac*(e_c-e_eff));
    
    [AbsL,AbsR] = calcAbsN(e_w, lambda, e_eff2, Length/2, Width/2, CL_modified,NoOdEllipsoids); 

    CD=AbsL-AbsR;
    CDall(count,:)=CD;
    Absall(count,:)=AbsL;
end

ABSscale = max(max(Absall));
CDscale = max(max(abs(CD)));

linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};
linStyles2 = {':o',':+',':*',':<',':s',':d',':^',':v',':>',':p',':x',':h'};

fracs_2 = [0,0.5,1,1.5,2];

figure
set(gcf,'units','inches','position',[0.5,0.5,6.5,2.75])


subplot(3,2,1)
for i=1:1:length(fracs)
    plot(lambda, Absall(i,:)/max(Absall(i,:)),linStyles{i},'color',colors(i,:),'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    hold on
    xlim([300,800])
    ylim([0,1.05])
    
    if i==1
        plot(ExpAbsDataArray0(:,1),ExpAbsDataArray0(:,2)/MAXEXPABS,':','color',colors(i,:),'LineWidth',1,'HandleVisibility','off')
    elseif i==2
        plot(ExpAbsDataArray05(:,1),ExpAbsDataArray05(:,2)/MAXEXPABS,':','color',colors(i,:),'LineWidth',1,'HandleVisibility','off')
    elseif i==3
        plot(ExpAbsDataArray10(:,1),ExpAbsDataArray10(:,2)/MAXEXPABS,':','color',colors(i,:),'LineWidth',1,'HandleVisibility','off')
    elseif i==4
        plot(ExpAbsDataArray15(:,1),ExpAbsDataArray15(:,2)/MAXEXPABS,':','color',colors(i,:),'LineWidth',1,'HandleVisibility','off')
    else
        plot(ExpAbsDataArray2(:,1),ExpAbsDataArray2(:,2)/MAXEXPABS,':','color',colors(i,:),'LineWidth',1,'HandleVisibility','off')
    end
end

leg = legend('units','inches','Position',[0.375, 2.2, 2.85, 0.5],'NumColumns',3)
leg.Position=[0.375, 2.2, 2.85, 0.5];

legend show
leg.ItemTokenSize = [16,10];

ylabel("Absorption (arb. units)")
xlabel("Wavelength (nm)")
box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')

for i=1:1:length(fracs)
    subplot(3,2,i+1)
    plot(lambda, -CDall(i,:)/CDscale,'-b','LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    hold on
    plot(lambda, CDall(i,:)/CDscale,'--r','LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat('Ag/Au=',num2str(fracs_2(i))));
    if i ~= 5
        set(gca,'xtick',[])
    end
    if i == 5
       xlabel("Wavelength (nm)")
    end
    ylim([-1.2,1.2])
    xlim([300,800])
    if i==3
        ylabel("CD (arb. units)")
    end 

    str=PlotLegends{i};
    text(500,0.6,str)
    set(gca, 'FontName', 'Arial')
    set(gca,'FontSize',8)
    box on
    set(gca,'Linewidth',1)
    
    if i==1
        plot(ExpDataArray0(:,1),ExpDataArray0(:,2)/MAXEXP, expcol{1},'LineWidth',1)
        plot(ExpDataArray0(:,1),-ExpDataArray0(:,2)/MAXEXP, expcol{2},'LineWidth',1)
    elseif i==2
        plot(ExpDataArray05(:,1),ExpDataArray05(:,2)/MAXEXP,expcol{1},'LineWidth',1)
        plot(ExpDataArray05(:,1),-ExpDataArray05(:,2)/MAXEXP,expcol{2},'LineWidth',1)
    elseif i==3
        plot(ExpDataArray10(:,1),ExpDataArray10(:,2)/MAXEXP,expcol{1},'LineWidth',1)
        plot(ExpDataArray10(:,1),-ExpDataArray10(:,2)/MAXEXP,expcol{2},'LineWidth',1)
    elseif i==4
        plot(ExpDataArray15(:,1),ExpDataArray15(:,2)/MAXEXP,expcol{1},'LineWidth',1)
        plot(ExpDataArray15(:,1),-ExpDataArray15(:,2)/MAXEXP,expcol{2},'LineWidth',1)
    else
        plot(ExpDataArray2(:,1),ExpDataArray2(:,2)/MAXEXP,expcol{1},'LineWidth',1)
        plot(ExpDataArray2(:,1),-ExpDataArray2(:,2)/MAXEXP,expcol{2},'LineWidth',1)
    end
end

box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')


subplot(3,2,6); set(gca,'units','inches','position',[3.65,0.35,2.65,0.425])
subplot(3,2,4); set(gca,'units','inches','position',[3.65,1.3,2.65,0.425])
subplot(3,2,5); set(gca,'units','inches','position',[3.65,0.825,2.65,0.425])
subplot(3,2,2); set(gca,'units','inches','position',[3.65,2.25,2.65,0.425])
subplot(3,2,3); set(gca,'units','inches','position',[3.65,1.775,2.65,0.425])
subplot(3,2,1);set(gca,'units','inches','position',[0.375,0.35,2.85,1.8])

annotation('textbox', [0, 1, 0, 0], 'string', '(A)')
annotation('textbox', [0.5, 1, 0, 0], 'string', '(B)')
