close all;
clear all;
update_files=true;
validate_results=false;

%% Analylitical Results
lambda = 400:1:800;
% Au permittivity
for j = 1:length(lambda)
    [eps1(j), eps2(j)] =  getEpsAuByLambda(lambda(j), 10e3);
end
e_Au = eps1 + 1i*eps2;

%Background permittivity
e_w_array =1.75:0.15:2.5;
% Array initialisation
Absall=zeros(length(e_w_array),length(lambda));
CDall=zeros(length(e_w_array),length(lambda));
% Outer dimensions of chiral Nanorod
Length = 30;
Width = 10;
NoOdEllipsoids=1e12;
% Chirality parameter
[CL,e_c] = calcChiralParam(lambda);
%maxwell garnett to get effective medium parameters
frac =0.05;
e_eff = e_Au .*(2*frac*(-e_Au+e_c)+e_c+2*e_Au)./(e_c+2*e_Au-frac*(e_c-e_Au));
CL_modified = frac*CL%3*frac*(CL.*e_Au./(e_c+2*e_Au-frac*(e_c-e_Au)));

linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};

for count=1:1:length(e_w_array)
    e_w = e_w_array(count);
    [AbsL,AbsR]= calcAbsN(e_w, lambda, e_eff, Length/2, Width/2, CL_modified,NoOdEllipsoids);
    Absall(count,:)=AbsL;
    
    CD = AbsR-AbsL;
    CDall(count,:)=CD;
end
ABSscale=max(max(Absall));
CDscale=max(max(CDall));
%% Abs vs medium permittivity plot
figure
for count=1:1:length(e_w_array)
    plot(lambda, Absall(count,:)/ABSscale,linStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(e_w_array(count))));
    hold on
end
leg = legend('Location','northwest','NumColumns',2);
legend show
leg.ItemTokenSize =  [16,10];
xlabel('Wavelength (nm)');
ylabel('Absorption (arb. units)');
ylim([0,1])
set(gcf,'units','inches','position',[0.5,0.5,3.33,2.05])
%set(gca,'units','inches','position',[0.325,0.3,2.85,1.7])%for fontsize 8
set(gca,'units','inches','position',[0.35,0.35,2.85,1.65])%for fontsize 9
box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
if(update_files)
    %saveas(gcf,'Fig_2_abs_medium','epsc')
    %print('abs_medium','-depsc', '-opengl', '-r600');
end

%% CD vs medium permittivity plot
figure
for count=1:1:length(e_w_array)
    plot(lambda, CDall(count,:)/CDscale,linStyles{count},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(e_w_array(count))));
    hold on
end
%ylim([-0.6,1])
leg = legend('Location','northwest','NumColumns',2);
legend show
leg.ItemTokenSize = [16,10];
xlabel('Wavelength (nm)');
ylabel('CD (arb. units)');
%set(gcf,'units','inches','position',[0.5,0.5,3.33,2.05])
%set(gca,'units','inches','position',[0.4,0.3,2.75,1.7])%for fontsize 8
%set(gca,'units','inches','position',[0.4,0.35,2.75,1.65])%for fontsize 9

box on
set(gca,'Linewidth',1)
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
if(update_files)
    %saveas(gcf,'Fig_3_cd_medium','epsc')
    %print('cd_medium','-depsc', '-opengl', '-r600');
end


%% Simulated results comsol
if (validate_results)
    % 30nm 10nm 10nm radii ellipsoid
    % Simulated Results Comsol
    Folderpath = "C:\PhD Research\Ellipsoid chiral\ExactChiralEllipsoid\newComsol\";
    Filenames = ["comsol_2.5_15_5.csv","comsol_2.35_15_5.csv","comsol_2.2_15_5.csv","comsol_2.05_15_5.csv","comsol_1.9_15_5.csv","comsol_1.75_15_5.csv"];
    %Filenames = ["comsol_5_15_2.35.csv","comsol_5_15_2.05.csv","comsol_5_15_1.75.csv"];

    h = zeros(1,2*length(Filenames));
    k = zeros(1,2*length(Filenames));

    e_medium = [2.5,2.35,2.2,2.05,1.9,1.75];
    %e_medium = [2.35,2.05,1.75];
    colors=[         
        0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
    linStyles = {'-o','-+','-*','-<','-s','-d','-^','-v','->','-p','-x','-h'};
    linStyles_dashes = {'--o','--+','--*','--<','--s','--d','--^','--v','-->','--p','--x','--h'};

    %ABSscale_comsol =1.3348e-21;
    %ABSscale_mie = 8.4156e-05;

    lambda_comsol=zeros(81,length(Filenames));
    ABSL_comsol=zeros(81,length(Filenames));
    ABSR_comsol=zeros(81,length(Filenames));
    CD_comsol=zeros(81,length(Filenames));

    for i=1:length(Filenames)
        Filepath = strcat(Folderpath,Filenames(i));
        comsolResultsTable = readtable(Filepath);
        comsolResultsArray = table2array(comsolResultsTable(1:81,1:5));
        lambda_comsol(:,i) =comsolResultsArray(:,1).*10^9;
        ABSL_comsol(:,i)=(2*comsolResultsArray(:,2)+comsolResultsArray(:,4))/3;%/max((comsolResultsArray(:,2)+comsolResultsArray(:,4)));
        ABSR_comsol(:,i)=(2*comsolResultsArray(:,3)+comsolResultsArray(:,5))/3;%/max((comsolResultsArray(:,3)+comsolResultsArray(:,5)));
        CD_comsol(:,i)=-ABSL_comsol(:,i)+ABSR_comsol(:,i);
    end

    ABSscale_comsol=max(max(ABSL_comsol));
    CDscale_comsol=max(max(abs(CD_comsol)));
    figure
    set(gcf,'units','inches','position',[0.5,0.5,3.33,2.9])
    subplot(2,1,2);
    for i=1:1:length(Filenames)
        h(i)=plot(lambda_comsol(:,i), ABSL_comsol(:,i)/ABSscale_comsol,strcat(linStyles{i}),'color',colors(i,:),'LineWidth',1,'MarkerIndices',1:10:length(lambda_comsol(:,i)),'DisplayName',strcat(char(949),'_m=',num2str(e_medium(i))));
        hold on%,
    end

    ylim([-0.05,1.05]);
    set(gca,'units','inches','position',[0.425,0.35,2.8,1])%for fontsize 9 [0.425,1.375,2.8,1]
    
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    %annotation('textbox', [0, 0.9, 0, 0], 'string', '(A)')
    text(410,0.9,'(B) Comsol')%'fontweight', 'bold'
    box on
    set(gca,'Linewidth',1)
    
    aax1=subplot(2,1,1);
    for count=length(e_w_array):-1:1
        h(length(e_medium)+count)=plot(lambda, Absall(count,:)/ABSscale,strcat(linStyles{length(e_w_array)-count+1}),'color',colors(length(e_w_array)-count+1,:),'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(e_w_array(count))));
        hold on
    end
    
    text(410,0.9,'(A) Ex. Mie Gans')%'fontweight', 'bold'
    leg = legend('units','inches','Position',[0.425, 2.4,2.8, 0.5],'NumColumns',3)
    %leg = legend('Location','northwest','NumColumns',2);
    legend show
    leg.ItemTokenSize = [20,30];
    leg.Position=[0.425, 2.4,2.8, 0.5];
    xticklabels(aax1,{});
    ylim([-0.05,1.05]);
    ylabel('Absorption')
    set(gca,'units','inches','position',[0.425,1.375,2.8,1])%for fontsize 9
    
   
    %ylim([0,1]);
    %set(gcf,'units','inches','position',[0.5,0.5,3.33,2.05])
    %set(gca,'units','inches','position',[0.34,0.34,2.85,1.65])
    %set(gcf,'units','inches','position',[0.5,0.5,3.33,2.05])
    %set(gca,'units','inches','position',[0.4,0.3,2.75,1.7])%for fontsize 8
    %set(gca,'units','inches','position',[0.4,0.35,2.75,1.65])%for fontsize 9
    %leg = legend('Location','northwest','NumColumns',2);

    
    box on
    set(gca,'Linewidth',1)
    set(findall(gcf,'-property','FontSize'),'FontSize',9)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    if(update_files)
        saveas(gcf,'Fig_2_abs_with_medium_new_val_3','epsc')
        print('abs_with_medium_validation','-depsc', '-opengl', '-r600');
    end
    %% Cd plots
    figure
    set(gcf,'units','inches','position',[0.5,0.5,3.33,2.9])
    %t = tiledlayout(2,1); % Requires R2019b or later
    %t.TileSpacing = 'none';
    %t.Padding = 'compact';
    %ax1 = nexttile;
    ax1=subplot(2,1,1);
    for i=1:length(e_medium)
        h(i)=plot(lambda, CDall(length(e_medium)-i+1,:)/CDscale,linStyles{i},'LineWidth',1,'MarkerIndices',1:60:length(lambda),'DisplayName',strcat(char(949),'_m=',num2str(e_medium(i))));
        hold on
    end
    text(410,0.9,'(A) Ex. Mie-Gans')%'fontweight', 'bold'
%     legend(h(4:1:6))
%     leg = legend('Location','northwest');
%     legend show
%     leg.ItemTokenSize = [20,30];
%     legend(k(1:length(e_medium)))
%     leg = legend('Location','northwest');
%     legend show
    ylim([-0.7,1.1])
    
    %ylabel('CD (Ex. Mie Gans) (arb. units)');
    fig2ylabel=ylabel('CD (arb. units)');
    
    %set(gcf,'units','inches','position',[0.5,0.5,5.5,2.5])
    set(gca,'units','inches','position',[0.425,1.375,2.8,1])
    xticklabels(ax1,{})
    box on
    
    
    ax2=subplot(2,1,2);
    for i=1:length(Filenames)
        k(i)=plot(lambda_comsol(:,i), CD_comsol(:,i)/CDscale_comsol,linStyles{i},'MarkerIndices',1:10:length(lambda_comsol(:,i)),'LineWidth',1,'DisplayName',strcat(char(949),'_m=',num2str(e_medium(i))));
        hold on%,
    end
    text(410,0.9,'(B) Comsol')%'fontweight', 'bold'
    legend(k(1:1:6))
    %leg = legend('Location','northoutside','NumColumns',3)
    leg = legend('units','inches','Position',[0.425, 2.4, 2.8, 0.5],'NumColumns',3)
    xlabel('Wavelength (nm)');
%legend('Position',[0 0 0.1 0.2],'NumColumns',3)
    legend show
    leg.ItemTokenSize = [20,30];
    leg.Position=[0.425, 2.4, 2.8, 0.5];
    %ylabel('CD (Comsol)(arb. units)');
    fig1ylabel=ylabel('CD (arb. units)');
 %   xlabel('Wavelength (nm)');
    
    %set(gca,'units','inches','position',[0.5,0.5,2.5,3])
    
    set(gca,'units','inches','position',[0.425,0.35,2.8,1])
    box on
    set(gca,'Linewidth',1)
    ylim([-0.45,1.1])
    %ylim([-0.7,1.1])

    fig1ylabel.Position(1) = fig2ylabel.Position(1);
    set(findall(gcf,'-property','Linewidth'),'Linewidth',1)
    set(findall(gcf,'-property','FontSize'),'FontSize',9)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')

    if(update_files)
        saveas(gcf,'Fig_3_cd_with_medium_new_val_2','epsc')
        print('cd_with_medium_validation','-depsc', '-opengl', '-r600');
    end
end
