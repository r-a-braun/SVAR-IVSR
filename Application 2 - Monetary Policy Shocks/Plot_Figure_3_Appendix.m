clear; clc;
addpath('data','functions','results')

finvars_name = {'real stock prices',' mortgage spread',' commercial paper spread',' excess bond premium'};
finvars_name_compact = {'SP','MS','CPS','EBP'};
n = 7; k =1 ; quants = [.05,.16,.5,.84,.95];
hor = 5*12; % 5 years of IRFs
finvar = 4;
nthin= 5;
IRFs = NaN(hor+1,size(quants,2),size(finvars_name_compact,2),2);
for i = 1:size(finvars_name_compact,2) 
    disp(finvars_name_compact{i})
    for restriction = 2:3
        load(strcat('ACRr',num2str(restriction),finvars_name_compact{i}),'output')
        irfir = zeros(hor+1,size(output.alphas,2)/nthin)  ;
        a = 1;
        for ii = 1:nthin:size(output.alphas,2) 
            irf =  IRF(output.alphas(:,ii),output.vecBs(:,ii),output.input.p,n+k,hor,n,n, output.input.c);
            irfir(:,a) = irf(:,finvar);
            a = a + 1;
        end
        IRFs(:,:,i,restriction-1) = squeeze(quantile(irfir,quants,2)); 
    end
end


disp('done')


hFig3 = figure(3);
horizon = (0:hor)';
xgraph = [horizon',fliplr(horizon')];
d1 = 80; d2 = 190;
color_1 = [80, 80, 80]./255;
color_2 = [190, 190, 190]./255;
color_3 = [240,240,240]./255;
alpha = 0.3; alpha2 = .6;
color_plot = [80,80,80]./255;
a = 1;
ylims_fig1 = [ [-4,3] ; [-0.2,0.2]; [-.2,.2]; [-.2,.2]]';  

model ={'R2','R3'};
for i = 1:size(finvars_name,2)
    IRFi = squeeze(IRFs(:,:,i,:));%IRFs(:,:,:,i);
    
    subplot(2,size(finvars_name,2),a)
    plot(horizon, IRFi(:,3,1), 'k','LineWidth', 1.3,'markersize',3)
    hold on; grid on;
    plot(horizon, IRFi(:,[1,2,4,5],1) , 'k--','LineWidth', .8,'markersize',3)
    Y_1 = [IRFi( :, 1,1)', fliplr(IRFi(:,5,1)')];
    bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
    set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha,'EdgeAlpha',alpha);
    Y_1 = [IRFi( :, 2,1)', fliplr(IRFi(:,4,1)')];
    bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
    set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha2,'EdgeAlpha',alpha2);
    xticks([0,12:12:hor])
    ylim([ylims_fig1(1,i),ylims_fig1(2,i)])
    xlim([horizon(1),horizon(end)])
    line([horizon(1),horizon(end)],[0,0],'color',[105,105,105]./255 )
    title(finvars_name{i},'Interpreter','latex')
    hold off
    if  i==1
        ylabel(model{1}, 'Interpreter', 'latex')
    end
    
    subplot(2,size(finvars_name,2),size(finvars_name,2)+a)
    plot(horizon, IRFi(:,3,2), 'k','LineWidth', 1.3,'markersize',3)
    hold on; grid on;
    plot(horizon, IRFi(:,[1,2,4,5],2) , 'k--','LineWidth', .8,'markersize',3)
    Y_1 = [IRFi( :, 1,2)', fliplr(IRFi(:,5,2)')];
    bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
    set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha,'EdgeAlpha',alpha);
    Y_1 = [IRFi( :, 2,2)', fliplr(IRFi(:,4,2)')];
    bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
    set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha2,'EdgeAlpha',alpha2);
    xticks([0,12:12:hor])
    ylim([ylims_fig1(1,i),ylims_fig1(2,i)])
    xlim([horizon(1),horizon(end)])
    line([horizon(1),horizon(end)],[0,0],'color',[105,105,105]./255 )
    title(finvars_name{i},'Interpreter','latex')
    
    
    a = a + 1;
    
    hold off
    if  i==1
        ylabel(model{2}, 'Interpreter', 'latex')
    end
end
set(gcf,'PaperPositionMode','auto') 
set(hFig3, 'Position', [30 50 1400 700])
hFig3 = tightfig(hFig3); 
print(hFig3,strcat('figures/IRF_finvars'), '-painters' ,'-dpdf') 
