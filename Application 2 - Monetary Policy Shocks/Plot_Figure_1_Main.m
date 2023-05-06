clear; clc;
n = 6; k =1 ; quants = [.05,.16,.5,.84,.95];
nthin = [5,50,5];
hor = 5*12; % 5 years of IRFs
addpath('data','functions','results')
load('ACRr1')
outputR1 = output;
irfR1 = zeros(hor+1,n,size(outputR1.alphas,2)/nthin(1))  ; 
a = 1;
for i = 1:nthin(1):size(outputR1.alphas,2)  
    Btil = reshape(outputR1.vecBs(:,i), n+1, n+1); 
    irf =  IRF(outputR1.alphas(:,i),vec(Btil(:,[n,1:n-1,n+1])),output.input.p,n+k,hor,n,n, output.input.c);
    irfR1(:,:,a) = irf(:,1:n); 
    a = a + 1;
end 
IRFs(:,:,:,1) = squeeze(quantile(irfR1,quants,3)); 
load('ACRr2') 
a = 1;
outputR2 = output;
irfR2    = zeros(hor+1,n,size(outputR2.alphas,2)/nthin(2))  ; 
for i = 1:nthin(2):size(outputR2.alphas,2)  
    irf =  IRF(outputR2.alphas(:,i),outputR2.vecBs(:,i),output.input.p,n+k,hor,n,n, output.input.c);
    irfR2(:,:,a) = irf(:,1:n);
    a = a + 1;
end 
IRFs(:,:,:,2) = squeeze(quantile(irfR2,quants,3));  
load('ACRr3')  
irfR3    = zeros(hor+1,n,size(output.alphas,2)/nthin(3))  ;  
a = 1; 
for i = 1:nthin(3):size(output.alphas,2)  
    irf =  IRF(output.alphas(:,i),output.vecBs(:,i),output.input.p,n+k,hor,n,n, output.input.c);
    irfR3(:,:,a) = irf(:,1:n); 
    a = a + 1;
end
IRFs(:,:,:,3) = squeeze(quantile(irfR3,quants,3)); 

%% Go Figure
hFig2 = figure(2);
horizon = (0:hor)';
xgraph = [horizon',fliplr(horizon')];
varnames = {' output', ' prices' , ' commodity prices', ' total reserves', ' nonborrowed reserves', ' federal funds rate'}; 
d1 = 80; d2 = 190;
color_1 = [80, 80, 80]./255;
color_2 = [190, 190, 190]./255;
color_3 = [240,240,240]./255;
alpha = 0.3; alpha2 = .6;
color_plot = [80,80,80]./255;
a = 1;
model ={'R1','R2','R3'};
ylims_fig1 = [ [-1.5,0.5] ; [-0.4,0.6]; [-4,4]; [-2,2]; [-2,2];[-.5,1]]';  

for i = 1:3
    IRFi = IRFs(:,:,:,i); 
    for j = 1:n
        subplot(3,n,a)
        plot(horizon, squeeze(IRFi(:,j,3)), 'k','LineWidth', 1.3,'markersize',3)
        hold on; grid on;
        plot(horizon, squeeze(IRFi(:,j,[1,2,4,5])), 'k--','LineWidth', .8,'markersize',3)  
        Y_1 = [IRFi( :, j,1)', fliplr(IRFi(:,j,5)')];
        bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
        set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha,'EdgeAlpha',alpha);
        Y_1 = [IRFi( :, j,2)', fliplr(IRFi(:,j,4)')];
        bounds = fill(xgraph,Y_1,color_1,'LineStyle','none');
        set(bounds,'FaceColor',color_1,'EdgeColor',color_1,'FaceAlpha',alpha2,'EdgeAlpha',alpha2); 
        ylim([ylims_fig1(1,j),ylims_fig1(2,j)])
        a = a + 1;
        xticks([0,12:12:hor])
        xlim([horizon(1),horizon(end)])
        line([horizon(1),horizon(end)],[0,0],'color',[105,105,105]./255 )
        title(varnames{j},'Interpreter','latex')
        hold off
        if j==1
            ylabel(model{i}, 'Interpreter', 'latex')
        end
    end 
end


set(gcf,'PaperPositionMode','auto') 
set(hFig2, 'Position', [30 50 1400 700])
hFig1a = tightfig(hFig2); 
print(hFig2,strcat('figures/IRF_baseline'), '-painters' ,'-dpdf') 



 