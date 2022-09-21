clear
clc

%Opening LCL heights from ThOR events
lcl = importdata('ThOR_SMERGE_LCL.csv'); lcl = lcl.data; 
%Opening PBL heights from ThOR events
pbl = importdata('ThOR_SMERGE_PBL.csv'); pbl = pbl.data;

%Open "flux" which contains all soil moisture and surface heat flux data
flux = importdata('ThOR_SMERGE_SurfaceFlux.csv'); flux = flux.data; 
lcl(lcl(:,8) > 18, :) = []; 
pbl(pbl(:,8) > 18, :) = []; 

%-------------------------------------
%Define the geographic regions
texas = [28 33; -102 -96]; 
midSouth = [33 38; -94 -89]; 
central = [34 40; -103 -95]; 

region = texas; %do each regional analysis separately
%Subset each dataset to that specific region
ctpHI = lcl(lcl(:,2) >= region(1,1) & lcl(:,2) <= region(1,2) & lcl(:,3) >= region(2,1) & lcl(:,3) <= region(2,2), :); 
pblSub = pbl(pbl(:,2) >= region(1,1) & pbl(:,2) <= region(1,2) & pbl(:,3) >= region(2,1) & pbl(:,3) <= region(2,2), :); 
evapFrac = flux(flux(:,2) >= region(1,1) & flux(:,2) <= region(1,2) & flux(:,3) >= region(2,1) & flux(:,3) <= region(2,2), :); 

%Retain only LCL and PBL heights 6 hours prior to initiation
lastVal = find_ndim(~isnan(ctpHI),2,'last'); %Find columns with hour right before initiation
for j = 1:size(ctpHI,1)
    lC(j,:) = ctpHI(j,lastVal(j)-6:lastVal(j)); 
    pC(j,:) = pblSub(j,lastVal(j)-6:lastVal(j)); 
    pC(pC<0) = 0; 
    LPdiff(j,:) = lC(j,:)-pC(j,:); 
end

%-----------------------------------------
%Only keep data where CTP is within reasonable range (-500 to 500 J/kg)
r1 = [find(ctpHI(:,10) < -500); find(ctpHI(:,10) > 500)]; 
evapFrac(r1,:) = []; 
ctpHI(r1,:) = []; 
pC(r1,:) = []; 
lC(r1,:) = []; 
clear r1


lC(lC<0) = 0; 
%Only include events with valid soil moisture and surface heat flux
evapFrac(isnan(ctpHI(:,9)), :) = []; 
pC(isnan(ctpHI(:,9)), :) = []; 
lC(isnan(ctpHI(:,9)), :) = []; 
ctpHI(isnan(ctpHI(:,9)), :) = []; 

w1 = find(ctpHI(:,9) > 0); %wet soil events
% r1 = find(ctpHI(:,9) > -0.01 & ctpHI(:,9) < 0.01); 
d1 = find(ctpHI(:,9) <= 0); %dry soil events

blues = importdata('Brew_blues.mat'); 
reds = importdata('Brew_reds.mat'); 

%------------------------------------------------------------------------
%CTP-HI Area and Histogram Plots
% figure
% group(w1,1) = 1; group(d1,1) = 2; 
% h = scatterhist(ctpHI(:,10),ctpHI(:,11),'Group',group,'kernel','on'); 
% clr = get(h(1),'colororder');
% 
% xgrid = linspace(-500,500,100); 
% ygrid = linspace(0,35,100); 
% [x1,y1] = meshgrid(xgrid,ygrid); 
% xi = [x1(:) y1(:)]; 
% 
% blues = importdata('Brew_blues.mat'); 
% reds = importdata('Brew_reds.mat'); 
% 
% [f,ep] = ksdensity([ctpHI(w1,10) ctpHI(w1,11)],xi); 
% X = reshape(ep(:,1),length(xgrid),length(ygrid)); 
% Y = reshape(ep(:,2),length(xgrid),length(ygrid)); 
% Z = reshape(f,length(xgrid),length(ygrid)); 
% 
% ax1 = axes; 
% [C1,h1] = contourf(X,Y,Z,2); 
% clear f ep X Y Z
% 
% [f,ep] = ksdensity([ctpHI(d1,10) ctpHI(d1,11)],xi); 
% X = reshape(ep(:,1),length(xgrid),length(ygrid)); 
% Y = reshape(ep(:,2),length(xgrid),length(ygrid)); 
% Z = reshape(f,length(xgrid),length(ygrid)); 
% 
% ax2 = axes; 
% [C2,h2] = contourf(X,Y,Z,2); 
% linkaxes([ax1,ax2]); 
% ax2.XTick = []; 
% ax2.YTick = []; 
% colormap(ax1,blues./255); 
% colormap(ax2,reds./255); 
% 
% hFills = h2.FacePrims;
% hFills(3).ColorType = 'truecoloralpha'; 
% hFills(3).ColorData(4) = 135; 
% hFills(2).ColorType = 'truecoloralpha'; 
% hFills(2).ColorData(4) = 135; 
% hFills(1).ColorType = 'truecoloralpha'; 
% hFills(1).ColorData(4) = 0; 
%--------------------------------------------------------------------
%Evaporative Fraction & PBL/LCL plot

% figure
% group(w1,1) = 1; group(d1,1) = 2; 
% h = scatterhist(evapFrac(:,19),mean(pC,2),'Group',group,'kernel','on');
% clr = get(h(1),'colororder');
% 
% xgrid = linspace(0,1,100); 
% ygrid = linspace(-500,2500,100); 
% [x1,y1] = meshgrid(xgrid,ygrid); 
% xi = [x1(:) y1(:)]; 
% 
% blues = importdata('Brew_blues.mat'); 
% reds = importdata('Brew_reds.mat'); 
% 
% [f,ep] = ksdensity([evapFrac(w1,19) mean(pC(w1,:),2)],xi); 
% X = reshape(ep(:,1),length(xgrid),length(ygrid)); 
% Y = reshape(ep(:,2),length(xgrid),length(ygrid)); 
% Z = reshape(f,length(xgrid),length(ygrid)); 
% 
% ax1 = axes; 
% [C1,h1] = contourf(X,Y,Z,2); 
% clear f ep X Y Z
% 
% [f,ep] = ksdensity([evapFrac(d1,19) mean(pC(d1,:),2)],xi); 
% X = reshape(ep(:,1),length(xgrid),length(ygrid)); 
% Y = reshape(ep(:,2),length(xgrid),length(ygrid)); 
% Z = reshape(f,length(xgrid),length(ygrid)); 
% 
% ax2 = axes; 
% [C2,h2] = contourf(X,Y,Z,2); 
% linkaxes([ax1,ax2]); 
% ax2.XTick = []; 
% ax2.YTick = []; 
% colormap(ax1,blues./255); 
% colormap(ax2,reds./255); 
% 
% hFills = h2.FacePrims;
% hFills(3).ColorType = 'truecoloralpha'; 
% hFills(3).ColorData(4) = 135; 
% hFills(2).ColorType = 'truecoloralpha'; 
% hFills(2).ColorData(4) = 135; 
% hFills(1).ColorType = 'truecoloralpha'; 
% hFills(1).ColorData(4) = 0; 
 
%----------------------------------------------------------------------
%PBL & LCL evolution figure

% figure
% count = 1; 
% for i = 1:2    
%     for p = 1:3
%         subplot(2,3,count)
%         if i == 1
%             rows = w1; 
%         else
%             rows = d1; 
%         end
%         pblSub = pC(rows,:); 
%         efSub = evapFrac(rows,:); 
%         lclSub = lC(rows,:); 
% 
%         if p == 1
%             r1 = find(efSub(:,19) < 0.33); 
%         elseif p == 2
%             r1 = find(efSub(:,19) >= 0.33 & efSub(:,19) <= 0.66); 
%         else
%             r1 = find(efSub(:,19) > 0.66);
%         end
%         
%         efSub = efSub(r1,:); 
%         pblSub = pblSub(r1,:); 
%         lclSub = lclSub(r1,:); 
% 
%         x = (1:7); 
%         ym = median(pblSub,'omitnan'); 
%         y1 = prctile(pblSub,25,1); 
%         y2 = prctile(pblSub,75,1); 
%         y3 = prctile(pblSub,10,1); 
%         y4 = prctile(pblSub,90,1); 
% 
%         plot(x,y3); hold on
%         xlim([1 7]); 
%         plot(x,y4); hold on
%         h(1) = patch([x fliplr(x)],[y3 fliplr(y4)],[0.85 0.50 0.60]); 
%         hold on
% 
%         plot(x,y1); hold on
%         plot(x,y2); hold on
%         h(2) = patch([x fliplr(x)],[y1 fliplr(y2)],[0.65 0.35 0.40]); 
%         hold on
% 
%         plot(ym,'Color','k','LineWidth',2); %daily median
% 
%         clear ym y1 y2 y3 y4 h
%         ym = median(lclSub,'omitnan'); 
%         y1 = prctile(lclSub,25,1); 
%         y2 = prctile(lclSub,75,1); 
%         y3 = prctile(lclSub,10,1); 
%         y4 = prctile(lclSub,90,1); 
% 
%         plot(x,y3); hold on
%         xlim([1 7]); 
%         plot(x,y4); hold on
%         h(1) = patch([x fliplr(x)],[y3 fliplr(y4)],[0.5 0.6 0.85]); 
%         alpha(h(1),0.5)
%         hold on
% 
%         plot(x,y1); hold on
%         plot(x,y2); hold on
%         h(2) = patch([x fliplr(x)],[y1 fliplr(y2)],[0.35 0.40 0.65]); 
%         alpha(h(2),0.5)
%         hold on
% 
%         plot(ym,'Color','k','LineWidth',2); %daily median
% 
%         clear ym y1 y2 y3 y4 h
% 
%         title([num2str(size(pblSub,1)) ' Events'],'FontName','Arial','FontSize',20); 
%         set(gca,'FontSize',20); 
%         xticks(1:6); 
%         xticklabels({'1' '2' '3' '4' '5' '6'});
%         if p == 1
%             ylabel('Height (m)'); 
%         end
%         if count == 5
%             xlabel('Hours Before Initiation'); 
%         end
%         if count == 1
%             xlabel('EF < 0.33'); 
%         end
%         if count == 2
%             xlabel('0.33 <= EF <= 0.66'); 
%         end
%         if count == 3
%             xlabel('EF > 0.66'); 
%         end
%         ylim([-500 2500])
%         count = count + 1;
%     end
% end
%----------------------------------------------------------------------
%Alternative PBL & LCL evolution figure - with boxplots


figure
count = 1; 
for i = 1:2    
    allSub = []; 
    for p = 1:3
        subplot(2,3,count)
        if i == 1
            rows = w1; 
        else
            rows = d1; 
        end
        pblSub = pC(rows,:); 
        efSub = evapFrac(rows,:); 
        lclSub = lC(rows,:); 

        if p == 1
            r1 = find(efSub(:,19) < 0.33); 
        elseif p == 2
            r1 = find(efSub(:,19) >= 0.33 & efSub(:,19) <= 0.66); 
        else
            r1 = find(efSub(:,19) > 0.66);
        end
        
        efSub = efSub(r1,:); 
        pblSub = pblSub(r1,:); 
        lclSub = lclSub(r1,:); 

        %-------------------------------------------
        %percent of events with PBL reaching LCL
        sum(pblSub(:,7)>lclSub(:,7))/size(lclSub,1)

        allSub = {pblSub, lclSub}; 
        g = boxplotGroup(allSub); 
        h = g.boxplotGroup(1).Children; 
        b = findobj(h,'Tag','Box'); 
        m = findobj(h,'Tag','Median'); 
        for j = 1:length(b)
            patch(get(b(j),'XData'),get(b(j),'YData'),[0.25 0.35 0.80]);            
            line(get(m(j),'XData'),get(m(j),'YData'),'Color','wh','LineWidth',2); 
        end
        clear h b m

        h = g.boxplotGroup(2).Children; 
        b = findobj(h,'Tag','Box'); 
        m = findobj(h,'Tag','Median'); 
        for j = 1:length(b)
            patch(get(b(j),'XData'),get(b(j),'YData'),[0.75 0.25 0.35]);            
            line(get(m(j),'XData'),get(m(j),'YData'),'Color','wh','LineWidth',2); 
        end
        clear h b m


        title([num2str(size(pblSub,1)) ' Events'],'FontName','Arial','FontSize',20); 
        set(gca,'FontSize',20); 
        xticks(1.5:3:21); 
        xticklabels({'7' '6' '5' '4' '3' '2' '1'});
        if p == 1
            ylabel('Height (m)'); 
        end
        if count == 5
            xlabel('Hours Before Initiation'); 
        end
        if count == 1
            xlabel('EF < 0.33'); 
        end
        if count == 2
            xlabel('0.33 <= EF <= 0.66'); 
        end
        if count == 3
            xlabel('EF > 0.66'); 
        end
        ylim([-500 3000])
        count = count + 1;
    end
end









   