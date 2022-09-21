clear
clc

R = georasterref;
R.RasterSize = [25,22];
R.Latlim = [25, 49];
R.Lonlim = [-109, -88];
R.ColumnsStartFrom = 'south';
R.RowsStartFrom = 'west';

thor = importdata('ThOR_SMERGE_Range_2005_2017.csv'); 
thor = thor.data; 
thor(thor(:,5) < 5 | thor(:,5) > 9, :) = []; 
thor(thor(:,8) < 12 | thor(:,8) > 20, :) = []; 

anom = importdata('Control_Anomalies.csv'); 
% anom(isnan(anom(:,9)), :) = []; 

%Use 2x2 degree boxes and assess significance of the event-control
%difference
tdiff = NaN(25,22); perc = NaN(25,22); bfCoef = NaN(25,22);
rowCount = 1; 
sig = []; 
for i = 25:49
    colCount = 1; 
    for ii = -109:-89
        ii
        tic
        r1 = find(anom(:,2) >= i-1 & anom(:,2) < i+1); 
        r2 = find(anom(:,3) >= ii-1 & anom(:,3) < ii+1); 
        r1 = r1(ismember(r1,r2));         

        if isempty(r1) == 0 && sum(~isnan(anom(r1,9))) > 50
            tdiff(rowCount,colCount) = mean(thor(r1,13),'omitnan')-mean(reshape(anom(r1,9:end),[403*size(r1,1), 1]),'omitnan');         
            samples = reshape(anom(r1,9:end),[403*size(r1,1),1]); 

            if i < 40
                if tdiff(rowCount,colCount) <= -0.003 || tdiff(rowCount,colCount) >= 0.003
                    [lat1,lon1] = pix2latlon(R,rowCount,colCount); 
                    sig = [sig; lat1 lon1]; 
                    clear lat1 lon1
                end

            end

            %--------------------------------------------------------------
            %Bootstrapping Bayes Factor Code Snippet
%             bsDiff = NaN(1000,1); 
%             for j = 1:1000                
%                 x = round(rand(size(r1,1),1).*size(samples,1));
%                 x(x==0) = 1; x(x>size(samples,1)) = size(samples,1); 
%                 Ye = samples(x,1); 
%                 Yc = samples; Yc(x, :) = []; 
%                 
%                 bsDiff(j,1) = mean(Ye,'omitnan')-mean(Yc,'omitnan'); 
%                 clear x Ye Yc                
%             end
%             [c,ind] = min(abs(tdiff(rowCount,colCount)-prctile(bsDiff,(0:100)')));
%             perc(rowCount,colCount) = ind; 
%             [bfCoef(rowCount,colCount),pval] = bf.ttest2(tdiff(rowCount,colCount),bsDiff); 
            %--------------------------------------------------------------
        
            clear samples c ind pval

        end
        clear r1 r2
        colCount = colCount + 1; 
        clc
        toc
    end
    rowCount = rowCount + 1;    

end
%-----------------------------------------------
tdiff(isnan(tdiff)) = 0; 
% perc(isnan(perc)) = 50; 
% bfCoef(isnan(bfCoef)) = 1; 
% 

red_blue = load('red_blue.mat'); 
red_blue = red_blue.red_blue; 

states = shaperead('usastatelo','UseGeoCoords',true); 
land = shaperead('landareas.shp','UseGeoCoords',true); 
canada = shaperead('C:\Users\twford\Box\Research\Datasets\GIS Data\Canada\PROVINCE.SHP','UseGeoCoords',true); 
% ocean = shaperead('C:\Users\twford\Box\Research\Datasets\GIS Data\Oceans\ne_110m_ocean.shp','UseGeoCoords',true);
% subplot(1,2,1); 
ax = usamap([25 50],[-109 -89]);
setm(ax,'MapProjection','lambert');
setm(ax,'MeridianLabel','off','ParallelLabel','off');
framem('off');
geoshow(tdiff,R,'DisplayType','TextureMap'); 
geoshow(states,'DisplayType','Polygon','FaceColor','none','EdgeColor','k');
geoshow(land,'DisplayType','Polygon','FaceColor','none','EdgeColor','k');
geoshow(canada,'DisplayType','Polygon','FaceColor','wh','EdgeColor','none');
% setm(gca,'ffacecolor',rgb('ocean blue'))

scatterm(sig(:,1),sig(:,2),50,'k','*'); 
colormap(ax,red_blue./255); 
colorbar(); 
caxis([-0.012 0.012]); 












