
function sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)

%% By groups
fieldNames = fieldnames(metTypes);
for f = 1:length(fieldNames)
    disp(['Generating figures for: ', char(fieldNames{f})])
    group = metTypes.(fieldNames{f}); groupname = fieldNames{f};
    groupSens = [];
    groupCod = [];
    groupCost = [];
    maxGroupSens = [];
    minGroupSens = [];
    N = [];
    N2 = [];
    for e = 1:length(ignoreMets)
        temp = strmatch(ignoreMets{e},mets);
        if ~isempty(temp)
            group(find(group == temp)) = [];
        end
    end
    if ~isempty(group)
        figure('visible','off');
        subplot(1,2,1)
        hold on
        title('Sensitivity')
        ylabel({'\Delta\mu/\Deltap'})
        
        colormat = round(rand(length(group),3),1);

        growth_change = [];
        parCoeff_change = [];
        maxGroupSens = [];
        minGroupSens = [];
        groupCod = [];
        groupSens = [];
        minPos = [];
        maxPos = [];
        for i = 1:length(group)
            k = group(i);            
            for j = 1:length(growth{k}(1,:))
                x = parMetCoeff{k}(:,j);
                y = growth{k}(:,j);
                
                poly = polyfit(x,y,1);
                dudp{k}(j) = poly(1);
            end
            
            [maxGroupSens(i),maxPos(i)] = max(dudp{k});
            plot(maxPos(i),maxGroupSens(i),'^','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))

            [minGroupSens(i),minPos(i)] = min(dudp{k});
            plot(minPos(i),minGroupSens(i),'s','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))
            
            groupCod(i) = abs(minGroupSens(i)-maxGroupSens(i));
            groupSens(i) = trimmean(dudp{k},10);
%             groupSens(i) = mean(dudp{k});
        end
        
        text(maxPos+0.6,maxGroupSens,metNames(group),'FontSize',8);
        text(minPos+0.6,minGroupSens,metNames(group),'FontSize',8);
        hold off

        subplot(1,2,2)
        barh(groupCod);
        set(gca,'ytick',[1:length(groupCod)]);
        set(gca,'yticklabel',metNames(group));
        xlabel({'\Delta(\Delta\mu/\Deltap)'});
        title('Parameter codependence');
        x0=10;
        y0=10;
        width=1000;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height]);
        saveas(gcf,['Sensitivity_',groupname,'_',ID],'svg');
        
        %% Correlation for Cost
        
        for j = 1:length(costRxn)
            groupCost = biosyntheticCost(model,mets(group),costRxn{j});
            
            figure('visible','off');
            subplot(1,2,1);
            plot(groupCost,groupSens,'*');
            title('Sensitivity');
            text(groupCost+0.1,groupSens,metNames(group),'FontSize');
            xlabel('Biosynthetic Cost');
            ylabel('\Delta\mu/\Deltap');
            
            subplot(1,2,2);
            plot(groupCost,groupCod,'*');
            title('Parameter codependence');
            text(groupCost+0.1,groupCod,metNames(group),'FontSize');
            xlabel('Biosynthetic Cost');
            ylabel('\Delta(\Delta\mu/\Deltap)');
            x0=10;
            y0=10;
            width=1000;
            height=400;
            set(gcf,'units','points','position',[x0,y0,width,height]);
            saveas(gcf,['[',costMet{j},']','Cost_',groupname,'_',ID],'svg');
        end

        %% Correlation with conectivity
        groupConectivity = calculateConectivity(model,mets(group));
        
        figure('visible','off');
        subplot(1,2,1)
        plot(groupConectivity,groupSens,'*')
        title('Sensitivity')
        text(groupConectivity+0.1,groupSens,metNames(group),'FontSize')%,8*2)
        xlabel('Conectivity')
        ylabel('\Delta\mu/\Deltap')
        
        subplot(1,2,2)
        plot(groupConectivity,groupCod,'*')
        title('Parameter codependence')
        text(groupConectivity+0.1,groupCod,metNames(group),'FontSize')%,8*2)
        xlabel('Conectivity')
        ylabel('\Delta(\Delta\mu/\Deltap)')
        x0=10;
        y0=10;
        width=1000;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height])

        saveas(gcf,['Conectivity_',groupname,'_',ID],'svg')
        
        if ~isempty(strmatch('allMets',groupname))
            figure('visible','off')
            plot3(groupConectivity,groupCost,groupSens,'o')
            xlabel('Conectivity')
            ylabel('Cost')
            zlabel('Sensitivity')
            
            figure('visible','off')
            plot3(groupConectivity,groupCost,groupCod,'o')
            xlabel('Conectivity')
            ylabel('Cost')
            zlabel('Codependence')
            save(['3D_',ID],'groupConectivity','groupCost','groupSens','groupCod','metTypes','group','mets','dudp');
        end
    end
end
end