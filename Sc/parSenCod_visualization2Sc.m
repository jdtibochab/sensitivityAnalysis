clear

load('paralleltestSc.mat')

model=model_an;
[metTypes,Met2] = findMetType(model,Met);

%% By groups


eliminateMets = {'h_c','pi_c','atp_c','h2o_c','starch'};

fieldNames = fieldnames(metTypes);
% fieldNames = {'TAG'}
for f = 1:length(fieldNames)
    
    group = metTypes.(fieldNames{f}); groupname = fieldNames{f};
    groupSens = [];
    groupCod = [];
    groupCost = [];
    growth_sens_max = [];
    growth_sens_min = [];
    N = [];
    N2 = [];
    
    for e = 1:length(eliminateMets)
        temp = strmatch(eliminateMets{e},Met);
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
        
        for i = 1:length(group)
            k = group(i);
            %     co2_uptake{k} = abs(co2_uptake{k});
            parMetCoeff{k} = abs(parMetCoeff{k});
            for j=1:it1 % For all BOFs
                growth_mean(i,j) = (max(growth{k}(:,j))+min(growth{k}(:,j)))/2; % Mean value of all points of 1 BOF
                growth_change(i,j)= max(growth{k}(:,j))-min(growth{k}(:,j)); % Change within 1 BOF
                parCoeff_change(i,j) = max(parMetCoeff{k}(:,j))-min(parMetCoeff{k}(:,j)); % Change within 1 BOF
                %         co2_change(i,j) = max(co2_uptake{k}(:,j))-min(co2_uptake{k}(:,j)); % Change within 1 BOF
                %         [poly] = polyfit(parMetCoeff{k}(:,j),co2_uptake{k}(:,j),1);
            end
            
            min_change = min(growth_change(i,:));
            max_change = max(growth_change(i,:));
            
            N_temp = find(growth_change(i,:)==max_change);
            N(i) = N_temp(1);
            
            N2_temp = find(growth_change(i,:)==min_change);
            N2(i) = N2_temp(1);
            
            growth_sens_vector(i,:) = growth_change(i,:)./parCoeff_change(i,:);
            growth_sens_vector(find(isnan(growth_sens_vector))) = 0;
            
            growth_sens_max(i) = growth_change(i,N(i))/parCoeff_change(i,N(i));
            plot(N(i),growth_sens_max(i),'^','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))
            
            
            growth_sens_min(i) = growth_change(i,N2(i))/parCoeff_change(i,N2(i));
            plot(N2(i),growth_sens_min(i),'s','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))
        end
        text(N+0.6,growth_sens_max,Met2(group),'FontSize',8)
        text(N2+0.6,growth_sens_min,Met2(group),'FontSize',8)
        hold off
        
        
        subplot(1,2,2)
        groupCod = abs(growth_sens_min-growth_sens_max);
        barh(groupCod)
        set(gca,'ytick',[1:length(groupCod)]);
        set(gca,'yticklabel',Met2(group));
        xlabel({'\Delta(\Delta\mu/\Deltap)'})
        title('Parameter codependence')
        x0=10;
        y0=10;
        width=1000;
        height=800;
        set(gcf,'units','points','position',[x0,y0,width,height])
        saveas(gcf,['parSenCod_Ch_',groupname],'svg')
        
        %% Correlation for A_nutrient
        groupCost = ATPcost(model,Met(group),'ATPM','adp_c');
        groupSens = mean([growth_sens_max;growth_sens_min]);
        % for i=1:length(group)
        %     metID = findMetIDs(model,Met(group(i)));
        %     groupCost(i) = numAtomsOfElementInFormula(char(model.metFormulas(metID)),'C');
        % end
        
        
        figure('visible','off');
        subplot(1,2,1)
        plot(groupCost,mean([growth_sens_max;growth_sens_min]),'*')
        title('Sensitivity')
        text(groupCost+0.1,mean([growth_sens_max;growth_sens_min]),Met2(group),'FontSize')%,8*2)
        xlabel('ATP cost')
        ylabel('\Delta\mu/\Deltap')
        set(gca,'FontSize')%,8*2)
        
        subplot(1,2,2)
        plot(groupCost,groupCod,'*')
        title('Parameter codependence')
        text(groupCost+0.1,groupCod,Met2(group),'FontSize')%,8*2)
        xlabel('ATP cost')
        ylabel('\Delta(\Delta\mu/\Deltap)')
        x0=10;
        y0=10;
        width=1000;
        height=800;
        set(gca,'FontSize')%,8*2)
        set(gcf,'units','points','position',[x0,y0,width,height])
        saveas(gcf,['parSenCodCorr_Ch_',groupname],'svg')
        
        %% Correlation with conectivity
        groupConectivity = calculateConectivity(model,Met(group));
        
        figure('visible','off');
        subplot(1,2,1)
        plot(groupConectivity,mean([growth_sens_max;growth_sens_min]),'*')
        title('Sensitivity')
        text(groupConectivity+0.1,mean([growth_sens_max;growth_sens_min]),Met2(group),'FontSize')%,8*2)
        xlabel('Conectivity')
        ylabel('\Delta\mu/\Deltap')
        % set(gca,'FontSize',8*2)
        
        subplot(1,2,2)
        plot(groupConectivity,groupCod,'*')
        title('Parameter codependence')
        text(groupConectivity+0.1,groupCod,Met2(group),'FontSize')%,8*2)
        xlabel('Conectivity')
        ylabel('\Delta(\Delta\mu/\Deltap)')
        x0=10;
        y0=10;
        width=1000;
        height=800;
        set(gcf,'units','points','position',[x0,y0,width,height])
        % set(gca,'FontSize',8*2)
        saveas(gcf,['parSenCodCorrConec_Ch_',groupname],'svg')
        
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
            
            save('3D_Sc.mat','groupConectivity','groupCost','groupSens','groupCod','metTypes','group','Met','growth_sens_vector');
        end
    end
end

