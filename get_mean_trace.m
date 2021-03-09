
minr = 5.9; % set to desired minimum log(r) 
maxr = 6.2; % set to desired maximum log(r) 
%%

%%

load(FILENAME) % result file (output from edge_analysis, colonies saved as fields in results structure)
results1=results;
load(FILENAME) % result file (output from edge_analysis, colonies saved as fields in results structure) for a second condition
results2=results;
%%
for cc = 1:2
    results = eval(['results',num2str(cc)]);
    colonies = fieldnames(results);
    log_radius_forplot = 0:0.001:8; % creates vector to have the same x axis for all traces
    %%
    for counter = 1:length(colonies) % loop over colonies
        results_colony = getfield(results,colonies{counter}); % get results for colony
        minimum_log_radius(counter) = results_colony(1).log_min_radius; % get the minimum radius for all traces in that colony
        traces_in_colony(counter) = length(results_colony.traces); % get number of traces in colony
        index_init(counter) = find(log_radius_forplot>round(minimum_log_radius(counter)*1000)/1000,1,'first'); % finds position in new x axis for the trace
        r{counter} = results_colony.r; % gets the radius of each point (before tracing)
        phi{counter} = results_colony.phi; % gets the angle for each point (before tracing)
        for counter1 = 1:length(results_colony) % loop over traces within colony
            traces{counter}{counter1} =  results_colony(counter1).traces; % get traces (radius and angle)
        end
        colony_name{counter}= colonies{counter};
    end
    
    % Builds matrix with all traces mapped to the new x axis
    phi_normalized_interp_smooth = nan(sum(traces_in_colony),length(log_radius_forplot));
    indexes_traces = [0,cumsum(traces_in_colony)];
    for counter = 1:length(colonies)
        results_colony = getfield(results,colonies{counter});
        temp = results_colony.phi_normalized_interp_smooth;
        phi_normalized_interp_smooth(indexes_traces(counter)+1:indexes_traces(counter+1),index_init(counter):size(temp,2)+index_init(counter)-1)=temp;
    end
    numberoftraces = sum(~isnan(phi_normalized_interp_smooth));
    numberoftraces(numberoftraces ==0)=nan;
    
    for counter = 1:length(colonies)
        results_colony = getfield(results,colonies{counter});
        temp = results_colony.phi_normalized_interp_smooth;
        temp_interp=nan(size(log_radius_forplot));
        temp_interp(index_init(counter):size(temp,2)+index_init(counter)-1)=nanmean(temp);
        phi_diff_temp = diff(temp_interp,[],2);
        phi_diff_temp(isnan(phi_diff_temp))=0;
        crop=log_radius_forplot>minr & log_radius_forplot<maxr;
        crop=crop(1:end-1);
        lr_diff = diff(log_radius_forplot);
        slope_colony(counter)=nanmean(atand(phi_diff_temp(:,crop)./lr_diff(crop)));
    end
    %% Get mean by integrating or just taking the mean
    phi_diff = diff(phi_normalized_interp_smooth,[],2);
    phi_diff_mean = nanmean(phi_diff);
    phi_diff_mean_temp = phi_diff_mean;
    phi_diff_mean_temp(isnan(phi_diff_mean_temp))=0;
    phi_cumsum_mean=cumsum(phi_diff_mean_temp);
    phi_cumsum_mean(isnan(phi_diff_mean))=nan;
    
    phi_mean = nanmean(phi_normalized_interp_smooth);
    
    %% bootstrapping
    allmeans=nan(200,size(phi_normalized_interp_smooth,2));
    allslopes=nan(200,size(phi_normalized_interp_smooth,2)-1);
    numtrajs=size(phi_normalized_interp_smooth,1);
    for i=1:200
        sampled_phi_normalized_interp_smooth=phi_normalized_interp_smooth(randi(numtrajs,numtrajs,1),:);
        allmeans(i,:)=nanmean(sampled_phi_normalized_interp_smooth);
        allslopes(i,:)=atand(diff(allmeans(i,:),1,2)./diff(log_radius_forplot,1,2));
        
        sampled_phi_diff = diff(sampled_phi_normalized_interp_smooth,[],2);
        sampled_phi_diff_mean = nanmean(sampled_phi_diff);
        sampled_phi_diff_mean_temp = sampled_phi_diff_mean;
        sampled_phi_diff_mean_temp(isnan(sampled_phi_diff_mean_temp))=0;
        sampled_phi_cumsum_mean=cumsum(sampled_phi_diff_mean_temp);
        sampled_phi_cumsum_mean(isnan(sampled_phi_diff_mean))=nan;
        
        allmeans_cumsum(i,:)=sampled_phi_cumsum_mean;
        allslopes_cumsum(i,:)=atand(diff(allmeans_cumsum(i,:),1,2)./diff(log_radius_forplot(1:end-1),1,2));
        i
    end
    ALL_allslopes_cumsum{cc} = allslopes_cumsum;
    ALL_allmeans_cumsum{cc} = allmeans_cumsum;
    ALL_log_radius_forplot{cc} = log_radius_forplot;
    
    x = ALL_log_radius_forplot{cc}(1:end-1);
    y = mean(ALL_allmeans_cumsum{cc},1);
    x(isnan(y))=[];
    minr=x(1)+0.2;
    maxr=x(1)+0.5;
    
    %% get distribution of slopes
    slopedistr=zeros(1,size(allmeans,1));
    for i=1:200
        crop=log_radius_forplot>minr & log_radius_forplot<maxr;
        slopedistr(i)=nanmean(allslopes(i,crop));
        slopedistr_cumsum(i)=nanmean(allslopes_cumsum(i,crop));
    end
    ALL_slopedistr{cc} = slopedistr;
   
end
