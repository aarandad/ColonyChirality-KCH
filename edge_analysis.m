function output = edge_analysis(BF_TIF,YFP_TIF,CFP_TIF,saveto)
% This function takes Bright Field and Fluorescence Images of colonies and
% performs analysis of the edges. Its output is an array witht the traces
% of each border.
% The functions were built on the ones provided by Kirill Korolev.
% BF_TIF: Bright Field image filename
% YFP_TIF: YFP image filename
% CFP_TIF: CFP image filename
% saveto: folder to save output images
%%
% Set resolution factor (will decrease resolution if >1)
resolution_factor = 1;

% Open images and get only one channel, for bright field I used blue
IMBF = imread(BF_TIF);
if size(IMBF,3)>1
    IMBF = IMBF(:,:,3);
end
IMBF = imresize(IMBF,size(IMBF)/resolution_factor);

IMY = imread(YFP_TIF);
if size(IMY,3)>1
    IMY= IMY(:,:,2);
end
IMY = imresize(IMY,size(IMY)/resolution_factor);

IMC = imread(CFP_TIF);
if size(IMC,3)>1
    IMC = IMC(:,:,3);
end
IMC = imresize(IMC,size(IMC)/resolution_factor);

% Create composite image that will show the sectors over the bright field
Composite = IMBF-IMC+IMY;

% Make the images matrices
IMY = double(IMY);
IMC = double(IMC);

%% Find edges
edges_IM=edge(imgaussfilt(IMY-IMC,10),'log',0.05*resolution_factor); % CHANGE THIS IF NEEDED
edges_IM(1:15,:)=0;
edges_IM(end-15:end,:)=0;
edges_IM(:,1:15)=0;
edges_IM(:,end-15:end)=0;
%% Find center: selecting 20 points in the halo in the center of the colony or in the colony border.

gotcenter = false;
while ~gotcenter
    figure
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.5 0.8])
    imagesc(Composite)
    colormap gray
    axis equal
    xlabel('X-axis')
    ylabel('Y-axis')
    title('Edges of colony')
    warning('Pick edges of center of colony')
    coords = ginput(20);
    % Fit a circle with the coordinates
    circle_coords = CircleFitByPratt(coords);
    viscircles(circle_coords(1:2),circle_coords(3));
    hold on
    plot(circle_coords(1),circle_coords(2),'rx','MarkerSize',10,'linewidth',4)
    hold off
    choice = questdlg('Did you get the center?','Circle fitting','Yes','No','No');
    switch choice
        case 'Yes'
            gotcenter = true;
        case 'No'
            gotcenter = false;
    end
    close
end

%%
% Get polar coordinates of sector edges and fitted circle
[x_edges,y_edges,phi,r]=AAD_polar_coords(IMY,IMC,edges_IM,circle_coords,Composite,saveto,resolution_factor);
%%
% Filters the data to get rid of poionts in center of colony
prompt = inputdlg('Enter threshold for radius. ONLY numbers [For 300 enter nothing] :',...
    'Minimum Radius', [1 50]);
min_radius = str2num(prompt{:});
close
if isempty(min_radius)
    min_radius = 300;
end

max_radius = max(r);
phi_clean = phi(r>=min_radius&r<=max_radius);
r_clean = r(r>=min_radius&r<=max_radius);

% Sorts clean data by radius
[r_sort,idx_temp] = sort(r_clean,1);
p_sort = phi_clean(idx_temp);

%%
% Get traces
% Max_distance from one point in a trace to the next
max_distance = 10;
% Minimum length (in radius) a trace has to have to be considered trace
min_len = 100/resolution_factor;
%%
% Find traces
traces = AAD_connect_dots(p_sort,r_sort,max_distance,min_len);
%%
traces_original = traces;
colors_traces = plotSpread_distinguishable_colors(length(traces));
%%
% Plot colony with traces on top
AAD_figure('landscape')
imagesc(Composite)
colormap gray
axis equal
hold on
plot(x_edges,y_edges,'y.','Markersize',1);
for m = 1: length(traces)
    p_ = traces{m}(:,1);
    r_ = traces{m}(:,2);
    x_ = r_.*cos(p_)+circle_coords(1);
    y_ = r_.*sin(p_)+circle_coords(2);
    plot(x_,y_,'color',colors_traces(m,:),'Linewidth',0.05)
    text(x_(ceil(length(x_)/2))+10,y_(ceil(length(x_)/2))+10,num2str(m),'color',colors_traces(m,:),'Linewidth',0.05)
end
hold off
saveas(gcf,strcat(saveto,'\','Traces_edges_overlaid'),'pdf')
close
%%
% Remove traces manually

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])

imagesc(Composite)
colormap gray
axis equal
hold on
middle_coords = nan(length(traces),2);
for m = 1: length(traces)
    p_ = traces{m}(:,1);
    r_ = traces{m}(:,2);
    x_ = r_.*cos(p_)+circle_coords(1);
    y_ = r_.*sin(p_)+circle_coords(2);
    plot(x_,y_,'color',colors_traces(m,:),'Linewidth',0.05)
    text(x_(ceil(length(x_)/2))+10,y_(ceil(length(x_)/2))+10,num2str(m),'color',colors_traces(m,:),'Linewidth',0.05)
    plot(x_(floor(length(x_)/2)),y_(floor(length(x_)/2)),'o','color',colors_traces(m,:),'markerfacecolor',colors_traces(m,:),'Linewidth',1)
    middle_coords(m,:)=[x_(floor(length(x_)/2)),y_(floor(length(x_)/2))];
end

choice = questdlg('Do you want to take any traces out?','Take traces out','Yes','No','No');
switch choice
    case 'Yes'
        takeout = true;
    case 'No'
        takeout = false;
end

if isempty(takeout)
    takeout = false;
end


plot(circle_coords(1),circle_coords(2),'xr','Markersize',10,'linewidth',4)
text(circle_coords(1),circle_coords(2),'Click here when done','Color','red','Fontsize',12)
hold off

if takeout
    coords_totakeout = [];
    prompt = 'Click in center when done:';
    stop_taking_out = false;
    while ~stop_taking_out
        coords_totakeout = ginput(1);
        if sqrt(sum((coords_totakeout-circle_coords(1:2)).^2))<0.5*circle_coords(3)
            stop_taking_out = true;
        else
            [~,takeout_index]=min(sqrt(sum((middle_coords - repmat(coords_totakeout,length(traces),1)).^2,2)));
            traces(takeout_index) =[];
            middle_coords(takeout_index,:)=[];
            colors_traces(takeout_index,:) = [];
            imagesc(Composite)
            colormap gray
            axis equal
            hold on
            for m = 1: length(traces)
                p_ = traces{m}(:,1);
                r_ = traces{m}(:,2);
                x_ = r_.*cos(p_)+circle_coords(1);
                y_ = r_.*sin(p_)+circle_coords(2);
                plot(x_,y_,'color',colors_traces(m,:),'Linewidth',0.05)
                text(x_(ceil(length(x_)/2))+10,y_(ceil(length(x_)/2))+10,num2str(m),'color',colors_traces(m,:),'Linewidth',0.05)
                plot(x_(floor(length(x_)/2)),y_(floor(length(x_)/2)),'o','color',colors_traces(m,:),'markerfacecolor',colors_traces(m,:),'Linewidth',1)
            end
            plot(circle_coords(1),circle_coords(2),'xr','Markersize',10,'linewidth',4)
            text(circle_coords(1),circle_coords(2),'Click here when done','Color','red','Fontsize',12)
            hold off
            
        end
    end
end
close
%%
% Plot traces and edges over image after removing
AAD_figure('landscape')
imagesc(Composite)
colormap gray
axis equal
hold on
plot(x_edges,y_edges,'y.','Markersize',1);
for m = 1: length(traces)
    p_ = traces{m}(:,1);
    r_ = traces{m}(:,2);
    x_ = r_.*cos(p_)+circle_coords(1);
    y_ = r_.*sin(p_)+circle_coords(2);
    plot(x_,y_,'color',colors_traces(m,:),'Linewidth',0.05)
    text(x_(ceil(length(x_)/2))+10,y_(ceil(length(x_)/2))+10,num2str(m),'color',colors_traces(m,:),'Linewidth',0.05)
end
hold off
saveas(gcf,strcat(saveto,'\','Traces_edges_overlaid_post_removal'),'pdf')
close
%%
% Plot traces over image after removing
AAD_figure('landscape')
imagesc(Composite)
colormap gray
axis equal
hold on
for m = 1: length(traces)
    p_ = traces{m}(:,1);
    r_ = traces{m}(:,2);
    x_ = r_.*cos(p_)+circle_coords(1);
    y_ = r_.*sin(p_)+circle_coords(2);
    plot(x_,y_,'color',colors_traces(m,:),'Linewidth',0.05)
    text(x_(ceil(length(x_)/2))+10,y_(ceil(length(x_)/2))+10,num2str(m),'color',colors_traces(m,:),'Linewidth',0.05)
end
hold off
saveas(gcf,strcat(saveto,'\','Traces_overlaid'),'pdf')
close
%%
if ~isempty(traces)
    % Corrects angle with initial angle, normalizes radius with initial radius
    % and interpolates to get averages over traces for a given radius
    %%
    phi_normalized = cell(length(traces),1);
    ratio_normalized = cell(length(traces),1);
    logratio_normalized = cell(length(traces),1);
    % Plots the corrected angle vs the log(r/ri)
    
    AAD_figure('landscape')
    hold on
    plot(0,0,'r','linewidth',0.05)
    plot(0,0,'r','linewidth',3)
    plot(0,0,'k','linewidth',5)
    for m = 1: length(traces)
        phi0 = median(traces{m}(1:5,1));
        phi_normalized{m} = traces{m}(:,1) - phi0;
        ratio0 = traces{m}(1,2);
        ratio_normalized{m} = traces{m}(:,2)./min_radius;%either r0 or min_radius
        logratio_normalized{m} = log(ratio_normalized{m});
        plot(logratio_normalized{m},phi_normalized{m},'color',colors_traces(m,:),'Linewidth',0.05)
    end
    %%
    %remove back steps and interpolate phi
    ratio_normalized_up = cell (length(traces),1);
    phi_normalized_up = cell (length(traces),1);
    logratio_normalized_interp = 0:0.001:max(cell2mat(cellfun(@max,logratio_normalized,'un',0)));
    phi_normalized_interp = nan(length(traces),length(logratio_normalized_interp));
    
    %%
    for m = 1 : length(traces)
        rn_up_temp = ratio_normalized{m};
        pn_up_temp = phi_normalized{m};
        while sum(diff(rn_up_temp)<=0)>0
            keep = logical([1;diff(rn_up_temp)>0 & rn_up_temp(2:end)>0]);
            rn_up_temp = rn_up_temp(keep);
            pn_up_temp = pn_up_temp(keep);
        end
        ratio_normalized_up{m} = rn_up_temp;
        phi_normalized_up{m} = pn_up_temp;
        if length(pn_up_temp)>10
            phi_normalized_interp(m,:) = interp1(log(rn_up_temp),pn_up_temp,logratio_normalized_interp);
        end
    end
    % calculate mean trace
    %%
    phi_normalized_interp_moved_to_mean=phi_normalized_interp;
    for m = 1 : length(traces)
        idx_init(m) = find(~isnan(phi_normalized_interp(m,:)),1,'first');
    end
    [idx_sorted,idx] = sort(idx_init);
    idx_sorted_unique = unique(idx_sorted);
    idx_sorted_unique = [1,idx_sorted_unique,length(logratio_normalized_interp)];
    %%
    phi_normalized_interp_mean=[];
    for i = 1:length(idx_sorted_unique)-1
        phi_normalized_interp_mean(1,idx_sorted_unique(i):idx_sorted_unique(i+1)) = ...
            nanmean(phi_normalized_interp(:,idx_sorted_unique(i):idx_sorted_unique(i+1)));
        traces_to_move_up = idx(idx_sorted_unique(i+1)==idx_sorted);
        for k = 1:length(traces_to_move_up)
            if idx_sorted_unique(i+1)>1
                phi_normalized_interp_moved_to_mean(traces_to_move_up(k),:) = phi_normalized_interp_moved_to_mean(traces_to_move_up(k),:) + phi_normalized_interp_mean(idx_sorted_unique(i+1)-1)-phi_normalized_interp(traces_to_move_up(k),idx_sorted_unique(i+1));
            end
        end
    end
    %%
    for counter = 1:length(traces)
        plot(logratio_normalized_interp,phi_normalized_interp(counter,:),'color',colors_traces(counter,:),'linewidth',3)
    end
    plot(logratio_normalized_interp,phi_normalized_interp_mean,'k','linewidth',5)
    hold off
    legend('Original traces', 'Interpolated traces', 'Mean trace')
    xlabel('log(r/ri)')
    ylabel('Phi-Phi_i')
    saveas(gcf,strcat(saveto,'\','Mean_trace'),'pdf')
    close
    %%
    
    %figure
    %set(gcf,'units','normalized','outerposition',[0 0 1 1])
    AAD_figure('landscape')
    hold on
    plot(0,0,'r','linewidth',1)
    plot(0,0,'k','linewidth',3)
    for counter = 1:length(traces)
        plot(logratio_normalized_interp,phi_normalized_interp(counter,:),'color',colors_traces(counter,:),'linewidth',1)
    end
    plot(logratio_normalized_interp,phi_normalized_interp_mean,'k','linewidth',3)
    hold off
    legend('Interpolated traces', 'Mean trace')
    xlabel('log(r/ri)')
    ylabel('Phi-Phi_i')
    saveas(gcf,strcat(saveto,'\','Mean_trace_nogray'),'pdf')
    close
    
    %%
    AAD_figure('landscape')
    hold on
    phi_normalized_interp_smooth = nan(length(traces),length(logratio_normalized_interp));
    for m =1:length(logratio_normalized)
        [lrn_sort,idx]=sort(logratio_normalized{m}+rand(size(logratio_normalized{m}))*1e-5);
        lrn_sorted{m} = lrn_sort;
        phi_sorted{m} = traces{m}(idx,1);
        phi_smooth{m} = smooth(lrn_sorted{m},phi_sorted{m},100);
        phi_normalized_smooth{m} = phi_smooth{m}-phi_smooth{m}(1);
        phi_normalized_interp_smooth(m,:) = interp1(lrn_sorted{m},phi_normalized_smooth{m},logratio_normalized_interp);
        phi_normalized_to_smooth_init{m}=traces{m}(:,1) -phi_smooth{m}(1);
        plot(logratio_normalized{m},phi_normalized_to_smooth_init{m},'color',[0.9 0.9 0.9])
        plot(lrn_sorted{m},phi_normalized_smooth{m},'color',colors_traces(m,:))
    end
    saveas(gcf,strcat(saveto,'\','Traces_smooth'),'pdf')
    close
    %%
    
    AAD_figure('landscape')
    hold on
    for m =1:length(logratio_normalized)
        plot(logratio_normalized{m},traces{m}(:,1),'.','color',colors_traces(m,:))
        plot(lrn_sorted{m},phi_smooth{m},'k')
        text(lrn_sorted{m}(end)+0.01,phi_smooth{m}(end),num2str(m),'color',colors_traces(m,:))
    end
    saveas(gcf,strcat(saveto,'\','Traces_smooth_non_normalized'),'pdf')
    close
    %%
    
    AAD_figure('landscape')
    hold on
    for m =1:length(logratio_normalized)
        plot(logratio_normalized{m},traces{m}(:,1),'.','color',colors_traces(m,:))
        plot(lrn_sorted{m},phi_smooth{m},'k')
    end
    saveas(gcf,strcat(saveto,'\','Traces_smooth_non_normalized_no_numbers'),'pdf')
    close
    %%
    
    
else
    logratio_normalized_interp = NaN;
    phi_normalized_interp = NaN;
    phi_normalized_interp_smooth = NaN;
end
log_min_radius = log(min_radius);


output = struct(...
    'traces',{traces},...
    'traces_original',{traces_original},...
    'log_ratio_normalized_interp',logratio_normalized_interp,...
    'phi_normalized_interp',{phi_normalized_interp},...
    'phi_normalized_interp_smooth',{phi_normalized_interp_smooth},...
    'phi_smooth',{phi_smooth},...
    'x_edges',x_edges,...
    'y_edges',y_edges,...
    'phi',phi,...
    'r',r,...
    'circle_coords',circle_coords,...
    'log_min_radius',log_min_radius,...
    'phi_normalized_interp_moved_to_mean',{phi_normalized_interp_moved_to_mean},...
    'phi_normalized_interp_mean',{phi_normalized_interp_mean});

end


function [x_edges,y_edges,phi,r]=AAD_polar_coords(IMY,IMC,edges_IM,circle_coords,Composite,saveto,resolution_factor)
%%
%This function obtains coordinates of sector boundaries in polar
%coordinates: p is the angle, r is the radius. They are sorted such that p
%increases from -pi to pi. (x0,y0) is the center of the circle.
x0 = circle_coords(1);
y0 = circle_coords(2);
radius = circle_coords(3);

%convert the image with boundaries into an array with the xy coordinates
%of the boundaries
[y_edges,x_edges]=ind2sub(size(edges_IM),find(edges_IM));

%plot the circle and the boundaries for a visual check
AAD_figure('landscape')
imagesc(IMC-IMY)
axis equal
colormap gray
hold on
plot(x_edges,y_edges,'r.','markersize',0.5);
axis equal;
hold off;
xlabel('x-axis')
ylabel('y-axis')
title('Image used for edging')
saveas(gcf,strcat(saveto,'\','Edges'),'pdf')
close
%%

AAD_figure('landscape')
imagesc(Composite)
colormap gray
axis equal
hold on
plot(x_edges,y_edges,'y.','Markersize',0.05);
axis equal;
plot(x0+radius*cos(0:0.01:2*pi),y0+radius*sin(0:0.01:2*pi),'r.');
plot(x0,y0,'rx','MarkerSize',10,'linewidth',4)
hold off;
xlabel('x-axis')
ylabel('y-axis')
title('Colony with edges and circle fitted to borders')
saveas(gcf,strcat(saveto,'\','Edges_circle_overlaid'),'pdf')
close
%%
%shift the origin to the center of the circle
%convert xy coordinates into polar coordinates using complex numbers

complex_number=(x_edges-x0)+(y_edges-y0)*1i;
phi=angle(complex_number);
radius=abs(complex_number);
%sort by angle and plot for a visual check
[phi,idx] = sort(phi,1);
r = radius(idx);
%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])
AAD_figure('landscape')
plot(phi,r,'.','markersize',0.5)
ylim([0 1000/resolution_factor])
xlabel('Phi angle')
ylabel('Radius')
title('Radius of edges as a function of polar angle')
saveas(gcf,strcat(saveto,'\','Polar_coords'),'pdf')

%%
end


function traces = AAD_connect_dots(p_sort,r_sort,max_distance,min_len)
%%
traces = cell(1);
j = 1;
% keeps track of used points
taken = zeros(size(r_sort));
taken(end) = true;
idx_next = length(r_sort);
idx_trace = length(r_sort);

h = waitbar(0/length(r_sort),'% points evaluated');
%%
while sum(~taken)>0
    %%
    % calculates all distances
    waitbar(sum(taken)/length(r_sort))
    distances = (r_sort(idx_next)^2 + r_sort.^2 - 2*r_sort(idx_next).*r_sort .* cos(p_sort - p_sort(idx_next))).^(1/2);
    % finds close points that are above the point that's being analyzed
    idx_close = find(distances < max_distance & ~taken);
    % from those, finds the closest in radius
    idx_closest = idx_close(distances(idx_close) == min(distances(idx_close)));
    
    if isempty(idx_closest)%
        
        taken(idx_next) = true;
        % if no more points, either "close" that trace or start over
        if abs(r_sort(idx_trace(1))-r_sort(idx_trace(end))) > min_len && length(idx_trace)>40
            traces{j} = flipud(traces{j});
            traces{j+1}=[];
            j = j+1;
        else
            traces{j} = [];
        end
        idx_next = find(max(r_sort(~taken))==r_sort,1,'last');
        idx_trace = idx_next;
        taken(idx_next)=true;
    else
        taken(idx_next) = true;
        traces{j}(end+1,:)=[p_sort(idx_next),r_sort(idx_next)];
        idx_next = idx_closest(1);
        idx_trace(end+1) = idx_next;
    end
    
    %%
end
if abs(r_sort(idx_trace(1))-r_sort(idx_trace(end))) > min_len && length(idx_trace)>40
    traces{j} = flipud(traces{j});
else
    traces{j} = [];
end
close(h)
traces(cellfun(@isempty,traces))=[];
%%
% Manually break traces if they go back and forth
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
trace_mid_coordinates = nan(length(traces),2);
colors = plotSpread_distinguishable_colors(length(traces));
for counter = 1:length(traces)
    if ~isempty(traces{counter})
        plot(traces{counter}(:,1),traces{counter}(:,2),'color',colors(counter,:))
        hold on
        plot(traces{counter}(ceil(size(traces{counter},1)/2),1),traces{counter}(ceil(size(traces{counter},1)/2),2),'x','color',colors(counter,:),'linewidth',3,'markersize',10)
        text(traces{counter}(ceil(size(traces{counter},1)/2),1)+0.05,traces{counter}(ceil(size(traces{counter},1)/2),2),num2str(counter),'color',colors(counter,:))
        trace_mid_coordinates(counter,:) = traces{counter}(ceil(size(traces{counter},1)/2),:);
    end
end
hold off

choice = questdlg('Break trace?','Breaking traces','Yes','No','No');
switch choice
    case 'Yes'
        break_trace = true;
    case 'No'
        break_trace = false;
        close
end


while break_trace == true
    %%
    title('Click close to x in the center of trace to be broken')
    coords_trace_to_break = ginput(1);
    close
    phi_to_break = coords_trace_to_break(1);
    r_to_break = coords_trace_to_break(2);
    phi_traces = trace_mid_coordinates(:,1);
    r_traces = trace_mid_coordinates(:,2);
    [~,trace_to_break] = ...
        min(...
        sqrt( r_to_break.^2 + r_traces.^2 - 2* r_to_break.*r_traces.*cos( phi_traces - phi_to_break ) )...
        );
    %%
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot(traces{trace_to_break}(:,1),traces{trace_to_break}(:,2))
    hold on
    prompt = inputdlg('Points to break:',...
        'Sample', [1 50]);
    points_to_break = str2num(prompt{:});
    %%
    coords_trace_to_break = ginput(points_to_break);
    close
    index_to_break=nan(1,points_to_break);
    %%
    for counter = 1:points_to_break
        phi_to_break = coords_trace_to_break(counter,1);
        r_to_break = coords_trace_to_break(counter,2);
        phi_trace = traces{trace_to_break}(:,1);
        r_trace = traces{trace_to_break}(:,2);
        [~,index_temp] = min(...
            sqrt( r_to_break.^2 + r_trace.^2 - 2* r_to_break.*r_trace.*cos( phi_trace - phi_to_break ) ) ...
            );
        index_to_break(counter) = index_temp;
    end
    index_to_break_sort = [1,sort(index_to_break),length(traces{trace_to_break})];
    %%
    trace_temp = traces{trace_to_break};
    if trace_temp(index_to_break_sort(1),2)<trace_temp(index_to_break_sort(2)-1,2)
        traces{trace_to_break} = trace_temp(index_to_break_sort(1):index_to_break_sort(2)-1,:);
    else
        traces{trace_to_break} = flipud(trace_temp(index_to_break_sort(1):index_to_break_sort(2)-1,:));
    end
    
    for counter = 1:points_to_break
        if trace_temp(index_to_break_sort(counter+1),2)<trace_temp(index_to_break_sort(counter+2)-1,2)
            traces{end+1} = trace_temp(index_to_break_sort(counter+1)+1:index_to_break_sort(counter+2)-1,:);
        else
            traces{end+1} = flipud(trace_temp(index_to_break_sort(counter+1)+1:index_to_break_sort(counter+2)-1,:));
        end
    end
    %%
    
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    trace_mid_coordinates = nan(length(traces),2);
    colors = plotSpread_distinguishable_colors(length(traces));
    for counter = 1:length(traces)
        if ~isempty(traces{counter})
            plot(traces{counter}(:,1),traces{counter}(:,2),'color',colors(counter,:))
            hold on
            plot(traces{counter}(ceil(size(traces{counter},1)/2),1),traces{counter}(ceil(size(traces{counter},1)/2),2),'x','color',colors(counter,:),'linewidth',3,'markersize',10)
            trace_mid_coordinates(counter,:) = traces{counter}(ceil(size(traces{counter},1)/2),:);
        end
    end
    hold off
    
    choice = questdlg('Break another trace?','Breaking traces','Yes','No','No');
    switch choice
        case 'Yes'
            break_trace = true;
        case 'No'
            break_trace = false;
            close
    end
end

%%
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
trace_mid_coordinates = nan(length(traces),2);
colors = plotSpread_distinguishable_colors(length(traces));
for counter = 1:length(traces)
    if ~isempty(traces{counter})
        plot(traces{counter}(:,1),traces{counter}(:,2),'color',colors(counter,:))
        hold on
        plot(traces{counter}(ceil(size(traces{counter},1)/2),1),traces{counter}(ceil(size(traces{counter},1)/2),2),'x','color',colors(counter,:),'linewidth',3,'markersize',10)
        text(traces{counter}(ceil(size(traces{counter},1)/2),1)+0.05,traces{counter}(ceil(size(traces{counter},1)/2),2),num2str(counter),'color',colors(counter,:))
        trace_mid_coordinates(counter,:) = traces{counter}(ceil(size(traces{counter},1)/2),:);
    end
end
hold off

choice = questdlg('Stitch traces?','Stitching traces','Yes','No','No');
switch choice
    case 'Yes'
        stitch_trace = true;
    case 'No'
        stitch_trace = false;
        close
end
%%
while stitch_trace == true
    
    number_traces_to_stitch = 2;
    
    title('Click close to x in the center of trace to be stitched')
    coords_trace_to_stitch = ginput(number_traces_to_stitch);
    close
    phi_to_stitch = coords_trace_to_stitch(:,1);
    r_to_stitch = coords_trace_to_stitch(:,2);
    phi_traces = trace_mid_coordinates(:,1);
    r_traces = trace_mid_coordinates(:,2);
    
    traces_to_stitch = nan(number_traces_to_stitch,1);
    r_min_traces_to_stitch = nan(number_traces_to_stitch,1);
    for counter = 1:number_traces_to_stitch
        [~,temp_trace_to_stitch] = ...
            min(...
            sqrt( r_to_stitch(counter).^2 + r_traces.^2 - 2* r_to_stitch(counter).*r_traces.*cos( phi_traces - phi_to_stitch(counter) ) )...
            );
        traces_to_stitch(counter)=temp_trace_to_stitch;
        r_min_traces_to_stitch(counter) = traces{traces_to_stitch(counter)}(1,2);
    end
    
    [~,indexes_traces_to_stitch] = sort(r_min_traces_to_stitch);
    
    
    traces{traces_to_stitch(1)} = vertcat(traces{traces_to_stitch(indexes_traces_to_stitch)});
    traces(traces_to_stitch(2:end))=[];
    
    %%
    
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    trace_mid_coordinates = nan(length(traces),2);
    colors = plotSpread_distinguishable_colors(length(traces));
    for counter = 1:length(traces)
        if ~isempty(traces{counter})
            plot(traces{counter}(:,1),traces{counter}(:,2),'color',colors(counter,:))
            hold on
            plot(traces{counter}(ceil(size(traces{counter},1)/2),1),traces{counter}(ceil(size(traces{counter},1)/2),2),'x','color',colors(counter,:),'linewidth',3,'markersize',10)
            text(traces{counter}(ceil(size(traces{counter},1)/2),1)+0.05,traces{counter}(ceil(size(traces{counter},1)/2),2),num2str(counter),'color',colors(counter,:))
            trace_mid_coordinates(counter,:) = traces{counter}(ceil(size(traces{counter},1)/2),:);
        end
    end
    hold off
    
    choice = questdlg('Stitch another trace?','Stitching traces','Yes','No','No');
    switch choice
        case 'Yes'
            stitch_trace = true;
        case 'No'
            stitch_trace = false;
            close
    end
end


%% Fix the ones that go around

for counter = 1:length(traces)
    diff_phi = diff(traces{counter}(:,1));
    if any(abs(diff_phi)>(2*pi*0.95))
        indexes=find(abs(diff_phi)>(2*pi*0.95));
        indexes_right = indexes;
        indexes_right(1:2:end) = indexes(1:2:end)+1;
        indexes_right(end+1)=length(traces{counter}(:,1));
        for counter1 = 1:2:length(indexes)
            traces{counter}(indexes_right(counter1):indexes_right(counter1+1),1) = ...
                traces{counter}(indexes_right(counter1):indexes_right(counter1+1),1) ...
                - 2*pi * sign(diff_phi (indexes_right(counter1)-1));
        end
    end
end


end
