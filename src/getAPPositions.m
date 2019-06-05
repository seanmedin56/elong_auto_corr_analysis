% finds center of the stripe and relative ap positions

%%% Smoothing Kernel
kernel_radius = 30; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 15; 
[x_ref_kernel, y_ref_kernel] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref_kernel = x_ref_kernel - kernel_radius - 1;
y_ref_kernel = y_ref_kernel - kernel_radius - 1;
r_mat = sqrt(x_ref_kernel.^2 + y_ref_kernel.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
g_kernel(r_mat>kernel_radius) = 0;
%%% fitting parameters
min_mat_time = 30*60; % Minimum time for appearance of mature stripe
xDim = 512; % FOV x size in pixels
yDim = 256;
search_kernel_width = 2.0; % Size of summation bin to use for optimization
new_trace_struct = []; % Store trace results 
new_nucleus_struct = []; % Store nucleus results
search_swath = 15; % Search space (degrees). Size of deviation from AP orthogonal permitted
search_struct = struct; % Structure to store results
set_ids = unique([nucleus_struct.setID]); 
for i = 1:length(set_ids)   
    set_trace_struct = nucleus_struct([nucleus_struct.setID] == set_ids(i));
    set_trace_struct = set_trace_struct(~isnan([set_trace_struct.ParticleID]));
    set_trace_struct = set_trace_struct([set_trace_struct.ncStart] == 14);
    CenterAngle = -round(set_trace_struct(1).APAngle); 
    
    % Allow for inferred center line to deviate by at most "search_swath"
    % degrees from perpendicular to AP axis
    theta_vec = (CenterAngle-search_swath):(CenterAngle+search_swath);    
    
    % Get position and production info
    xp_all = [set_trace_struct.xPosParticle]; % X positions
    xp_all = xp_all(~isnan(xp_all));
    yp_all = [set_trace_struct.yPosParticle]; % Y positions    
    yp_all = yp_all(~isnan(yp_all));
    fluo_all = [set_trace_struct.fluo]; % Fluorescence
    time_all = [set_trace_struct.time]; % Time
    time_all = time_all(~isnan(fluo_all));
    fluo_all = fluo_all(~isnan(fluo_all));
    ap_all = [set_trace_struct.apPosParticle]; % AP postions    
    
    % Filter for times later than specified maturation time. Also remove
    % obs with nonpositive fluorescence
    time_used = min(min_mat_time, max(time_all) - 600);
    filter = (time_all>time_used)&(fluo_all>0);
    fluo_all = fluo_all(filter);      
    time_all = time_all(filter);
    xp_all = xp_all(filter);
    yp_all = yp_all(filter);    
    ap_all = ap_all(filter);
    
    % Array to store Fluorescence values per pixel
    set_frame = zeros(yDim,xDim);   
    for j = 1:length(fluo_all) % Could vectorize this...
        row = round(yDim - yp_all(j) + 1);
        set_frame(row,ceil(xp_all(j))) ...
            = set_frame(row,ceil(xp_all(j))) + fluo_all(j);  
    end    
    norm_ref_array = ones(size(set_frame));
    norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
    gauss_array = conv2(set_frame,g_kernel,'same');
    gauss_array = gauss_array./norm_ref_array;    
    % Calculate AP-per-pixel calibration
    xMat = abs(repmat(xp_all,length(xp_all),1)' - repmat(xp_all,length(xp_all),1)).^2;
    yMat = abs(repmat(yp_all,length(yp_all),1)' - repmat(yp_all,length(yp_all),1)).^2;
    rMat = sqrt(yMat + xMat);
    apMat = abs(repmat(ap_all,length(ap_all),1)-repmat(ap_all,length(ap_all),1)');
    apMat = apMat(rMat>100); % Guard against potential anomalies with smal p separations
    rMat = rMat(rMat>100);
    APperPixel = max(apMat./rMat)*100; % Max value should occur when spot separation is orthogonal to AP axis

    % Calculate Fluorecence Profile for each search angle
    f_sums = []; % Track fluorescence
    t_vec = []; % Track Angles
    p_vec = []; % Projected Position
    index_vec = []; % Convenience vector for indexing
    for t = 1:length(theta_vec)
        theta = theta_vec(t);
        projection_mat = zeros(yDim,xDim); % Store projected positions
        for m = 1:size(projection_mat,1)
            for n = 1:size(projection_mat,2)                
                row = size(projection_mat,1) - m + 1; % Flip Y direction
                projection_mat(row,n) = round(cosd(atand(row/n)-theta)*(sqrt(row^2 + n^2)));
            end
        end        
        search_struct(t).p_mat = projection_mat;
        search_struct(t).theta = theta;
        % Get unique projection values 
        unique_p = unique(projection_mat);        
        % We will use mean fluo per pixel for each position along
        % projection axis
        projection_means = zeros(1,length(unique(projection_mat)));        
        for p = 1:length(unique_p)
            projection_means(p) = mean(gauss_array(projection_mat==unique_p(p)));            
        end
        projection_means = projection_means/sum(projection_means); % Normalize
        search_struct(t).p_means = projection_means;
        search_struct(t).p_index = unique_p;
  
        % Find total share of fluorescence captured by prescribed window size 
        % centered at each point along projection axis
        for o = 1:length(unique_p)                        
            center = unique_p(o);
            w = search_kernel_width/APperPixel;
            f_sums = [f_sums sum(projection_means((unique_p>=(center-w))&(unique_p<=(center+w))))];            
            t_vec = [t_vec theta];
            p_vec = [p_vec center];
            index_vec = [index_vec t];
        end       
    end    
    %find best radius from mean set after using screening for well
    %populated orientations
    best_f = max(f_sums);  
    if sum(f_sums==best_f) > 1
        warning('Degenerate Solutions...Taking Median Index')
    end
    candidate_t = t_vec(f_sums==best_f);
    candidate_p = p_vec(f_sums==best_f);
    min_t = min(candidate_t);
    max_t = max(candidate_t);    
    candidate_t = sort(candidate_t+90)-90;
    %If there are degenerate solutions, take median angle
    med_index = ceil(length(candidate_t)/2);
    best_angle = candidate_t(med_index);
    best_center = candidate_p(candidate_t==best_angle);
    best_index = unique(index_vec(t_vec==best_angle));
    if length(best_center) > 1 || length(best_index) > 1
        error('Degenerate Centers or Indices')
    end       
    best_projection_mat = search_struct(best_index).p_mat-best_center;
    plot_mat = abs(best_projection_mat*APperPixel);
    plot_mat(plot_mat>8) = Inf;
    PixelperAP = APperPixel.^-1; 
    unique_p = sort(p_vec(t_vec==best_angle));
    %%% Make Figure      
    stripe_fig = figure('Position',[100 100 1024 512]);%,'Visible','off');        
    map = flipud(jet(64));
    colormap(map)    
    hold on    
    norm_array = 8 - 4*gauss_array/max(gauss_array(:));     
    im = imagesc(plot_mat); 
    h = colorbar;
    set(im,'AlphaData',.5);  
%     im2 = imagesc(Rnorm_array);
%     set(im2,'AlphaData',.8);  
%     colormap(parula(64))
    plot_times = time_all - min_mat_time;
    ss = scatter(xp_all,yDim - yp_all + 1,(fluo_all+15)/15,8-plot_times/max(plot_times)*4,'filled');    
    set(ss,'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    axis([0 xDim 0 yDim]);
    grid on
    title(['Estimated Stripe Position with Fluorescence Emissions, Set ' num2str(set_ids(i))]);
    xlabel('X Dim (Pixels)')
    ylabel('Y Dim (Pixels)')    
    ylabel(h,'distance from stripe (AP)')    
    
    %Assign Relative AP Position to Each Particle
    for j = 1:length(set_trace_struct)
        x_vec = set_trace_struct(j).xPosParticle;
        x_vec = x_vec(~isnan(x_vec));
        y_vec = yDim - set_trace_struct(j).yPosParticle + 1;
        y_vec = y_vec(~isnan(y_vec));
        ap_pos_vec = [];
        for k = 1:length(x_vec)
            ap_pos_vec = [ap_pos_vec APperPixel*best_projection_mat(y_vec(k),x_vec(k))];
        end
        set_trace_struct(j).rel_ap_vector = ap_pos_vec;
        set_trace_struct(j).MeanAP = mean(ap_pos_vec);
    end    
    % Update set struct with search info
    for j = 1:length(set_trace_struct)                
        set_trace_struct(j).PixelperAP = PixelperAP;        
        set_trace_struct(j).stripe_angle = best_angle + 90;
        set_trace_struct(j).deviation = CenterAngle-best_angle;
        set_trace_struct(j).projection_center = best_center;
        set_trace_struct(j).search_radius = search_kernel_width;        
%         set_trace_struct(j).MeanAP = mean(set_trace_struct(j).ap_vector);               
        set_trace_struct(j).search_swath = search_swath;
    end
    new_trace_struct = [new_trace_struct set_trace_struct]; 
    disp(['Completed ' num2str(i) ' of ' num2str(length(set_ids))])
end     
nucleus_struct = new_trace_struct;
%save('../dat/3_4kb_final/nucleus_struct2.mat', 'nucleus_struct');