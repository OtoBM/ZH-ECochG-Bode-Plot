%% Define insertion angle and tonotopic region of cochlear implant electrodes

% This MATLAB script is designed to define the insertion angle and
% tonotopic region of cochlear implant electrodes, facilitating
% visualization of electrocochleography (ECochG) data using the ZH-ECochG
% bode plot. It has been tested using MATLAB R2022b.
% 
% The algorithm embedded within this script allows for precise mapping of
% recording electrode placement relative to the cochlea's tonotopic
% organization, aiding the interpretation of ECochG responses.
% 
% % If you use this script, please cite the following publication: Geys,
% M.; Sijgers, L.; Dobrev, I.; Dalbert, A.; Röösli, C.; Pfiffner, F.;
% Huber, A. ZH-ECochG Bode Plot: A Novel Approach to Visualize
% Electrocochleographic Data in Cochlear Implant Users. J. Clin. Med. 2024,
% 13, 3470. https://doi.org/10.3390/jcm13123470
% 
% We encourage users to read the publication to understand the context,
% implementation, and potential applications of the methods employed in
% this script.
% 
% This script was developed by Ivo Dobrev and Leanne Sijgers, Department of
% Otolaryngology, Head & Neck Surgery, University Hospital Zurich,
% University of Zurich, Switzerland

clear all; close all; clc;

% calls the script 'parameters.m' to ask for user input and load parameters
parameters

if individualized_axes

    disp(' ')
    disp('Select the part of the image that is of interest, then right click and select "crop image"')
    img_raw1 = imread(axes_params.image_name);
    
    
    %% Crop image
    
    figure('name', 'Image cropping', 'units', 'normalized', ...
        'outerposition', [0 0.05 1 0.95], 'PaperPositionMode', 'auto', ...
        'Color', 'w')
    img_raw = imcrop(img_raw1); % ! right mouse click, "Crop Image"
    
    
    %% Show filtered image
    
    %convert top gray scale
    img_gray = rgb2gray(img_raw);
    
    %filter image - only contrast adjustment for now  
    Img_filt = img_gray;
    
    Img_filt = mag2db(  double(Img_filt)  );
    
    
    [y_n, x_n] = size(img_gray);
    x_vec  = axes_params.pixel_size.*((1:x_n) -1 );
    y_vec  = axes_params.pixel_size.*((1:y_n) -1 );
    
    figure
    subplot(1,2,1)
    imagesc(x_vec,y_vec,img_raw);
    axis image
    title('Original')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap bone
    
    subplot(1,2,2)
    imagesc(x_vec,y_vec,Img_filt);
    title('Filtered')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    axis image
    colormap jet
    
    % use a filtered image 
    img =  Img_filt;
    
    %%  Points selection
    
    f = figure('name','CI electrode location','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc( img );
    colormap turbo
    % xlabel('X (mm)')
    % ylabel('Y (mm)')
    axis image
    
    
    
    % Initialize an empty matrix to store the coordinates of the special points
    special_points = [];
    
    % Allow the user to select 2 special points using the mouse
    hold on;
    for i = 1:2
        if i == 1
            title('Select RW')
            
        else
             title('Select cochlear center')
        end
        
        
        [x, y] = ginput(1);
        plot(x, y, 'mo','markersize',20,'linewidth',2); % plot special points in green
        plot(x, y, 'm+','markersize',10,'linewidth',1); % plot special points in green
        special_points = [special_points; x y];
    end
    
    % Initialize an empty matrix to store the coordinates of the control points
    control_points = [];

    disp(' ')
    disp(['Select points along the electrode array, starting with the apical electrode and ending with the basal electrode.' ...
         ' Apart from the first and last electrode, the points selected along the array do not need to be exactly at the electrode locations, ' ...
         'and the number of selected points does not need to match the number of electrodes.'])
    
    % Allow the user to select control points using the mouse
     title(['Selection of electrode points, starting with the apical electrode and ending with the basal electrode' ...
         ' - left-click selects/ right-click removes last selection/ ENTER ends selection']);
    while true
       
        [x, y, button] = ginput(1);
        if isempty(button)
            break;
        end
        if button == 1 % left click
    %         plot(x, y, 'w+'); % plot control points
            plot(x, y, 'wo'); % plot control points
            control_points = [control_points; x y];
        elseif button == 3 % right click
            if ~isempty(control_points)
                plot(control_points(end,1), control_points(end,2), 'k*'); % remove marker
                control_points(end, :) = []; % remove last control point
            end
        end
    end
    
    % Ask the user if they want to correct the points
    choice = questdlg('Do you want to correct the points?', 'Control Points Correction', 'Yes', 'No', 'No');
    
    
    while strcmp(choice, 'Yes')
        % Clear the figure and show the image again
        clf;
        imagesc(img);
        colormap jet
    %     xlabel('X (mm)')
    %     ylabel('Y (mm)')
        axis image
        hold on;
         title('New selection of electrode points - left-click selects/ right-click removes last selection/ ENTER ends selection');
        
        % Plot the special points
        for i = 1:size(special_points, 1)
                plot(x, y, 'mo','markersize',20,'linewidth',2); % plot special points in green
                plot(x, y, 'm+','markersize',10,'linewidth',1); % plot special points in green
        end
        
        % Allow the user to select control points again
        control_points = [];
        while true
            [x, y, button] = ginput(1);
            if isempty(button)
                break;
            end
            if button == 1 % left click
                plot(x, y, 'wo'); % plot control points
                control_points = [control_points; x y];
            elseif button == 3 % right click
                if ~isempty(control_points)
                    plot(control_points(end,1), control_points(end,2), 'k*'); % remove marker
                    control_points(end, :) = []; % remove last control point
                end
            end
        end
        
        % Ask the user again if they want to correct the points
        choice = questdlg('Do you want to correct the points again?', 'Control Points Correction', 'Yes', 'No', 'No');
    end
    
    
    %% Show selection results 
    
    % Draw a polyline that connects the control point
    
    Coch_center = (special_points(2,:)-1).*axes_params.pixel_size;
    RW_center = (special_points(1,:)-1).*axes_params.pixel_size;
    Electrode_points = (control_points-1).*axes_params.pixel_size;
    
    f = figure('name','CI electrode location','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc(x_vec,y_vec,img);
    axis image
    title('Original')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap jet
    axis image
    
    hold on
    
    % % Plot the special points
    % for i = 1:size(special_points, 1)
    plot(Coch_center(1), Coch_center(2), 'mo','markersize',20,'linewidth',2); % plot special points in green
    plot(Coch_center(1), Coch_center(2), 'm+','markersize',10,'linewidth',1); % plot special points in green
    % end
    
    plot(RW_center(1), RW_center(2), 'wo','markersize',20,'linewidth',2); % plot special points in green
    plot(RW_center(1), RW_center(2), 'w+','markersize',10,'linewidth',1); % plot special points in green
    text(RW_center(1), RW_center(2) ,[' --- RW'],'color','w','fontsize',20)
    
    if size(Electrode_points,1) > 1
        plot(Electrode_points(:,1), Electrode_points(:,2), 'b-', 'LineWidth', 2);
    end
    
    elec_n = size(Electrode_points,1);
    
    
    
    for i = 1:elec_n
        line([Coch_center(1) Electrode_points(i,1)],[Coch_center(2) Electrode_points(i,2)],'color','w')
        text(Electrode_points(i,1), Electrode_points(i,2) ,['\leftarrow',num2str(i)],'color','k','fontsize',16)
        
    end
    
    line([Coch_center(1) RW_center(1)],[Coch_center(2) RW_center(2)],'color','w','linewidth',3,'LineStyle',':')
    
    
    hold off;
    
    
    
    %% calculate insertion angle
    
    start_point = RW_center;
    
    % Array of points (each row represents a point [x, y])
    points = Electrode_points; % Replace with your points
    
    % Central point
    central_point = Coch_center; % Replace with your central point
    
    % Number of points
    num_points = size(points, 1);
    
    % Initialize array to store angles
    angles = zeros(num_points, 1);
    
    % This forms a line along the x-axis
    reference_vector = start_point- central_point;
        
    points = flip(points);    
    
    for i = 1:num_points
        % Get the vector from central point to current point
        vector = points(i, :) - central_point;
    
        % Get the vector from central point to a reference point
        % Calculate the angle
        if axes_params.side == string('right')
            angle_rad = atan2(vector(1), vector(2)) - atan2(reference_vector(1), reference_vector(2));
        elseif axes_params.side == string('left')
            angle_rad = atan2(vector(2), vector(1)) - atan2(reference_vector(2), reference_vector(1));
        end
    
        % Convert angle to degrees
        angles(i) = rad2deg(angle_rad);
    
        % Adjust for negative angles
        if angles(i) < 0
            angles(i) = angles(i) + 360;
        end
    end
    
    angles_unwrapped = rad2deg(unwrap(deg2rad(angles)));
    
    % Display the angles
    disp('Angles (in degrees) between each point and the central point along the x-axis:');
    disp(angles_unwrapped);
    
    
    len_array = insdepth(end)-insdepth(1);
    
    if angles_unwrapped(1) > 0 
        angles_linspace = linspace(angles_unwrapped(1), angles_unwrapped(end), nrEL);
    else
        angles_linspace = linspace(0, angles_unwrapped(end), nrEL);
    end
    
    
    %% extract insertion depth
    
    %in direction from inside to outside - last is RW location - low to high freq'
    electrode_points_re_RW = flipud(  [Electrode_points; RW_center] ).';  %flip the electrode order before doing the distance calcualtion
    
    % Calculate a distance between for each two consecutive points - then make the cumulative
    % distance
    elec_dept = [0, cumsum(sqrt(sum(diff( electrode_points_re_RW , [], 2).^2)))];
    
    
    electrode_vec = [ 0 1:elec_n];
    electrode_vec_rev = fliplr(electrode_vec);
    
    
    %%  interpolate over the control points for better estimation of the depth
    
    interp_points_n = 1e6; % adapted from 1000
    interp_method = 'spline';  %  'cubic'   'spline'
    
    elec_dept_interp = linspace(elec_dept(1),elec_dept(end),interp_points_n); % Parameter for evaluating the spline
    
    x = elec_dept;
    
    X = electrode_points_re_RW(1,:);
    Y = electrode_points_re_RW(2,:);
    
    
    electrode_points_re_RW_interp = NaN(2,interp_points_n);
    
    electrode_points_re_RW_interp(1,:) = interp1(x,electrode_points_re_RW(1,:),elec_dept_interp,interp_method);
    electrode_points_re_RW_interp(2,:) = interp1(x,electrode_points_re_RW(2,:),elec_dept_interp,interp_method);
    
    
    figure
    plot(electrode_points_re_RW(1,:),electrode_points_re_RW(2,:),'bo')
    hold 
    plot(electrode_points_re_RW_interp(1,:),electrode_points_re_RW_interp(2,:),'r-')
    grid on
    
    xlabel('X (mm)')
    ylabel('Y (mm)')
    axis equal
    axis ij
    
    
    if elec_dept_interp(end)-len_array > 0 
        elecs_depth = linspace(elec_dept_interp(end)-len_array, elec_dept_interp(end), nrEL);
    else
        elecs_depth = linspace(0, len_array, nrEL);
    end
    
    
    %% recalculate electrode depth based on the spline
    
    max_depth = axes_params.cochlear_length; % mm;

    f_vec = Greenwood(elec_dept, max_depth);
    
    f_vec2 = Greenwood(elecs_depth, max_depth);

    
    %% plot frequency against electrode depth
    
    figure
    plot(elec_dept,f_vec,'o-')
    hold on
    plot(elecs_depth,f_vec2,'-')
    xlabel('Depth [mm]')
    ylabel('Greenwood frequency')
    set(gca, 'YScale', 'log')
    grid on
    
    
    %% Show final results 
    
    
    f = figure('name','CI electrode location','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc(x_vec,y_vec,img_raw);
    axis image
    title('Final result')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap gray
    axis image
    
    hold on
    
    
    plot(Coch_center(1), Coch_center(2), 'mo','markersize',20,'linewidth',2); % plot special points in green
    plot(Coch_center(1), Coch_center(2), 'm+','markersize',10,'linewidth',1); % plot special points in green
    
    plot(RW_center(1), RW_center(2), 'yo','markersize',20,'linewidth',2); % plot special points in green
    plot(RW_center(1), RW_center(2), 'y+','markersize',10,'linewidth',1); % plot special points in green
        
    text(RW_center(1), RW_center(2) ,[' --- RW'],'color','y','fontsize',20)
    
    
    if size(Electrode_points,1) > 1
        plot(Electrode_points(:,1), Electrode_points(:,2), 'b-', 'LineWidth', 2);
    end
    
    elec_n = size(Electrode_points,1);
    
    
    
    for i = 1:elec_n
        line([Coch_center(1) Electrode_points(i,1)],[Coch_center(2) Electrode_points(i,2)],'color','y')
        text(Electrode_points(i,1), Electrode_points(i,2) ,['\leftarrow',num2str(electrode_vec_rev(i))],'color','k','fontsize',18)
        
    end
    
    line([Coch_center(1) RW_center(1)],[Coch_center(2) RW_center(2)],'color','r','linewidth',3,'LineStyle',':')
    
    
    plot(electrode_points_re_RW_interp(1,:),electrode_points_re_RW_interp(2,:),'r-','linewidth',2)

    set(gca, 'FontSize', 20)
    
    
    hold off;

else 

    %% define axes to be saved 
    
    max_depth = axes_params.cochlear_length;
    elecs_depth = insdepth;
    angles_linspace = insangle;
    f_vec2 = Greenwood(elecs_depth, max_depth);

end

%% save sweep file

d = date;
t = timeofday(datetime);

file_out = strcat('axes_sweep_',  string(d), '_', string(t), '.xlsx');
file_out = strrep(file_out, ':', '-');

xlswrite(file_out, [string('Electrode'), 'Depth (mm)', 'Angle (degrees)', 'Frequency (Hz)'])
xlswrite(file_out, [ELs', elecs_depth', angles_linspace', f_vec2'], strcat('A2:D', num2str(length(ELs)+1)));

%% calculate monitoring axes

d_inter_EL = elecs_depth(2)-elecs_depth(1);
ang_inter_EL = angles_linspace(2)-angles_linspace(1);
ELs_at_base = floor(elecs_depth(1)/d_inter_EL);

elecs_depth_m = elecs_depth;
angles_linspace_m = angles_linspace;

ELs_m = [];

if ELs_at_base > 0
    for EL = 1:ELs_at_base
        elecs_depth_m = [elecs_depth(1)-d_inter_EL*EL, elecs_depth_m];
        angles_linspace_m = [angles_linspace(1)-ang_inter_EL*EL, angles_linspace_m];
        ELs_m = [ELs_m, string('full')];
    end
end

ELs_m = [string(1:length(ELs)), ELs_m];

f_vec2_m = Greenwood(elecs_depth_m, max_depth);

%% save monitoring file

file_out = strcat('axes_monitoring_',  string(d), '_', string(t), '.xlsx');
file_out = strrep(file_out, ':', '-');

xlswrite(file_out, [string('Nr. electrodes inserted'), 'Depth (mm)', 'Angle (degrees)', 'Frequency (Hz)'])
xlswrite(file_out, [ELs_m', elecs_depth_m', angles_linspace_m', f_vec2_m'], strcat('A2:D', num2str(length(ELs_m)+1)));

disp(' ')
disp('The excel files with the insertion angle and tonotopic frequency of the recording electrodes are now generated and saved in the folder.')

% exportgraphics(figure(7), 'output.pdf');

%% FUNCTIONS

function f = Greenwood(elec_depth, max_depth)
% https://en.wikipedia.org/wiki/Greenwood_function

    A_Greenwood = 165.4;
    a_Greenwood = 2.1; % if x is relative to the cochlea length
    K_Greenwood = 0.88;

    f = A_Greenwood.*(10.^(a_Greenwood.*(max_depth - elec_depth) ./max_depth)  - K_Greenwood);
    
end