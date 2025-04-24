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
    
    f = figure('name','Cropping ROI selection','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    sgtitle('Select ROI, right mouse click, "Crop Image" ')
    img_raw = imcrop(img_raw1); % right mouse click, "Crop Image"
    

    
    %% Show filtered image
    
    % convert top gray scale
    img_gray = rgb2gray(img_raw);
    
    % filter image - only contrast adjustment
    Img_filt = img_gray;
    
    % protect against bad data 
    Img_filt(isinf(Img_filt)) = 1;
    Img_filt(Img_filt == 0) = 1;

    Img_filt = mag2db(  double(Img_filt)  );
    
    pix2mm =  axes_params.pixel_size;
    mm2pix = 1./pix2mm;
    
    f = figure('name','Filtered image','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    subplot(1,2,1)
    imagesc(img_raw);
    axis image
    title('Original')
    xlabel('X (pix)')
    ylabel('Y (pix)')
    colormap bone
    grid on
    
    subplot(1,2,2)
    imagesc(Img_filt);
    title('Filtered')
    xlabel('X (pix)')
    ylabel('Y (pix)')
    axis image
    colormap jet
    
    % use a filtered image
    img =  Img_filt;
    grid on
    
    

    %%  Selection of points for coordinate system - RW, Cochlear center, Electrode/Cochlear entrance
    
    f = figure('name','Selection of points for coordinate system','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc( img );
    
    try
        colormap turbo % this is avaialble only as of MATLAB R2020b - use jet for older versions
        disp('Current MATLAB version does not support color map "turbo". Using "jet" instead')
    catch
        colormap jet  % use jet for older versions
    end
    xlabel('X (pix)')
    ylabel('Y (pix)')
    axis image
    grid on
    
    
    % Initialize an empty matrix to store the coordinates of the special points
    special_points = [];
    
    Points_list = {'RW','Cochlear center','Electode entrance'};
    
    
    % Allow the user to select 3 special points using the mouse - RW, electrode entrance
    % and cochlear center
    hold on;
    for i = 1:numel(Points_list)
        
        title(['Select ',Points_list{i}]);
        
        [x, y] = ginput(1);
        plot(x, y, 'mo','markersize',20,'linewidth',2); % plot special points  
        plot(x, y, 'm+','markersize',10,'linewidth',1); % plot special points  
        special_points = [special_points; x y];
    end
    
    RW_pos_raw_pix = special_points(1,:);
    Coch_center_pos_raw_pix = special_points(2,:);
    Electrode_entrance_pos_raw_pix = special_points(3,:);
    

    
    %% Image transformation
     
    % Get image dimensions
    [h, w, ~] = size(img);
    
    % Compute translation amounts to move center_pt to image center
    img_center_x = w / 2;
    img_center_y = h / 2;
    translation_x = img_center_x - Coch_center_pos_raw_pix(1);
    translation_y = img_center_y - Coch_center_pos_raw_pix(2);
    
    % Create translation matrix
    translation_matrix = [1, 0, translation_x;
        0, 1, translation_y;
        0, 0, 1];
    


    %% Translation 

    % Apply translation to points
    Coch_center_pos_translated_pix = apply_transformation(Coch_center_pos_raw_pix, translation_matrix);
    RW_pos_translated_pix = apply_transformation(RW_pos_raw_pix, translation_matrix);
    Electrode_entrance_pos_translated_pix = apply_transformation(Electrode_entrance_pos_raw_pix, translation_matrix);
    
    % Translate the image
    translated_img = imtranslate(img, [translation_x, translation_y]);
    
    % Get the color range of the original image
    originalMin = min(img(:));
    originalMax = max(img(:));
    

    
    %% Rotation
    
    % Compute the angle for rotation 

    % center - RW
    dx = Coch_center_pos_translated_pix(1) - RW_pos_translated_pix(1);
    dy = Coch_center_pos_translated_pix(2) - RW_pos_translated_pix(2);
    
    if axes_params.side == string('right')
         rot_angle_img_deg = atan2d(dy, dx);  % Angle in degrees
    elseif axes_params.side == string('left')
        rot_angle_img_deg = atan2d(-dy, -dx);  % Angle in degrees
    end 

    % Rotate the image
    rotated_img = imrotate(translated_img, rot_angle_img_deg, 'bilinear', 'crop');  % ''bicubic''  bilinear
    
    % Create rotation matrix about the new center
    rot_angle_points_rad = -deg2rad(rot_angle_img_deg); % the rotation for carthesian system is opposite direction to the rotation for image coordite system becaue the image system is left hand
    
    % this is the rotation matirx but putiplied by the rotaion point : Rv = R*(point_pos_vec)
    rotation_matrix = [cos(rot_angle_points_rad), -sin(rot_angle_points_rad), img_center_x * (1 - cos(rot_angle_points_rad)) + img_center_y * sin(rot_angle_points_rad);
        sin(rot_angle_points_rad), cos(rot_angle_points_rad), img_center_y * (1 - cos(rot_angle_points_rad)) - img_center_x * sin(rot_angle_points_rad);
        0, 0, 1];
    
    % Apply rotation to second and third points
    RW_pos_rot_pix = apply_transformation(RW_pos_translated_pix, rotation_matrix);
    Electrode_entrance_pos_rot_pix = apply_transformation(Electrode_entrance_pos_translated_pix, rotation_matrix);
    Coch_center_pos_rot_pix = apply_transformation(Coch_center_pos_translated_pix, rotation_matrix); % Should  be  the same - just for checking
    


    %% Plotting
    
    marker_size = 10;
    
    % Show the final result
    f = figure('name','Image transformation','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    subplot(1,3,1)
    imagesc( ( img) );
    colormap jet
    axis image
    
    hold on
    plot(Coch_center_pos_raw_pix(1),Coch_center_pos_raw_pix(2),'+w','markersize',500)
    plot(img_center_x, img_center_y, 'wo', 'MarkerSize', marker_size, 'LineWidth', 2); % Center stays at image center
    plot(RW_pos_raw_pix(1), RW_pos_raw_pix(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(Electrode_entrance_pos_raw_pix(1), Electrode_entrance_pos_raw_pix(2), 'ko', 'MarkerSize', marker_size, 'LineWidth', 2);
    title('Raw position');
    xlabel('X (pix)')
    ylabel('Y (pix)')
    grid on
    
    
    subplot(1,3,2)
    imagesc( ( translated_img) );
    axis image
    colormap jet
    caxis([originalMin, originalMax]); % Set the color range
    hold on
    plot(Coch_center_pos_translated_pix(1),Coch_center_pos_translated_pix(2),'+w','markersize',500)
    plot(img_center_x, img_center_y, 'wo', 'MarkerSize', marker_size, 'LineWidth', 2); % Center stays at image center
    plot(RW_pos_translated_pix(1), RW_pos_translated_pix(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(Electrode_entrance_pos_translated_pix(1), Electrode_entrance_pos_translated_pix(2), 'ko', 'MarkerSize', marker_size, 'LineWidth', 2);
    title('Center of cochlea is translated to center of image');
    xlabel('X (pix)')
    ylabel('Y (pix)')
    grid on
    
        
    subplot(1,3,3)
    imagesc(rotated_img);
    hold on;
    plot(RW_pos_rot_pix(1), RW_pos_rot_pix(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(Electrode_entrance_pos_rot_pix(1), Electrode_entrance_pos_rot_pix(2), 'ko', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(img_center_x, img_center_y, 'wo', 'MarkerSize', marker_size, 'LineWidth', 2); % Center stays at image center
    plot(Coch_center_pos_translated_pix(1),Coch_center_pos_translated_pix(2),'+w','markersize',500)
    
    
    hold off;
    axis image
    xlabel('X (pix)')
    ylabel('Y (pix)')
    % end
    caxis([originalMin, originalMax]); % Set the color range
    colormap jet
    grid on
    
    title('Line of RW-Center is now horizontal');
    legend('RW','Entrance','Cochlear center','image center','Color', [1 1 1]*0.7)
    
    
    img_final = rotated_img;
    
    img_final(img_final < originalMin) = originalMin;
    img_final(img_final > originalMax) = originalMax;
    


    %%  Selection of points electrode array estimation - 1st and last electrode plus any point along the fit

    f = figure('name','CI electrode location','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    
    imagesc(img_final);
    hold on;
    plot(RW_pos_rot_pix(1), RW_pos_rot_pix(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(Electrode_entrance_pos_rot_pix(1), Electrode_entrance_pos_rot_pix(2), 'ko', 'MarkerSize', marker_size, 'LineWidth', 2);
    plot(img_center_x, img_center_y, 'wo', 'MarkerSize', marker_size, 'LineWidth', 2); % Center stays at image center
    plot(Coch_center_pos_translated_pix(1),Coch_center_pos_translated_pix(2),'+w','markersize',500)
    
    axis image
    xlabel('X (pix)')
    ylabel('Y (pix)')
    
    colormap jet
    grid on
    
    
    % Initialize an empty matrix to store the coordinates of the control points
    control_points = [];
    
    disp(' ') % the message was split into individual lines for better readability when the command window is not too wide
    disp(['Select points along the electrode array , starting with the apical electrode ']);
    disp(['and ending with the basal electrode. ']);
    disp([' Apart from the first and last electrode, the points selected along the array  ']);
    disp(['do not need to be exactly at the electrode pad locations, ']);
    disp(['and the number of selected points does not need to match the number of electrodes.']);
    
    % Allow the user to select control points using the mouse
    title(['Selection of electrode trajectory, starting with the apical electrode and ending with the basal electrode' ...
        ' - left-click selects/ right-click removes last selection/ ENTER ends selection']);
    
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
    
    % Ask the user if they want to correct the points
    choice = questdlg('Do you want to correct the points?', 'Control Points Correction', 'Yes', 'No', 'No');
    
    
    while strcmp(choice, 'Yes')
        % Clear the figure and show the image again
        clf;

        imagesc(img_final);
        hold on;
        plot(RW_pos_rot_pix(1), RW_pos_rot_pix(2), 'ro', 'MarkerSize', marker_size, 'LineWidth', 2);
        plot(Electrode_entrance_pos_rot_pix(1), Electrode_entrance_pos_rot_pix(2), 'ko', 'MarkerSize', marker_size, 'LineWidth', 2);
        plot(img_center_x, img_center_y, 'wo', 'MarkerSize', marker_size, 'LineWidth', 2); % Center stays at image center
        plot(Coch_center_pos_translated_pix(1),Coch_center_pos_translated_pix(2),'+w','markersize',500)
        
        axis image
        xlabel('X (pix)')
        ylabel('Y (pix)')
        title('New selection of electrode points - left-click selects/ right-click removes last selection/ ENTER ends selection');
        
        colormap jet
        grid on
        hold on
        
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



    %% Show selection results for electrode array and cochlea
    
    % Draw a polyline that connects the control points
    
    % convert from pixels to mm for all points of interest
    Coch_center_mm =Coch_center_pos_translated_pix.*axes_params.pixel_size;
    RW_center_mm = RW_pos_rot_pix.*axes_params.pixel_size;
    Electrode_entrance_mm = Electrode_entrance_pos_rot_pix.*axes_params.pixel_size;
    
    Electrode_points = (control_points-1).*axes_params.pixel_size;
    
    img_mm_x_vec = ((1:size(img_final,2)) -1).*axes_params.pixel_size;
    img_mm_y_vec = ((1:size(img_final,1)) -1).*axes_params.pixel_size;
    


    %% Selection of outer wall of cochlea  - only the section along the electrode array
    
    disp(['Select points along the outer cochlear wall, starting apically (approx. as deep as the electrodes) ']);
    disp(['and ending with the cochlear exit ']);
    
    f = figure('name','Cochlear outer wall','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc(img_mm_x_vec,img_mm_y_vec,img_final);
    
    grid on
 
    colormap jet  % - use jet for older versions of matlab
    
    axis image

    title('Outer wall selection')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    
    hold on
    
    % Cochlear center
    plot(Coch_center_mm(1), Coch_center_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
    plot(Coch_center_mm(1), Coch_center_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
    % end
    
    % RW
    plot(RW_center_mm(1), RW_center_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
    plot(RW_center_mm(1), RW_center_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
    text(RW_center_mm(1), RW_center_mm(2) ,['RW --'],'color','w','fontsize',20,'HorizontalAlignment', 'right')
    
    % electrode entrance
    plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
    plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
    text(Electrode_entrance_mm(1), Electrode_entrance_mm(2) ,['Entrance--'],'color','w','fontsize',20,'HorizontalAlignment', 'right')
    
    % electrode trajectory lines
    if size(Electrode_points,1) > 1
        plot(Electrode_points(:,1), Electrode_points(:,2), 'w-', 'LineWidth', 2);
    end
    
    elec_n = size(Electrode_points,1);
    
    
    % lines from the center to the electode trajectory fit
    for i = 1:elec_n
        line([Coch_center_mm(1) Electrode_points(i,1)],[Coch_center_mm(2) Electrode_points(i,2)],'color','w')        
    end

    % zero degree line
    line([Coch_center_mm(1) RW_center_mm(1)],[Coch_center_mm(2) RW_center_mm(2)],'color','w','linewidth',3,'LineStyle',':')
    
   
    % Initialize an empty matrix to store the coordinates of the control points
    Cochlear_out_wall = [];
    
    disp(' ') % the message was split into individual lines for better readability when the command window is not too wide
    disp(['Select points along the outside cochlear wall near the electrode array ,  ']);
    disp(['starting apically (approx. as deep as the electrodes) and ending with the entrance. ']);
    
    % Allow the user to select control points using the mouse
    title(['Selection of outer wall, starting  apically (at least as deep as the electrode) and ending with the entrace' ...
        ' - left-click selects/ right-click removes last selection/ ENTER ends selection']);
 
    while true
        
        [x, y, button] = ginput(1);
        if isempty(button)
            break;
        end
        if button == 1 % left click
            %         plot(x, y, 'w+'); % plot control points
            plot(x, y, 'wo'); % plot control points
            Cochlear_out_wall = [Cochlear_out_wall; x y];
        elseif button == 3 % right click
            if ~isempty(Cochlear_out_wall)
                plot(Cochlear_out_wall(end,1), Cochlear_out_wall(end,2), 'k*'); % remove marker
                Cochlear_out_wall(end, :) = []; % remove last control point
            end
        end
    end
    
    % Ask the user if they want to correct the points
    choice = questdlg('Do you want to correct the points?', 'Control Points Correction', 'Yes', 'No', 'No');
    
    
    while strcmp(choice, 'Yes')
        % Clear the figure and show the image again
        clf;
    
        f = figure('name','Cochlear outer wall','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
        imagesc(img_mm_x_vec,img_mm_y_vec,img_final);
        
        grid on
     
        colormap jet  % - use jet for older versions of matlab
        
        axis image
        title('Outer wall selection')
        xlabel('X (mm)')
        ylabel('Y (mm)')
        title('New selection of outer wall - left-click selects/ right-click removes last selection/ ENTER ends selection');
        %     colormap jet
        axis image
        
        hold on
        
        % Cochlear center
        plot(Coch_center_mm(1), Coch_center_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
        plot(Coch_center_mm(1), Coch_center_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
        % end
        
        % RW
        plot(RW_center_mm(1), RW_center_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
        plot(RW_center_mm(1), RW_center_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
        text(RW_center_mm(1), RW_center_mm(2) ,['RW --'],'color','w','fontsize',20,'HorizontalAlignment', 'right')
        
        % electrode entrance
        plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
        plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
        text(Electrode_entrance_mm(1), Electrode_entrance_mm(2) ,['Entrance--'],'color','w','fontsize',20,'HorizontalAlignment', 'right')
        
        % electrode trajectory lines
        if size(Electrode_points,1) > 1
            plot(Electrode_points(:,1), Electrode_points(:,2), 'w-', 'LineWidth', 2);
        end
        
        elec_n = size(Electrode_points,1);
        
        
        % lines from the center to the electode trajectory
        for i = 1:elec_n
            line([Coch_center_mm(1) Electrode_points(i,1)],[Coch_center_mm(2) Electrode_points(i,2)],'color','w')            
        end

        % zero degree line
        line([Coch_center_mm(1) RW_center_mm(1)],[Coch_center_mm(2) RW_center_mm(2)],'color','w','linewidth',3,'LineStyle',':')

        Cochlear_out_wall = [];

        while true
            
            [x, y, button] = ginput(1);
            if isempty(button)
                break;
            end
            if button == 1 % left click
                %         plot(x, y, 'w+'); % plot control points
                plot(x, y, 'wo'); % plot control points
                Cochlear_out_wall = [Cochlear_out_wall; x y];
            elseif button == 3 % right click
                if ~isempty(Cochlear_out_wall)
                    plot(Cochlear_out_wall(end,1), Cochlear_out_wall(end,2), 'k*'); % remove marker
                    Cochlear_out_wall(end, :) = []; % remove last control point
                end
            end
        end
        
        % Ask the user again if they want to correct the points
        choice = questdlg('Do you want to correct the points again?', 'Control Points Correction', 'Yes', 'No', 'No');
    end
    
    % electrode trajectory lines
    if size(Cochlear_out_wall,1) > 1
        plot(Cochlear_out_wall(:,1), Cochlear_out_wall(:,2), 'w-', 'LineWidth', 2);
    end
    


    %% Show selection results for electrode array and cochlea and outer wall of otic capsule
    
    f = figure('name','Selected points on the electrode array and outer wall','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc(img_mm_x_vec,img_mm_y_vec,img_final); 
    grid on
    axis image
    title('Selected points on the electrode array and outer wall')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap jet
    axis image
    
    hold on
    
    % Plot the special points
    plot(Coch_center_mm(1), Coch_center_mm(2), 'mo','markersize',20,'linewidth',2); % plot special points  
    plot(Coch_center_mm(1), Coch_center_mm(2), 'm+','markersize',10,'linewidth',1); % plot special points  
    
    plot(RW_center_mm(1), RW_center_mm(2), 'wo','markersize',20,'linewidth',2); % plot special points  
    plot(RW_center_mm(1), RW_center_mm(2), 'w+','markersize',10,'linewidth',1); % plot special points  
    text(RW_center_mm(1), RW_center_mm(2) ,['RW --'],'color','k','fontsize',20,'HorizontalAlignment', 'right')
    
    plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'ko','markersize',20,'linewidth',2); % plot special points  
    plot(Electrode_entrance_mm(1), Electrode_entrance_mm(2), 'k+','markersize',10,'linewidth',1); % plot special points  
    text(Electrode_entrance_mm(1), Electrode_entrance_mm(2) ,['Entrance--'],'color','k','fontsize',20,'HorizontalAlignment', 'right')
    
    elec_n = size(Electrode_points,1);
    
    
    
    for i = 1:elec_n
        line([Coch_center_mm(1) Electrode_points(i,1)],[Coch_center_mm(2) Electrode_points(i,2)],'color','w')
        text(Electrode_points(i,1), Electrode_points(i,2) ,['\leftarrow',num2str(i)],'color','k','fontsize',16)
        
    end
    
    %  Cochlear_out_wall
    plot(Cochlear_out_wall(:,1), Cochlear_out_wall(:,2), 'ws','markersize',20,'linewidth',2); % plot special points  
    
    line([Coch_center_mm(1) RW_center_mm(1)],[Coch_center_mm(2) RW_center_mm(2)],'color','w','linewidth',3,'LineStyle',':')
    
    
    hold off;



    %% take the outer wall XY coordinates and fit a polynomial in polar coordinates
    
    % Convert Cartesian coordinates to polar coordinates relative to the center
    Cochlear_out_wall_centered = Cochlear_out_wall - Coch_center_mm;
    
    CI_spline =  [Electrode_points];
    CI_spline_centered = CI_spline - Coch_center_mm;
    
    [wall_theta_wrapped, wall_r] = cart2pol(Cochlear_out_wall_centered(:, 1), ...
        Cochlear_out_wall_centered(:, 2));
    
    [CI_theta_wrapped, CI_r] = cart2pol(CI_spline_centered(:, 1), ...
        CI_spline_centered(:, 2));
    
    
    % unwrap the angle - very important ! - angle based on image horizonal
    % line - this us used so that we can plot back to the original image. Alternatively
    % we can rotate the image so that the RW-Cochlear_center line is horizontal
    wall_theta = unwrap(wall_theta_wrapped);
    CI_theta = unwrap(CI_theta_wrapped);
    
    % Fit a polynomial to the radius as a function of theta
    wall_poly_polar_n  = 6; % Degree of the polynomial, adjust as needed
    wall_poly_polar = polyfit(wall_theta, wall_r, wall_poly_polar_n);
    CI_poly_polar = polyfit(CI_theta, CI_r, wall_poly_polar_n);
    
    % Define the fitted spiral
    wall_r_fit = polyval(wall_poly_polar,wall_theta);
    CI_r_fit = polyval(CI_poly_polar,CI_theta);
    
    % convert to cartesian space purely for plotting purposes
    [wall_x_fit,wall_y_fit] = pol2cart( wall_theta,wall_r_fit);
    [CI_x_fit,CI_y_fit] = pol2cart( CI_theta,CI_r_fit);
    
    % calculate CI spline control points on the cochlear wall
    wall_r_fit_CI_points = polyval(wall_poly_polar,CI_theta);
    [wall_x_fit_CI_points,wall_y_fit_CI_points] = pol2cart( CI_theta,wall_r_fit_CI_points);
    

    % dummy for manufacturer specified electrode pad positions relative to the length of
    % the electrode starting from the apex
    electrode_pos_manufacturer_relative2apical_dummy = linspace(0, 1, nrEL);
    
    % intergate along the samples
    num_samples = 10000; % number of samples for numerical integration
    % use the apex angle from the CI spline, and the base angle from the cochlear wall,
    % but use the cochlear wall polynomial - then scale everything with 0.9 to account
    % for the discrepancy between the organ of corti and outer wall length
    % define theta and r with file resolution - it matters for the integration
    wall_theta_fine = linspace(CI_theta(1), CI_theta(end), num_samples).';
    wall_r_fit_fine = polyval(wall_poly_polar,wall_theta_fine);
    
    CI_theta_fine = linspace(CI_theta(1), CI_theta(end), num_samples).';
    CI_r_fit_fine = polyval(CI_poly_polar,CI_theta_fine);
    
    % calcualte the total length to be used as a refernce value for the relative lengths
    wall_length_total = trapz(wall_theta_fine,wall_r_fit_fine);
    
    CI_length_total = trapz(CI_theta_fine,CI_r_fit_fine);
    
    
    wall_length = zeros(num_samples,1);
    CI_length = zeros(num_samples,1);
    
    % doing cumulative sum across the length of the wall
    for i =2:num_samples
        wall_length(i) =  trapz(wall_theta_fine(1:i),wall_r_fit_fine(1:i));
        CI_length(i) =  trapz(CI_theta_fine(1:i),CI_r_fit_fine(1:i));
    end
    
    
    clear electrode_pos_r electrode_pos_theta electrode_length ...
        electrode_pos_r_CI_based  electrode_pos_theta_CI_based electrode_length_CI_based
    % find the length down the cochear wall for every relative position along the electrode
    for i = 1:numel(electrode_pos_manufacturer_relative2apical_dummy)
        
        % this will keep the length from apex to any specified electode pad - based on the electrode_pos_manufacturer_relative2apical_dummy
        % this could be subsituted with absolute positions of electrode pads based on the manufacturer spec
        cur_length = wall_length_total.*electrode_pos_manufacturer_relative2apical_dummy(i);
        cur_length_CI = CI_length_total.*electrode_pos_manufacturer_relative2apical_dummy(i);
        
        % cur_length_CI = wall_length_total.*electrode_pos_manufacturer_relative2apical_dummy(i);
        
        if cur_length >= 0
            cur_ind = find(cur_length <= wall_length, 1);
        else
            cur_ind = find(cur_length >= wall_length, 1);
        end
        
        if cur_length_CI >= 0
             cur_ind_CI = find(cur_length_CI <= CI_length, 1);
        else
             cur_ind_CI = find(cur_length_CI >= CI_length, 1);
        end

        electrode_pos_r(i) = wall_r_fit_fine(cur_ind);
        electrode_pos_theta(i) = wall_theta_fine(cur_ind);
        electrode_length(i) = wall_length(cur_ind); % scale this with a factor 0.9 to account for the offset between the OOC location and the cochlear wall 
        
        
        electrode_pos_r_CI_based(i) = CI_r_fit_fine(cur_ind_CI);
        electrode_pos_theta_CI_based(i) = CI_theta_fine(cur_ind_CI);
        electrode_length_CI_based(i) = CI_length(cur_ind_CI); % scale this with a factor 0.9 to account for the offset between the OOC locations and the cochlear wall 
        
        
        
    end
    
    % calculate electrode pad positions on the electrode spline - assuming the same angle
    % postions as on the outer wall
    electrode_pos_CI_r = polyval(CI_poly_polar,electrode_pos_theta);
    electrode_pos_CI_r_CI_based = polyval(CI_poly_polar,electrode_pos_theta_CI_based);
    
    
    % purely for plotting purposes
    [electrode_pos_x,electrode_pos_y] = pol2cart( electrode_pos_theta,electrode_pos_r);
    [electrode_pos_CI_x,electrode_pos_CI_y] = pol2cart( electrode_pos_theta,electrode_pos_CI_r);
    
    [electrode_pos_x_CI_based,electrode_pos_y_CI_based] = pol2cart( electrode_pos_theta_CI_based,electrode_pos_r_CI_based);
    [electrode_pos_CI_x_CI_based,electrode_pos_CI_y_CI_based] = pol2cart( electrode_pos_theta_CI_based,electrode_pos_CI_r_CI_based);
  
    image_name = axes_params.image_name; % e.g., 'S26.jpg'
    
    
    % Create the full filename for saving the plot
    file_out = strcat(patient_number, '_final.jpg');
    
    f = figure('name','CI and wall selected and interpolated points','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    imagesc(img_mm_x_vec-Coch_center_mm(1),img_mm_y_vec-Coch_center_mm(2),(img_final));
    grid on
    title('Selected points and interpolated points on the electrode array and outer wall')
    hold on
    plot(0,0,'+k','markersize',50, 'DisplayName', 'cochlear center');
    

    % special points
    plot(RW_center_mm(1) - Coch_center_mm(1), RW_center_mm(2) - Coch_center_mm(2), '*k','markersize',20, 'DisplayName', 'RW'); % plot special points  
  
    plot(Electrode_entrance_mm(1) - Coch_center_mm(1), Electrode_entrance_mm(2) - Coch_center_mm(2), 'dk','markersize',20, 'DisplayName', 'Electrode entrance'); % plot special points  

    plot(CI_spline_centered(:,1),CI_spline_centered(:,2), ...
          'o', 'Color', 'b', ...
          'DisplayName', 'selected points on electrode array');
    
%    plot(electrode_pos_CI_x,electrode_pos_CI_y, ...
%         'o', 'Color', 'b', ...
%         'markersize',20,'linewidth',2, ...
%          'DisplayName', 'electrode points on electrode array');

    hw = plot(wall_x_fit_CI_points,wall_y_fit_CI_points, ...
          'o', 'Color', 'k', ...
          'DisplayName', 'selected points on outer wall');
        
    plot(electrode_pos_x,electrode_pos_y, ...
        'o', 'Color', 'k', ...
        'markersize',20,'linewidth',2, ...
         'DisplayName', 'electrode points on outer wall');
        
    axis image
    colormap jet
    
    xlabel('X - mm')
    ylabel('Y - mm')
    
    % add legend
    legend show


    saveas(f, file_out);


    
    %% calculate insertion angle  
    
    ang_deg = CI_theta * (180 / pi);

    if axes_params.side == string('right')
    
        if any(ang_deg < 0)
            angles_unwrapped = 180 - ang_deg;
        else
            angles_unwrapped = 540 - ang_deg;
        end

    elseif axes_params.side == string('left')

        if any(ang_deg < 0)
            angles_unwrapped = ang_deg + 360;
        else
            angles_unwrapped = ang_deg
        end
    end

    % Display the angles
    disp(' ')
    disp('Angles (in degrees) between each point and the central point along the x-axis (RW-Center):');

    angles_unwrapped = flip(angles_unwrapped);

    disp(angles_unwrapped);
    
    
    len_array = insdepth(end)-insdepth(1);

    angles_linspace = linspace(angles_unwrapped(1), angles_unwrapped(end), nrEL);
   

    
    %% extract insertion depth based on electrode spline 
    
    % Evaluate fitted polynomial for upsampled points
    CI_theta_upsampled = linspace(CI_theta(1), CI_theta(end), 1e6);  % Define theta values 
    wall_r_fit_CI_points_upsampled = polyval(wall_poly_polar, CI_theta_upsampled);

    % Convert to Cartesian Coordinates
    [wall_x_fit_CI_points_u, wall_y_fit_CI_points_u] = pol2cart(CI_theta_upsampled, wall_r_fit_CI_points_upsampled);

    % Compute Arc Length (Sum of Euclidean segment distances)
    dx = diff(wall_x_fit_CI_points_u);
    dy = diff(wall_y_fit_CI_points_u);
    segment_lengths = sqrt(dx.^2 + dy.^2);
    arc_length_CI = sum(segment_lengths);

    % Display result
    fprintf('Approximated arc length of the fitted CI electrode array on the outer wall: %.4f\n', arc_length_CI);


    % compute angle between electrode insertion point and end of fitted wall 

    % coordinates at the end of the fitted wall
    x_w = wall_x_fit(end);
    y_w = wall_y_fit(end);
    
    % Compute the standard atan2 angle
    theta_wall = atan2(y_w, x_w); % Angle in radians
    
    % Adjust the angle based on the closest x-axis
    if x_w < 0
        theta_wall = pi - theta_wall;  % Measure from the negative x-axis
    end

    if theta_wall > pi
        theta_wall =  theta_wall - 2 * pi;
    end
   
    % Convert to degrees
    theta_wall_deg = rad2deg(theta_wall);
    
    % Display result
    fprintf('Angle of the electrode entrance point with the RW-cochlear center axis: %.2f radians (%.2f degrees)\n', theta_wall, theta_wall_deg);


    % compute angle for insertion point vs round window

    % compute angle between RW and cochleostomy 
    vector_RW_to_center = RW_center_mm - Coch_center_mm;
    vector_cochl_to_center = Electrode_entrance_mm - Coch_center_mm;

    % to make sure it works for both sides
    if axes_params.side == string('right')
        angle_coch = atan2(vector_cochl_to_center(1), vector_cochl_to_center(2)) - ...
            atan2(vector_RW_to_center(1), vector_RW_to_center(2));    
    elseif axes_params.side == string('left')
        angle_coch = atan2(vector_cochl_to_center(2), vector_cochl_to_center(1)) - ...
            atan2(vector_RW_to_center(2), vector_RW_to_center(1));   
    end


    % calculate total angle to be added to cochlear wall
    angle_to_add = theta_wall - angle_coch; % in radians



    % Evaluate fitted polynomial for upsampled points
%     basal_wall_theta_upsampled = linspace(angle_coch, CI_theta(1), 1e5);  % Define theta values 
    if ((CI_theta(end) > 0) && (wall_theta(end) > 0)) % && (angle_to_add > 0)
        basal_wall_theta_upsampled = linspace(CI_theta(end), wall_theta(end) + angle_to_add, 1e5);  % Define theta values 
    elseif ((CI_theta(end) < 0) && (wall_theta(end) < 0)) % && (angle_to_add > 0)
        basal_wall_theta_upsampled = linspace(CI_theta(end), wall_theta(end) - angle_to_add, 1e5);  % Define theta values 
    else
        disp('error with electrode insertion point angle calculation')
    end

    basal_wall_r_fit_upsampled = polyval(wall_poly_polar, basal_wall_theta_upsampled);

    % Convert to Cartesian Coordinates
    [basal_wall_x_fit_u, basal_wall_y_fit_u] = pol2cart(basal_wall_theta_upsampled, basal_wall_r_fit_upsampled);

    % Compute Arc Length (Sum of Euclidean segment distances)
    dx = diff(basal_wall_x_fit_u);
    dy = diff(basal_wall_y_fit_u);
    segment_lengths = sqrt(dx.^2 + dy.^2);
    arc_length_basal = sum(segment_lengths);

    % Display result
    fprintf('Approximated arc length of the basal part on the outer wall: %.4f\n', arc_length_basal);


    % plot result

    figure('name','Outer wall fit','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    plot(wall_x_fit_CI_points, wall_y_fit_CI_points, 'bo')
    hold all
    plot(wall_x_fit_CI_points_u, wall_y_fit_CI_points_u,'r-')
    plot(basal_wall_x_fit_u, basal_wall_y_fit_u,'g-')
    grid on

    legend('projected electrodes', 'wall fit', ['electrode entrance' newline 'to 1st projection'], 'Location', 'northwest')
    
    xlabel('X (mm)')
    ylabel('Y (mm)')
    title('Electrodes projected on outer wall fit')
    axis equal
    axis ij

    

    %%  interpolate over the control points across the electrode for better estimation of the depth 

    % calculate inter-electrode distances
    elecs_depth = linspace(arc_length_basal, arc_length_CI + arc_length_basal, nrEL);
    elecs_depth_ooc = elecs_depth * 0.9; % defined at organ of corti
    


    %% recalculate electrode depth based on the spline - TO BE MODIFFIED
    
    max_depth = axes_params.cochlear_length; % mm;
    
    f_vec2 = Greenwood(elecs_depth_ooc, max_depth); % uses estimated locations of electrodes at organ of corti
        

    
    %% plot frequency against electrode depth
    
    figure('name','Greenwood frequency against electrode depth (mm)','units','normalized','outerposition',[0 0.05 1 0.95],'PaperPositionMode','auto','Color','w');
    plot(elecs_depth,f_vec2,'-')
    xlabel('Depth [mm]')
    ylabel('Greenwood frequency')
    set(gca, 'YScale', 'log')
    grid on
    
    
else
    

    %% define axes to be saved
    
    max_depth = axes_params.cochlear_length;
    elecs_depth_ooc = insdepth;
    angles_linspace = insangle;
    f_vec2 = Greenwood(elecs_depth_ooc, max_depth);
    
end

%% save sweep file

d = date;
t = timeofday(datetime);

file_out = strcat(string(patient_number),'_final_axes_sweep_',  string(d), '_', string(t), '.xlsx');
file_out = strrep(file_out, ':', '-');

xlswrite(file_out, [string('Electrode'), 'Depth (mm)', 'Angle (degrees)', 'Frequency (Hz)'])
xlswrite(file_out, [ELs', elecs_depth_ooc', angles_linspace', f_vec2'], strcat('A2:D', num2str(length(ELs)+1)));



%% calculate monitoring axes

d_inter_EL = elecs_depth_ooc(2)-elecs_depth_ooc(1);
ang_inter_EL = angles_linspace(2)-angles_linspace(1);
ELs_at_base = floor(elecs_depth_ooc(1)/d_inter_EL);

elecs_depth_m = elecs_depth_ooc;
% angles_linspace = angles_linspace';
angles_linspace_m = angles_linspace;

ELs_m = [];

if ELs_at_base > 0
    for EL = 1:ELs_at_base
        elecs_depth_m = [elecs_depth_ooc(1)-d_inter_EL*EL, elecs_depth_m];
        angles_linspace_m = [angles_linspace(1)-ang_inter_EL*EL, angles_linspace_m];
        ELs_m = [ELs_m, string('full')];
    end
end

ELs_m = [string(1:length(ELs)), ELs_m];

f_vec2_m = Greenwood(elecs_depth_m, max_depth);



%% save monitoring file

file_out = strcat(string(patient_number), '_final_axes_monitoring_',  string(d), '_', string(t), '.xlsx');
file_out = strrep(file_out, ':', '-');

xlswrite(file_out, [string('Nr. electrodes inserted'), 'Depth (mm)', 'Angle (degrees)', 'Frequency (Hz)'])
xlswrite(file_out, [ELs_m', elecs_depth_m', angles_linspace_m', f_vec2_m'], strcat('A2:D', num2str(length(ELs_m)+1)));

disp(' ')
disp('The excel files with the insertion angle and tonotopic frequency of the recording electrodes are now generated and saved in the folder.')



%% FUNCTIONS

function f = Greenwood(elec_depth, max_depth)
% https://en.wikipedia.org/wiki/Greenwood_function

A_Greenwood = 165.4;
a_Greenwood = 2.1; % if x is relative to the cochlea length
K_Greenwood = 0.88;

f = A_Greenwood.*(10.^(a_Greenwood.*(max_depth - elec_depth) ./max_depth)  - K_Greenwood);

end

function new_pt = apply_transformation(pt, transformation_matrix)
    % Convert point to homogeneous coordinates
    pt_homogeneous = [pt(1); pt(2); 1];
    
    % Apply transformation
    new_pt_homogeneous = transformation_matrix * pt_homogeneous;
    
    % Extract new coordinates
    new_pt = new_pt_homogeneous(1:2)';
end

function rotatedImage = rotateAroundCenter(grayImage, center, angle)
% ROTATEAROUNDCENTER Rotates a grayscale image around a specified center.
%   rotatedImage = ROTATEAROUNDCENTER(grayImage, center, angle) rotates the
%   input grayscale image grayImage around the specified center of rotation
%   by the specified angle (in degrees). The center remains in the same
%   location after rotation.

%   Inputs:
%       grayImage - Grayscale image matrix (2D array).
%       center    - [x, y] coordinates of the center of rotation.
%       angle     - Rotation angle in degrees (positive for counterclockwise).
%   Outputs:
%       rotatedImage - Rotated image with the same center location.

    % Ensure the image is grayscale
    if size(grayImage, 3) ~= 1
        error('Input image must be grayscale (2D matrix).');
    end

    % Extract center coordinates
    xCenter = center(1);
    yCenter = center(2);

    % Get the size of the original image
    [rows, cols] = size(grayImage);

    % Create padding to avoid cutting off image data during rotation
    padX = max(rows, cols) - rows;
    padY = max(rows, cols) - cols;
    paddedImage = padarray(grayImage, [padX, padY], 'both');
    
%     paddedImage =grayImage;

    % Calculate the translation to shift the custom center to the image center
    xShift = (size(paddedImage, 2) / 2) - xCenter;
    yShift = (size(paddedImage, 1) / 2) - yCenter;

    % Translate the image to center the desired rotation point
    translatedImage = imtranslate(paddedImage, [xShift, yShift], 'OutputView','full');  % 'OutputView','full'

    % Rotate the image
    rotatedPaddedImage = imrotate(translatedImage, angle, 'bilinear', 'crop');
      rotatedPaddedImage = imrotate(translatedImage, angle, 'bilinear', 'loose');

%     % Translate back to restore the original center location
    restoredImage = imtranslate(rotatedPaddedImage, [-xShift, -yShift]);

% restoredImage = rotatedPaddedImage;

% 
%     % Crop the image to the original size
%     rotatedImage = imcrop(restoredImage, [padY + 1, padX + 1, cols - 1, rows - 1]);

rotatedImage = rotatedPaddedImage;
end