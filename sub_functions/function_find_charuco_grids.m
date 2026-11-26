function [x_cam, x_cam_points] = function_find_charuco_grids(pathname, filename,file_ext, maxWorkers)
% % 
% save temp_charuco
% moo=1
% return

% clear
% close all
% load temp_charuco

% 
% Parse video file

if strcmp(file_ext,'.cihx')
    video = quick_parse_cihx(fullfile(pathname, filename));
    num_im = video.TotalFrame;
    
else
    video = VideoReader(fullfile(pathname,filename));
    num_im=video.NumFrames;
end

f_start = 1;

% --- ChArUco PARAMETERS ---
% patternDims = [9 10]; 
patternDims = [14 14]; 

family = "DICT_5X5_1000";
checkerSize = 100; 
markerSize = 75;   

% Number of internal corners for a 9x10 square board is 8x9 = 72 corners
num_corners = (patternDims(1) - 1) * (patternDims(2) - 1);

% Get the physical world coordinates of the ideal board.
% This M-by-2 matrix matches the output of detectCharucoBoardPoints row-for-row.
worldPts = patternWorldPoints("charuco-board", patternDims, checkerSize);
world_X = worldPts(:, 1); % Physical columns
world_Y = worldPts(:, 2); % Physical rows

% Initialize outputs
x_cam = repmat({NaN(2, num_corners)}, 1, num_im);
x_cam_points = NaN(1, num_im);

warning off;

% for f = f_start
parfor (f = f_start:num_im, maxWorkers)

    if strcmp(file_ext,'.cihx')
        im = readmraw(video, f);
    else
        im=read(video,f);
    end

    if isa(im, 'uint16')
        im = uint8(double(im) / 16);
    end
    im = im(:, :, 1); 

    im=255-im;
    
    % Initialize default empty array for this frame
    frame_pts = NaN(2, num_corners);
    pts_count = NaN;
    
    % ---------------------------------------------------------
    % DETECT CHARUCO
    % imagePoints: A 72x2 matrix perfectly aligned to the world points. 
    % Missing corners are automatically padded with NaN.
    % ---------------------------------------------------------
    
    imagePoints = detectCharucoBoardPoints(im, patternDims, family, checkerSize, markerSize);

    if ~isempty(imagePoints)
            
        % Find the row indices of successfully detected points
        valid_mask = ~isnan(imagePoints(:, 1));
        
        % RULE 1: At least 9 points total (absolute minimum for 3 rows x 3 pts + 3 cols x 3 pts)
        if sum(valid_mask) >= 9
            
            % Extract the physical world X and Y of the detected points
            detected_world_X = world_X(valid_mask);
            detected_world_Y = world_Y(valid_mask);
            
            % Count how many points exist in each unique physical column (X)
            unique_X = unique(detected_world_X);
            pts_per_col = arrayfun(@(val) sum(detected_world_X == val), unique_X);
            valid_cols = sum(pts_per_col >= 3);
            
            % Count how many points exist in each unique physical row (Y)
            unique_Y = unique(detected_world_Y);
            pts_per_row = arrayfun(@(val) sum(detected_world_Y == val), unique_Y);
            valid_rows = sum(pts_per_row >= 3);
            
            % RULE 2 & 3: At least 3 columns must have 3+ points, AND 3 rows must have 3+ points
            if valid_cols >= 3 && valid_rows >= 3
                
                % Assign directly: it is already the correctly formatted 72x2 matrix
                frame_pts = imagePoints'; 
                pts_count = sum(valid_mask);
                
            end
        end
    end
    
    % Drop the local variables into the sliced parfor output array
    x_cam{f} = frame_pts;
    x_cam_points(f) = pts_count;
        
    warning on;

    %% just if doing plotting
    % 
    % num_rows=patternDims(1);
    % num_cols=patternDims(2);
    % 
    % I = generateCharucoBoard([num_rows num_cols]*100,patternDims,family,checkerSize,markerSize);
    % imagePoints = detectCharucoBoardPoints(I, patternDims, family, checkerSize, markerSize);
    % 
    % x_cam0=reshape(x_cam{f}(1,:),num_rows-1,num_cols-1);
    % y_cam0=reshape(x_cam{f}(2,:),num_rows-1,num_cols-1);
    % 
    % px0=reshape(imagePoints(:,1),num_rows-1,[]);
    % py0=reshape(imagePoints(:,2),num_rows-1,[]);
    % 
    % px00=px0;
    % py00=py0;
    % 
    % px00(isnan(x_cam0))=NaN;
    % py00(isnan(x_cam0))=NaN;
    % 
    % subplot(1,2,1)
    % hold off
    % imagesc2(im); hold on
    % plot(x_cam0,y_cam0,'.-','markersize',15,'linewidth',2)
    % % plot(px_mid2,py_mid2,'.r','markersize',15)
    % % plot(px_ctrl(1),py_ctrl(1),'.r','markersize',16);
    % % plot(px_ctrl(2),py_ctrl(2),'.g','markersize',16);
    % 
    % for i=1:numel(x_cam0)
    % 
    %     text(x_cam0(i)+10,y_cam0(i)-20,num2str(i),'Color','r')
    % 
    % end
    % 
    % title(num2str(f))
    % axis off
    % 
    % subplot(1,2,2)
    % imagesc2(I); hold on
    % plot(px00,py00,'.-','markersize',15,'linewidth',2)
    % 
    % axis off
    % 
    % set(gcf,'color','w','position',[100 100 900 450])
    % 
    % for i=1:numel(px00)
    % 
    %     text(px00(i)+.1,py00(i)-.2,num2str(i),'Color','r')
    % 
    % end
    % 
    % drawnow
    % return
end

warning on;
