% function [x_cam, x_cam_points] = function_find_checkerboard_grids(pathname, filename, maxWorkers)
% save temp_function_find_checkerboard_grids

clear
close all
load temp_function_find_checkerboard_grids

% Parse video file
video = quick_parse_cihx(fullfile(pathname, filename));
n_im = video.TotalFrame;

% Board parameters based on your original script
boardSize_real = [8, 9]; % size of board in rows x columns
ctrl_marker_pos = [4, 4]; % position of white ctrl marker (rows x column)

% Initialize outputs
x_cam = repmat({NaN(2, boardSize_real(1) * boardSize_real(2))}, 1, n_im);
x_cam_points = NaN(1, n_im);

warning off;

f_start = 1;

% parfor (f = f_start:n_im, maxWorkers)
for f=f_start:n_im
    
    % Initialize default empty values for this specific frame
    x_cam0 = NaN(boardSize_real);
    y_cam0 = NaN(boardSize_real);
    x_cam_f = [x_cam0(:) y_cam0(:)]';
    x_cam_points_f = 0;
    
    im = readmraw(video, f);
    
    if isa(im, 'uint16')
        im = uint8(double(im) / 16);
    end
    % Ensure single channel if RGB
    im = im(:, :, 1); 
    
    [imagePoints, boardSize] = detectCheckerboardPoints(im, 'MinCornerMetric', 0.15);
 
    % imagesc2(im); hold on
    % plot(imagePoints(:,1),imagePoints(:,2),'-.g')
    % return

    if min(boardSize) > 4 && max(boardSize) <= max(boardSize_real + 1)
        
        % Trim points too close to the image edge
        imagePoints(imagePoints(:,1) < 10, :) = NaN;
        imagePoints(imagePoints(:,2) < 10, :) = NaN;
        imagePoints(imagePoints(:,1) > size(im, 2) - 9, :) = NaN;
        imagePoints(imagePoints(:,2) > size(im, 1) - 9, :) = NaN;
        
        if sum(~isnan(imagePoints(:,1))) > 20
            
            px = reshape(imagePoints(:,1), boardSize(1)-1, []);
            py = reshape(imagePoints(:,2), boardSize(1)-1, []);
            
            % imagesc2(im); hold on; plot(px,py,'.-g')
            % return
            % =========================================================
            % STRAIGHTNESS CHECK
            % =========================================================
            tolerance_pixels = 5.0; 
            max_deviation = 0;
            
            % Check Rows
            for r = 1:size(px, 1)
                validX = px(r, ~isnan(px(r, :)));
                validY = py(r, ~isnan(py(r, :)));
                if length(validX) > 2
                    p1 = [validX(1), validY(1)];
                    p2 = [validX(end), validY(end)];
                    line_len = norm(p2 - p1);
                    if line_len > 0
                        for k = 2:length(validX)-1
                            pk = [validX(k), validY(k)];
                            dist = abs((p2(1)-p1(1))*(p1(2)-pk(2)) - (p1(1)-pk(1))*(p2(2)-p1(2))) / line_len;
                            max_deviation = max(max_deviation, dist);
                        end
                    end
                end
            end
            
            % Check Columns
            for col = 1:size(px, 2)
                validX = px(~isnan(px(:, col)), col);
                validY = py(~isnan(py(:, col)), col);
                if length(validX) > 2
                    p1 = [validX(1), validY(1)];
                    p2 = [validX(end), validY(end)];
                    line_len = norm(p2 - p1);
                    if line_len > 0
                        for k = 2:length(validX)-1
                            pk = [validX(k), validY(k)];
                            dist = abs((p2(1)-p1(1))*(p1(2)-pk(2)) - (p1(1)-pk(1))*(p2(2)-p1(2))) / line_len;
                            max_deviation = max(max_deviation, dist);
                        end
                    end
                end
            end
            
            % Proceed only if lines are straight enough
            if max_deviation <= tolerance_pixels
                
                % get mid-points of each square
                px_mid = 0.5 * (px(:,1:end-1) + px(:,2:end));
                py_mid = 0.5 * (py(:,1:end-1) + py(:,2:end));
            
                px_mid2 = 0.5 * (px_mid(1:end-1,:) + px_mid(2:end,:));
                py_mid2 = 0.5 * (py_mid(1:end-1,:) + py_mid(2:end,:));
            
                px_mid2(isnan(py_mid2)) = NaN;
                py_mid2(isnan(px_mid2)) = NaN;
        
                % blur image
                imb = imgaussfilt(im, 5);
            
                % get values of square centres
                px_mid22 = px_mid2;
                py_mid22 = py_mid2;
                
                px_mid22(isnan(px_mid2)) = 1;
                py_mid22(isnan(py_mid2)) = 1;
        
                pind_mid2 = sub2ind(size(im), round(py_mid22), round(px_mid22));
                
                point_int = double(imb(pind_mid2));
                point_int(isnan(py_mid2)) = NaN;
        
                point_int = point_int - mean(point_int(:), 'omitnan');
            
                % create blank mesh
                BW01 = repmat([-1 1; 1 -1], ceil(0.5*size(point_int,1)), ceil(0.5*size(point_int,2)));
                BW01 = BW01(1:size(point_int,1), 1:size(point_int,2));
                
                R = corrcoef(point_int, BW01, 'rows', 'complete');
        
                % imagesc2(im); hold on;
                % plot(px,py,'.-g')
                % plot(px',py','.-g')
                % % plot(px_mid,py_mid,'.r')
                % return
                
                if abs(R(2)) > 0.5
                
                    if R(2) < 0
                        BW01 = BW01 * -1;
                    end
                
                    point_match = BW01 .* point_int;
                    ctrl_points = find(point_match < 0);
        
                    if length(ctrl_points) == 2
                
                        if point_int(ctrl_points(1)) < point_int(ctrl_points(2))
                            ctrl_points = flip(ctrl_points);
                        end
                    
                        [py_ctrl_ind, px_ctrl_ind] = ind2sub(size(point_match), ctrl_points);
                        
                        px_ctrl = px_mid2(ctrl_points);
                        py_ctrl = py_mid2(ctrl_points);
            
                        % map found points onto original grid
                        rotated = 0;
            
                        % rotate grid lines if running wrong way
                        if diff(py_ctrl_ind) ~= 0
                            px = px';
                            py = py';
                            
                            temp_offset = py_ctrl_ind;
                            py_ctrl_ind = px_ctrl_ind;
                            px_ctrl_ind = temp_offset;
                            
                            rotated = 1;
                        end
                        
                        if ~any(size(px)-boardSize_real > 0) && ~any(px_ctrl_ind > 5) && ...
                           ~any(size(px,2)-px_ctrl_ind > 5) && ~any(py_ctrl_ind > 4) && ...
                           ~any(size(py,1)-py_ctrl_ind > 4)
                            
                            is_valid_orientation = true;
                            
                            if diff(py_ctrl_ind) == 0
                                
                                if diff(px_ctrl_ind) == 1
                                    py_ctrl_offset = ctrl_marker_pos(1) - py_ctrl_ind(1) + 1;
                                    px_ctrl_offset = ctrl_marker_pos(2) - px_ctrl_ind(1) + 1;
                                    
                                    x_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1, px_ctrl_offset:px_ctrl_offset+size(px,2)-1) = px;
                                    y_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1, px_ctrl_offset:px_ctrl_offset+size(px,2)-1) = py;
                                    
                                else
                                    py_ctrl_offset = ctrl_marker_pos(1) - (size(px,1) - py_ctrl_ind(1)) + 1;
                                    px_ctrl_offset = ctrl_marker_pos(2) - (size(px,2) - px_ctrl_ind(1)) + 1;
                
                                    % rotate 180 deg
                                    px = flip(px, 2);
                                    py = flip(py, 2);
                                    px = flip(px, 1);
                                    py = flip(py, 1);
                                                        
                                    x_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1, px_ctrl_offset:px_ctrl_offset+size(px,2)-1) = px;
                                    y_cam0(py_ctrl_offset:py_ctrl_offset+size(px,1)-1, px_ctrl_offset:px_ctrl_offset+size(px,2)-1) = py;
                                end
                                
                            else
                                % Replaces the 'return' from original script to prevent loop crashing
                                is_valid_orientation = false; 
                            end
                            
                            
                            if is_valid_orientation
                                % ensure correct orientation
                                P0 = [px_ctrl(2), py_ctrl(2)];
                                P1 = [px_ctrl(1), py_ctrl(1)];
                                P2 = [px(:,1), py(:,1)];
                
                                angle4 = NaN(size(px,1), 1);
                
                                for i = 1:size(P2, 1)
                                    angle4(i) = angle3Points(P1, P0, P2(i,:));
                                end
                
                                angle4 = unwrap(angle4);
                
                                if mean(diff(angle4), 'omitnan') > 0
                                    x_cam0 = flip(x_cam0, 1);
                                    y_cam0 = flip(y_cam0, 1);
                                end
            
                                % Assign to local temporary variables for parfor slicing
                                x_cam_f = [x_cam0(:) y_cam0(:)]';
                                x_cam_points_f = sum(~isnan(x_cam0(:)));
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Drop the local variables into the sliced output array
    x_cam{f} = x_cam_f;
    x_cam_points(f) = x_cam_points_f;


    %% just if doing plotting
    num_rows=boardSize_real(1);
    num_cols=boardSize_real(2);
    
    [xc0, yc0]=pol2cart(linspace(0,2*pi),0.25);

    x_cam0=reshape(x_cam{f}(1,:),num_rows,num_cols);
    y_cam0=reshape(x_cam{f}(2,:),num_rows,num_cols);

    squareSize=1;
    worldPoints = squareSize+0.5+generateCheckerboardPoints([num_rows+1 num_cols+1], squareSize);
    I_check=(checkerboard(squareSize,ceil(0.5*(num_rows+1)),ceil(0.5*(num_cols+1)))>0.5);

    if mod(num_rows, 2) == 0
        I_check(end,:)=[];
    end

    if mod(num_cols, 2) == 0
        I_check(:,end)=[];
    end
        
    px0=reshape(worldPoints(:,1),num_rows,[]);
    py0=reshape(worldPoints(:,2),num_rows,[]);

    px00=px0;
    py00=py0;
    
    px00(isnan(x_cam0))=NaN;
    py00(isnan(x_cam0))=NaN;

    subplot(1,2,1)
    hold off
    imagesc2(im); hold on
    plot(x_cam0,y_cam0,'.-','markersize',15,'linewidth',2)
    % plot(px_mid2,py_mid2,'.r','markersize',15)
    % plot(px_ctrl(1),py_ctrl(1),'.r','markersize',16);
    % plot(px_ctrl(2),py_ctrl(2),'.g','markersize',16);

    for i=1:numel(x_cam0)

        text(x_cam0(i)+10,y_cam0(i)-20,num2str(i),'Color','r')

    end

    title(num2str(f))
    axis off

    subplot(1,2,2)
    imagesc2(I_check); hold on
    plot(px00,py00,'.-','markersize',15,'linewidth',2)
    patch(xc0+ctrl_marker_pos(1)+1,yc0+ctrl_marker_pos(2)+1,'w')
    patch(xc0+ctrl_marker_pos(1)+2,yc0+ctrl_marker_pos(2)+1,'k')
    axis off

    set(gcf,'color','w','position',[100 100 900 450])

    for i=1:numel(px00)

        text(px00(i)+.1,py00(i)-.2,num2str(i),'Color','r')

    end
    
    drawnow
    % return

end
