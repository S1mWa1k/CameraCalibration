function [x_cam, x_cam_points]=function_find_dot_grids(pathname,filename,file_ext,maxWorkers,num_cols,num_rows,marker_grid_cols,marker_grid_rows)
% save temp_function_find_grids

% clear
% close all
% load temp_function_find_grids

if strcmp(file_ext,'.cihx')
    video = quick_parse_cihx(fullfile(pathname, filename));
    num_im = video.TotalFrame;
    
else
    video = VideoReader(fullfile(pathname,filename));
    num_im=video.NumFrames;
end

f_start=1;

marker_grid_cols=marker_grid_cols(:);
marker_grid_rows=marker_grid_rows(:);

% Position of marker points
% marker_grid_col = [8 11 11 8 9]';
% marker_grid_row = [8 8 11 11 8]'; 
% marker_grid_cols = [12 15 15 12 13]';
% marker_grid_rows = [12 12 15 15 12]'; 

x_cam=repmat({NaN(2,num_cols*num_rows)},1,num_im);
x_cam_points=NaN(1,num_im);

parfor (f=f_start:num_im,maxWorkers)
% for f=f_start:num_im

    if strcmp(file_ext,'.cihx')
        im = readmraw(video, f);
    else
        im=read(video,f);
    end
        
    if isa(im,'uint16')
        im=uint8(double(im)/16);
    end
    
    im=im(:,:,1);
             
    if mean(im(:))>1

%         level = graythresh(im(im<240));
        % turn into a binary image
        
%         im2=imbinarize(im,level);
        im2=imbinarize(im,'adaptive','ForegroundPolarity','dark','Sensitivity',0.5);
%         im22(im<10)=0;

        if mean(im(:))>10
            im2(im<10)=0;
        else
            im2(im<3)=0;
        end
%         
        % close
        % imagesc2(im);
        % return
   
        im3=edge(im2,'canny',0.9,1);
        
        % imagesc2(im3); return

        % label each edge found in the binary image
        imL = logical(im3);

        % Find centre and coordinates of each edge
        s = regionprops(imL,'Centroid','PixelList');

        c=round(cat(1,s.Centroid));

        if ~isempty(c)
      
            imD=bwdist(im2);

            c2=imD(sub2ind(size(im),c(:,2),c(:,1)));
            c3=find(c2<1.5);
        %         return
            im4=zeros(size(im3));

            % Separate any edges the are brighter in the centre than outside
            for i=c3'
                c=s(i).PixelList(:,1);
                r=s(i).PixelList(:,2);
                im3(sub2ind(size(im),r,c))=0;
                im4(sub2ind(size(im),r,c))=1;
            end

            % relabel remaining edges
            s = regionprops(im3,'Area');
            L=(cat(1,s.Area));
        
            % Remove small edges that are likely to be noise
            L(L<10)=[];
            % Find the most common edge length, should be calibration spots
            Lm=median(L);

            % Find edges that are likely to be calibration spots
            imLp=imL;

            imLp2=imfill(imLp,'holes');


            % Find how far away they are from each other and remove ones that
            % aren't near others - must near at least 3 others

            s = regionprops(imLp2,255-im, 'WeightedCentroid','Extrema','Circularity','Area');

            % s = regionprops(imLp2,'Centroid','Extrema','Circularity','Area');

            circ=(cat(1,s.Circularity));

            s(circ<0.5)=[];

            e=(cat(1,s.Extrema));
            % c=(cat(1,s.Centroid));
            c=(cat(1,s.WeightedCentroid));

            a=(cat(1,s.Area));

            % imagesc2(im); hold on
            % plot(c(:,1),c(:,2),'.')
            % return
            % 
            % close
            % return

            if ~isempty(c)

                ex=reshape(e(:,1),8,[]);
                ey=reshape(e(:,2),8,[]);

                exmin=min(ex);
                exmax=min(ex);
                eymin=min(ey);
                eymax=max(ey);

                nn=max([exmin<2;exmax>(size(im,2)-1);eymin<2;eymax>(size(im,1)-1)]);
                c(nn,:)=[];
                a(nn)=[];

                n=find(isnan(c(:,1)));
        %         s(n)=[];

        
                if size(c,1)<8
            %         fprintf(['\n Error C' num2str(cam) ' fr' num2str(ii) ': No calibration points not found \n'])
                else

                    c(n,:)=[]; 
                    a(n)=[];

                    calib_row=c(:,2);
                    calib_col=c(:,1);

                    d=NaN(length(calib_row),2);
                    d3=NaN(length(calib_row),2);

                    for i=1:length(calib_row)

                        dd=sqrt((calib_row(i)-calib_row).^2+(calib_col(i)-calib_col).^2);
                        dd(dd==0)=[];
                        dd=sort(dd);
                        d(i,:)=dd(1:2);
                        d3(i,:)=dd(3);
                    end

                    c_spacing=median(d(:));

                    d2=[find(d(:,2)>3.5*c_spacing)' find(d(:,2)<0.4*c_spacing)'];

                    calib_col(d2)=[];
                    calib_row(d2)=[];

                    d(d2,:)=[];
                    a(d2)=[];

                    c_area=median(a);
                    a2=[find(a>1.5*c_area)' find(a<0.75*c_area)'];

                    d(a2,:)=[];
                    % a(a2)=[];
                    calib_col(a2)=[];
                    calib_row(a2)=[];
                    
                    % close
                    % imagesc2(im); hold on
                    % plot(calib_col,calib_row,'.g'); return    


                    c_spacing2=median(d(:));

                    d3(d2,:)=[];
                    c_spacing=median(d3(:));

                    % close
                    % imagesc2(im); hold on
                    % plot(calib_col,calib_row,'.g'); return        

                    %% Now label edges that aren't calibration points
                    [imLm,~] = bwlabel(im4);

                    imLm2=imfill(imLm,'holes');
        
                    s = regionprops(imLm2,255-im,'WeightedCentroid','Circularity','Perimeter');
        
                    circ=(cat(1,s.Circularity));
        
                    s(circ<0.4)=[];

                    c=(cat(1,s.WeightedCentroid));

                    % imagesc2(im4); hold on; plot(c(:,1),c(:,2),'.r')
                    % return

                    D = pdist2(c, c);

                    close_neighbors = (D < 3) & ~eye(size(c,1));

                    neighbor_counts = sum(close_neighbors, 2);

                    s(neighbor_counts~=1)=[];
                    c(neighbor_counts~=1,:)=[];

                    LLp=(cat(1,s.Perimeter));

                    c(LLp<0.9*Lm | LLp>2.2*Lm,:)=[];

                    % imagesc2(im); hold on; plot(c(:,1),c(:,2),'.r','markersize',14)
                    % return
                    % Now check centroids and that each is part of a pair.
                   
                    if size(c,1)<5

                    else

                        marker_col=c(:,1);
                        marker_row=c(:,2);

                        % close
                        % imagesc2(im); hold on
                        % plot(calib_col,calib_row,'.r');
                        % plot(marker_col,marker_row,'.g');
                        % drawnow
                        % return
                    
                        D = pdist2([marker_col marker_row], [calib_col calib_row]);
                        
                        c_spacing_h=hypot(2*c_spacing,2*c_spacing2);

                        close_neighbors = (D < 1.5*c_spacing_h);
                        neighbor_counts = sum(close_neighbors, 2);
                        
                        marker_col(neighbor_counts<8)=[];
                        marker_row(neighbor_counts<8)=[];
                        
                        % close
                        % imagesc2(im); hold on
                        % plot(calib_col,calib_row,'.g');
                        % plot(marker_col,marker_row,'.r');
                        % plot(marker_col(1),marker_row(1),'og');
                        % drawnow
                        % return 

                        % From here there should be five markers found. If less or more
                        % then terminate                    
                        if length(marker_col)~=5
        %                     fprintf(['\n Error C' num2str(cam) ' fr' num2str(ii) ': Incorrect number of marker points found \n'])
                        else

                            true_marker_coords = [marker_grid_cols marker_grid_rows];
                                                            
                            % Your 5 extracted, unsorted image points [col, row]
                            markers = [marker_col, marker_row];

                            % Generate all 120 possible orderings of the numbers 1 through 5
                            all_orders = perms(1:5);
                            best_error = inf;
                            best_order = [];
                            
                            % Suppress warnings temporarily for the invalid permutations
                            warning('off', 'images:geotrans:matrixBadlyScaled');

                            
                            
                            % Loop through every possible arrangement
                            for i = 1:size(all_orders, 1)
                                test_order = all_orders(i, :);
                                test_markers = markers(test_order, :);
                                
                                % Assume the first 4 points map to the 4 square corners
                                moving_pts = test_markers(1:4, :);
                                fixed_pts  = true_marker_coords(1:4, :);
                                
                                try
                                    % Calculate the projective warp
                                    tform = fitgeotform2d(moving_pts, fixed_pts, 'affine');
                                    
                                    % Push the assumed 5th point (The Key) through the warp
                                    [proj_key_x, proj_key_y] = transformPointsForward(tform, ...
                                                                test_markers(5, 1), test_markers(5, 2));
                                    
                                    % Calculate how far the projected point is from [1, 0]
                                    err = norm([proj_key_x, proj_key_y] - true_marker_coords(5, :));
                                    
                                    % If this is the lowest error so far, save it
                                    if err < best_error
                                        best_error = err;
                                        best_order = test_order;
                                    end
                                catch
                                    % Invalid geometry (e.g., self-intersecting), skip it
                                    continue; 
                                end
                            end
                            warning('on', 'images:geotrans:matrixBadlyScaled');
                            
                            
                            % --- Apply the Winning Order ---
                            % best_order contains the exact indices mapping to true_marker_coords
                            sorted_markers = markers(best_order, :);

                            marker_col_ordered = sorted_markers(:, 1);
                            marker_row_ordered = sorted_markers(:, 2);

                            % close
                            % imagesc2(im); hold on
                            % plot(calib_col,calib_row,'.g');
                            % plot(marker_col_ordered,marker_row_ordered,'.r');
                            % plot(marker_col_ordered(1),marker_row_ordered(1),'.r','markersize',16);
                            % plot(marker_col_ordered(1:4),marker_row_ordered(1:4),'-r','linewidth',2);
                            % plot(marker_col_ordered(5),marker_row_ordered(5),'.m','markersize',16);
                            % drawnow
                            % return

                            %% New code (from Gemini) to map points using homography
                            % 1. Define the true grid indices of your 4 known markers.
                            % Replace these arrays with the actual (column, row) logical indices 
                            % in the 13x10 grid. For example, if they sit near the middle:
                            % marker_grid_col = [7 10 7 10 8]'; % The 1-13 column positions of your L-shape
                            % marker_grid_row = [5 5 8 8 5]'; % The 1-10 row positions of your L-shape
                            
                            % return
                            % 2. Calculate the projective transformation.
                            % This anchors your distorted image points to an ideal flat grid.
                            movingPoints = [marker_col_ordered, marker_row_ordered]; % Your known image coordinates
                            fixedPoints  = [marker_grid_cols, marker_grid_rows]; % Their mathematical positions
                            tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
                           
                            % 3. Project ALL detected points into the ideal grid space.
                            % [marker_ideal_col, marker_ideal_row] = transformPointsForward(tform, marker_col_ordered, marker_row_ordered);

                            [ideal_cols, ideal_rows] = transformPointsForward(tform, calib_col, calib_row);

                            
                            % plot(marker_grid_col, marker_grid_row,'.'); hold on
                            % plot(marker_ideal_col, marker_ideal_row,'.'); hold on
                            % 
                            % plot(ideal_cols, ideal_rows,'.');axis image
                            % return

                            
                            % 4. Round to the nearest integer to snap them to exact grid indices.
                            grid_cols_idx = round(ideal_cols);
                            grid_rows_idx = round(ideal_rows);

                            % This measures how far each point had to be "dragged" to snap to the grid
                            grid_deviation = sqrt((ideal_cols - grid_cols_idx).^2 + (ideal_rows - grid_rows_idx).^2);

                            % Define a tolerance limit (in grid units). 
                            tolerance = 0.30;
                            % return

                            grid_image_col = NaN(num_rows, num_cols);
                            grid_image_row = NaN(num_rows, num_cols);
                            best_deviation = inf(num_rows, num_cols); % Tracks the best fit for each slot

                            % 6. Populate the grid with the original image coordinates.
                            for i = 1:length(calib_col)
                                c = grid_cols_idx(i);
                                r = grid_rows_idx(i);
                                dev = grid_deviation(i);

                                % Verify the projected point is inside bounds AND within the snap tolerance
                                if (c >= 1 && c <= num_cols) && (r >= 1 && r <= num_rows) && (dev <= tolerance)
                                    
                                    % Collision Resolution: Only assign this point if the slot is empty, 
                                    % OR if this point is a closer match to the exact integer than the one currently there.
                                    if dev < best_deviation(r, c)
                                        % Place the original sub-pixel image coordinate into the right slot
                                        grid_image_col(r, c) = calib_col(i);
                                        grid_image_row(r, c) = calib_row(i);
                                        
                                        % Update the best deviation so future points are compared against this one
                                        best_deviation(r, c) = dev;
                                    end
                                end
                            end

                            % --- NEW STEP 7: CONNECTED COMPONENT FILTERING (ISLAND REMOVAL) ---
                            % We want to remove any isolated "islands" of points that are physically 
                            % separated from the main calibration grid.
                            
                            % 1. Create a binary mask of the current grid (1 = point, 0 = empty)
                            valid_mask = ~isnan(grid_image_col);
                            
                            % 2. Instantly erase any contiguous island with fewer than 5 points.
                            % The '8' tells it to count diagonal connections as part of the same island.
                            % This completely wipes out 1x1, 1x2, 2x1, and 2x2 clusters in one shot.
                            min_cluster_size = 5; 
                            cleaned_mask = bwareaopen(valid_mask, min_cluster_size, 8);
                            
                            % 3. Find exactly which points the filter decided to erase
                            isolated_mask = valid_mask & ~cleaned_mask;
                            
                            % 4. Erase them from the actual coordinate grid
                            grid_image_col(isolated_mask) = NaN;
                            grid_image_row(isolated_mask) = NaN;

                            % --- NEW STEP 7: ITERATIVE STRUCTURAL OUTLIER REMOVAL ---
                            % % Iteratively removes points with < 2 neighbors. 
                            % % Capped at 5 iterations to prevent infinite loops.
                            % 
                            % max_iterations = 5; 
                            % neighbor_kernel = [1 1 1; 
                            %                    1 0 1; 
                            %                    1 1 1];
                            % 
                            % for iter = 1:max_iterations
                            %     % 1. Assess the current state of the grid
                            %     valid_mask = ~isnan(grid_image_col);
                            % 
                            %     % 2. Count neighbors for every point
                            %     neighbor_counts = conv2(double(valid_mask), neighbor_kernel, 'same');
                            % 
                            %     % 3. Identify targets (Exists AND has < 2 neighbors)
                            %     isolated_mask = valid_mask & (neighbor_counts < 2);
                            % 
                            %     % 4. Equilibrium Check: Did we find any points to remove this time?
                            %     % If not, the grid is perfectly solid. Break out of the loop early!
                            %     if ~any(isolated_mask(:))
                            %         break; 
                            %     end
                            % 
                            %     % 5. Erase the isolated points and loop again
                            %     grid_image_col(isolated_mask) = NaN;
                            %     grid_image_row(isolated_mask) = NaN;
                            % end

                            %% proceed with rest of grid assignment
                            marker_grid_ind=sub2ind(size(grid_image_col),marker_grid_rows,marker_grid_cols);

                            grid_image_col(marker_grid_ind)=marker_col_ordered;
                            grid_image_row(marker_grid_ind)=marker_row_ordered;

                            final_grid_x=grid_image_col;
                            final_grid_y=grid_image_row;
                            
                            x_cam{f}=[final_grid_x(:) final_grid_y(:)]';
                            x_cam_points(f)=sum(~isnan(final_grid_x(:)));

                        end
                    end
                end
            end
        end
    end

end
