function [T_cams_2_ref, om_cams_2_ref,om_grids_2_ref,T_grids_2_ref, valid_frames, opt_ref_cam]=est_extrinsics_new2(om,T)
% save temp_est_extrinsics
% Entirely new version from Gemini
% Now removed ref_cam and instead select the camera with the most number of
% shared views

% clear
% close all
% load temp_est_extrinsics
% ref_cam=1;

num_cameras=length(om);
num_images = size(om{1}, 2);

% --- 1. FILTER VALID FRAMES ---
seen_matrix = false(num_cameras, num_images);
for c = 1:num_cameras
    seen_matrix(c, :) = ~isnan(T{c}(1, :));
end

cameras_per_frame = sum(seen_matrix, 1);
valid_frames = (cameras_per_frame >= 2);

% fprintf('Found %d valid frames (out of %d) where >= 2 cameras saw the target.\n', ...
%     sum(valid_frames), num_images);


% --- 2. PRECOMPUTE PAIRWISE CONNECTIONS (Using only valid frames) ---
overlaps = zeros(num_cameras, num_cameras);
rel_R = cell(num_cameras, num_cameras);
rel_T = cell(num_cameras, num_cameras);

% Extract just the indices of the valid frames to speed up the loop
valid_indices = find(valid_frames);

for cam_A = 1:num_cameras
    for cam_B = 1:num_cameras
        if cam_A == cam_B, continue; end
        
        om_candidates = [];
        T_candidates = [];
        
        for i = valid_indices
            % Check if BOTH cameras saw the grid in this specific valid frame
            if seen_matrix(cam_A, i) && seen_matrix(cam_B, i)
                R_A = rodrigues(om{cam_A}(:, i));
                R_B = rodrigues(om{cam_B}(:, i));
                T_A = T{cam_A}(:, i);
                T_B = T{cam_B}(:, i);
                
                % Transform from A to B
                R_B_A = R_B * R_A';
                T_B_A = T_B - R_B_A * T_A;
                
                om_candidates = [om_candidates, rodrigues(R_B_A)];
                T_candidates  = [T_candidates, T_B_A];
            end
        end
        
        overlaps(cam_A, cam_B) = size(T_candidates, 2);
        
        if overlaps(cam_A, cam_B) > 0
            rel_R{cam_A, cam_B} = rodrigues(median(om_candidates, 2));
            rel_T{cam_A, cam_B} = median(T_candidates, 2);
        end
    end
end

% --- 2.5 DYNAMICALLY SELECT BEST REFERENCE CAMERA ---
% Sum the rows of the overlaps matrix to find the total number of 
% shared views each camera has with the entire network.
total_connections = sum(overlaps, 2);

% Find the camera with the absolute maximum number of shared views
[max_connections, best_ref_cam] = max(total_connections);
opt_ref_cam = best_ref_cam;

% fprintf('\nAutomatically selected Cam %d as World Origin (Strongest network connections: %d shared views)\n', ...
%     ref_cam, max_connections);

% --- 3. BUILD THE MAXIMUM SPANNING TREE ---
% Start the tree at the designated reference camera
solved_cams = opt_ref_cam;
unsolved_cams = setdiff(1:num_cameras, opt_ref_cam);

global_R = cell(1, num_cameras);
global_T = cell(1, num_cameras);

% The reference camera is the origin (Rotation = Identity, Translation = 0)
global_R{opt_ref_cam} = eye(3);
global_T{opt_ref_cam} = zeros(3, 1);

% fprintf('\n--- Building Camera Network (Origin: Cam %d) ---\n', opt_ref_cam);

while ~isempty(unsolved_cams)
    best_overlap = -1;
    best_A = -1;
    best_B = -1;
    
    for s = solved_cams
        for u = unsolved_cams
            if overlaps(s, u) > best_overlap
                best_overlap = overlaps(s, u);
                best_A = s;
                best_B = u;
            end
        end
    end
    
    if best_overlap == 0
        warning('Cameras [%s] share no overlapping views with the network!', num2str(unsolved_cams));
        break; 
    end
    
    % fprintf('Linking Cam %d to Cam %d (Based on %d shared views)\n', best_B, best_A, best_overlap);
    
    % Daisy-chain the transformations mathematically
    R_B_A = rel_R{best_A, best_B};
    T_B_A = rel_T{best_A, best_B};
    
    global_R{best_B} = R_B_A * global_R{best_A};
    global_T{best_B} = R_B_A * global_T{best_A} + T_B_A;
    
    % Update lists
    solved_cams = [solved_cams, best_B];
    unsolved_cams(unsolved_cams == best_B) = [];
end


% --- 4. FORMAT FINAL OUTPUT ---
om_extrinsic = zeros(3, num_cameras);
T_extrinsic  = zeros(3, num_cameras);

for c = 1:num_cameras
    if ismember(c, solved_cams)
        om_extrinsic(:, c) = rodrigues(global_R{c});
        T_extrinsic(:, c)  = global_T{c};
    else
        om_extrinsic(:, c) = NaN;
        T_extrinsic(:, c)  = NaN;
    end
end

om_grids_2_ref=om{opt_ref_cam};
T_grids_2_ref=T{opt_ref_cam};

% --- 5. INTERPOLATE MISSING REFERENCE CAMERA POSES ---
% Find frames where the reference camera missed the grid, 
% but at least 2 other cameras saw it (valid_frames).
missing_ref_frames = find(~seen_matrix(opt_ref_cam, :) & valid_frames);

% fprintf('\n--- Estimating Missing Grid Poses for Ref Cam (Cam %d) ---\n', opt_ref_cam);
% fprintf('Found %d frames to estimate.\n', length(missing_ref_frames));

for i = missing_ref_frames
    om_estimates = [];
    T_estimates = [];
    
    for c = 1:num_cameras
        % Only use cameras that actually saw the grid in frame i,
        % AND successfully got a global pose in Step 3.
        if c ~= opt_ref_cam && seen_matrix(c, i) && ismember(c, solved_cams)
            
            % 1. Global Extrinsics of Camera 'c' (from Ref Cam to c)
            R_ext = rodrigues(om_extrinsic(:, c));
            T_ext = T_extrinsic(:, c);
            
            % 2. Local Pose of the grid relative to Camera 'c'
            R_grid_c = rodrigues(om{c}(:, i));
            T_grid_c = T{c}(:, i);
            
            % 3. Back-calculate the grid's pose relative to the Ref Cam
            R_grid_ref = R_ext' * R_grid_c;
            T_grid_ref = R_ext' * (T_grid_c - T_ext);
            
            % Store the candidates
            om_estimates = [om_estimates, rodrigues(R_grid_ref)];
            T_estimates  = [T_estimates, T_grid_ref];
        end
    end
    
    % If we generated valid estimates from the network, take the median
    if ~isempty(T_estimates)
        % Overwrite the NaNs with our new estimated seed values
        om_grids_2_ref(:, i) = median(om_estimates, 2);
        T_grids_2_ref(:, i)  = median(T_estimates, 2);
        
        % Update the logical matrix so the rest of your pipeline knows 
        % the ref_cam now has valid data for this frame
        seen_matrix(opt_ref_cam, i) = true;
    end
end

% ---  6. Initialize arrays for the inverted (Camera-to-World) poses ---
om_cams_2_ref = zeros(3, num_cameras);
T_cams_2_ref  = zeros(3, num_cameras);

for c = 1:num_cameras
    if ~any(isnan(om_extrinsic(:, c)))
        
        % 1. Invert the rotation
        % A simple negation of the Rodrigues vector gives the inverse rotation
        om_cams_2_ref(:, c) = -om_extrinsic(:, c);
        
        % 2. Invert the translation
        % T_c2w = -R^T * T_w2c
        R_inv = rodrigues(om_cams_2_ref(:, c)); % This is R^T
        T_cams_2_ref(:, c) = -R_inv * T_extrinsic(:, c);
        
    else
        om_cams_2_ref(:, c) = NaN;
        T_cams_2_ref(:, c)  = NaN;
    end
end

