%Function function_multi_cam runs collinearity equations for 4-camera setup

%Function four_cams2 takes as its arguments:
   
    % - x (an 8-by-n matrix containing pixel coordinates of the n points in the four cameras in the order
        %[x1,y1,x2,y2,x3,y3,x4,y4])
    % - the 3-by-18+n matrix chi, concatenating:
        % - Rodrigues vectors for cameras 2-4 (:,1:3)
        % - Pose vectors for cameras 2-4 (:,4-6)
        % - A 3-by-3 matrix for each of cameras 1-4 containing:
            % - the focal length in horizontal (1,1) and vertical (2,1)
            % pixels
            % - the principal point coordinates (3,1) to (1,2)
            % - the tangential distortion coefficients P1 and P2 (2:3,2)
            % - the radial distortion coefficents K1, K2, K3 (:,3)    
    % - X (a 3-by-n matrix containing the X, Y, Z coordinates of the n
        %points in the axis system of camera 1)
   
            
function F=function_multi_cam_full(xx,chi,XX,num_cam)

num_grids = size(XX, 2);

    % --- 1. SPLIT EXTRINSICS AND INTRINSICS ---
    om_cams_2_ref  = chi(:, 1:num_cam);
    T_cams_2_ref   = chi(:, num_cam+1:2*num_cam);
    om_grids_2_ref = chi(:, 2*num_cam+1:2*num_cam+num_grids);
    T_grids_2_ref  = chi(:, 2*num_cam+num_grids+1:end);
    
    % The intrinsic parameters matrix starts exactly after the grid extrinsics
    chi_int = chi(:, 2*num_cam+2*num_grids+1 : end);

    % --- 2. CALCULATE 3D POINTS (X) ---
    X = NaN(3, size(xx,2));
    ii = 1;
    for i = 1:num_grids
        X_temp = XX{i};
        cols = size(X_temp,2);
        
        % OPTIMIZATION: Implicit expansion handles vector addition natively.
        X(:, ii:ii+cols-1) = rodrigues(om_grids_2_ref(:,i)) * X_temp + T_grids_2_ref(:,i);
        ii = ii + cols;
    end

    x_all = xx(1:2:end, :);
    y_all = xx(2:2:end, :);
    n = size(xx, 2);
    F = NaN(size(xx));

    % --- 3. CALCULATE RESIDUALS (FORWARD MODEL) ---
    for i = 1:num_cam
        
        % Unpack parameters (Renamed variables to prevent fc/cc confusion)
        col = 2 * i - 1;
        fc_val = chi_int(1, col); % Focal length (Enforcing square pixels fx = fy)
        cc_x   = chi_int(2, col); % Principal Point X
        cc_y   = chi_int(3, col); % Principal Point Y
        
        K = [chi_int(1:2, col+1); 0]; % Radial [k1; k2; k3]
        P = [0; 0];                   % Tangential [p1; p2]
        
        R = rodrigues(om_cams_2_ref(:,i))';
        T_vec = T_cams_2_ref(:,i); 
        
        % 3D Projection
        dX = T_vec - X; 
        proj_pts = R * dX; 
        
        % Normalize to Z=1 plane
        xn = proj_pts(1,:) ./ proj_pts(3,:);
        yn = proj_pts(2,:) ./ proj_pts(3,:);
        
        % Calculate Radius from the mathematically perfect 3D projections
        r2 = xn.^2 + yn.^2; 
        
        % Forward Radial Distortion Multiplier (Plus signs!)
        dist_factor = ones(1,n) + K(1).*r2 + K(2).*(r2.^2) + K(3).*(r2.^3);
        
        % Forward Tangential Distortion Math
        tang_x = 2*P(1).*xn.*yn + P(2).*(r2 + 2*xn.^2);
        tang_y = P(1).*(r2 + 2*yn.^2) + 2*P(2).*xn.*yn;
        
        % Apply Full Forward Distortion
        xd = xn .* dist_factor + tang_x;
        yd = yn .* dist_factor + tang_y;
        
        % Scale to final Distorted Pixels
        proj_pixel_x = fc_val * xd + cc_x;
        proj_pixel_y = fc_val * yd + cc_y;
            
        % Calculate True Pixel Residual (Model Projection - Raw Observation)
        F(2*i-1:2*i, :) = [
            proj_pixel_x - x_all(i,:);
            proj_pixel_y - y_all(i,:)
        ];
    end
    
    Fnan = isnan(F);
    F(Fnan) = 0;
    
end

% % % Older version with incorrect, inverse distortion model
% num_grids = size(XX, 2);
% 
% %% Optimised using Gemini to make it a bit more efficient
% om_cams_2_ref = chi(:, 1:num_cam);
% T_cams_2_ref  = chi(:, num_cam+1:2*num_cam);
% om_grids_2_ref= chi(:, 2*num_cam+1:2*num_cam+size(XX,2));
% T_grids_2_ref = chi(:, 2*num_cam+size(XX,2)+1:end);
% 
% % The intrinsic parameters matrix starts exactly after the grid extrinsics
% chi_int = chi(:, 2*num_cam+2*num_grids+1 : end);
% 
% % 1. Calculate 3xn matrix X
% X = NaN(3, size(xx,2));
% ii = 1;
% for i = 1:size(XX,2)
%     X_temp = XX{i};
%     cols = size(X_temp,2);
% 
%     % OPTIMIZATION: Removed 'ones()' allocation. 
%     % Implicit expansion handles vector addition natively.
%     X(:, ii:ii+cols-1) = rodrigues(om_grids_2_ref(:,i)) * X_temp + T_grids_2_ref(:,i);
%     ii = ii + cols;
% end
% 
% x_all = xx(1:2:end, :);
% y_all = xx(2:2:end, :);
% n = size(xx, 2);
% F = NaN(size(xx));
% 
% for i = 1:num_cam
% 
%     col = 2 * i - 1;
% 
%     cx = chi_int(1, col);
%     cy = chi_int(1, col);
% 
%     cc=chi_int(2:3, col);
% 
%     K=[chi_int(1:2, col+1);0];
%     P=[0 0];
% 
%     R = rodrigues(om_cams_2_ref(:,i))';
%     T_vec = T_cams_2_ref(:,i); % Just the 3x1 vector, no repmat
% 
%     dx = x_all(i,:) - cc(1);
%     dy = y_all(i,:) - cc(2);
% 
%     dx_n = dx / cx;
%     dy_n = dy / cy;
%     r_n2 = dx_n.^2 + dy_n.^2; % Avoided doing sqrt() early to save ops on r_n.^2
% 
%     % OPTIMIZATION: Calculate (T - X) exactly once per camera using implicit expansion
%     dX = T_vec - X; 
% 
%     % OPTIMIZATION: Matrix multiply all 3 rows simultaneously
%     proj_pts = R * dX; 
% 
%     % Extract for division
%     proj_X = proj_pts(1,:);
%     proj_Y = proj_pts(2,:);
%     proj_Z = proj_pts(3,:);
% 
%     % Pre-calculate distortion blocks to keep the final assignment clean
%     dist_rad = ones(1,n) - (K(1).*r_n2 + K(2).*(r_n2.^2) + K(3).*(r_n2.^3));
%     dist_tang_x = 2*P(1).*dx_n.*dy_n + P(2).*(r_n2 + 2*dx_n.^2);
%     dist_tang_y = 2*P(2).*dx_n.*dy_n + P(1).*(r_n2 + 2*dy_n.^2);
% 
%     F(2*i-1:2*i, :) = [
%         cx * (proj_X ./ proj_Z) - (dx .* dist_rad + dist_tang_x);
%         cy * (proj_Y ./ proj_Z) - (dy .* dist_rad + dist_tang_y)
%     ];
% end
% 
% Fnan = isnan(F);
% F(Fnan) = 0;