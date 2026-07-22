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
   
            
function F=function_multi_cam_ext(xx,chi,XX,fc,cc,kc,num_cam)

%% Optimised using Gemini to make it a bit more efficient
om_cams_2_ref = chi(:, 1:num_cam);
T_cams_2_ref  = chi(:, num_cam+1:2*num_cam);
om_grids_2_ref= chi(:, 2*num_cam+1:2*num_cam+size(XX,2));
T_grids_2_ref = chi(:, 2*num_cam+size(XX,2)+1:end);

% 1. Calculate 3xn matrix X
X = NaN(3, size(xx,2));
ii = 1;
for i = 1:size(XX,2)
    X_temp = XX{i};
    cols = size(X_temp,2);
    
    % OPTIMIZATION: Removed 'ones()' allocation. 
    % Implicit expansion handles vector addition natively.
    X(:, ii:ii+cols-1) = rodrigues(om_grids_2_ref(:,i)) * X_temp + T_grids_2_ref(:,i);
    ii = ii + cols;
end

x_all = xx(1:2:end, :);
y_all = xx(2:2:end, :);
n = size(xx, 2);
F = NaN(size(xx));

for i = 1:num_cam
    
    cx = fc(1,i);
    cy = fc(2,i);
    
    R = rodrigues(om_cams_2_ref(:,i))';
    T_vec = T_cams_2_ref(:,i); % Just the 3x1 vector, no repmat
    
    % Only using the 2 radial terms as requested
    k1 = kc(1,i);
    k2 = kc(2,i);
    
    % 1. Calculate 3D Pinhole Projections
    dX = T_vec - X; 
    proj_pts = R * dX; 
    
    % 2. Normalize to Z=1 plane
    xn = proj_pts(1,:) ./ proj_pts(3,:);
    yn = proj_pts(2,:) ./ proj_pts(3,:);
    
    % 3. Calculate Radius from the mathematically perfect 3D projections!
    r_n2 = xn.^2 + yn.^2; 
    
    % 4. Apply FORWARD Radial Distortion (Plus signs, not minus!)
    dist_factor = ones(1,n) + k1.*r_n2 + k2.*(r_n2.^2);
    
    xd = xn .* dist_factor;
    yd = yn .* dist_factor;
    
    % 5. Scale to final Distorted Pixels
    proj_pixel_x = cx * xd + cc(1,i);
    proj_pixel_y = cy * yd + cc(2,i);
        
    % 6. Calculate True Pixel Residual (Model Projection - Raw Observation)
    F(2*i-1:2*i, :) = [
        proj_pixel_x - x_all(i,:);
        proj_pixel_y - y_all(i,:)
    ];

    % % Older version with inverse distortion model
    % % K = kc(1:3,i);
    % % P = kc(4:5,i);
    % % 
    % % dx = x_all(i,:) - cc(1,i);
    % % dy = y_all(i,:) - cc(2,i);
    % % 
    % % dx_n = dx / cx;
    % % dy_n = dy / cy;
    % % r_n2 = dx_n.^2 + dy_n.^2; % Avoided doing sqrt() early to save ops on r_n.^2
    % % 
    % % % OPTIMIZATION: Calculate (T - X) exactly once per camera using implicit expansion
    % % dX = T_vec - X; 
    % % 
    % % % OPTIMIZATION: Matrix multiply all 3 rows simultaneously
    % % proj_pts = R * dX; 
    % % 
    % % % Extract for division
    % % proj_X = proj_pts(1,:);
    % % proj_Y = proj_pts(2,:);
    % % proj_Z = proj_pts(3,:);
    % % 
    % % % Pre-calculate distortion blocks to keep the final assignment clean
    % % dist_rad = ones(1,n) - (K(1).*r_n2 + K(2).*(r_n2.^2) + K(3).*(r_n2.^3));
    % % dist_tang_x = 2*P(1).*dx_n.*dy_n + P(2).*(r_n2 + 2*dx_n.^2);
    % % dist_tang_y = 2*P(2).*dx_n.*dy_n + P(1).*(r_n2 + 2*dy_n.^2);
    % % 
    % % F(2*i-1:2*i, :) = [
    % %     cx * (proj_X ./ proj_Z) - (dx .* dist_rad + dist_tang_x);
    % %     cy * (proj_Y ./ proj_Z) - (dy .* dist_rad + dist_tang_y)
    % % ];
end

Fnan = isnan(F);
F(Fnan) = 0;