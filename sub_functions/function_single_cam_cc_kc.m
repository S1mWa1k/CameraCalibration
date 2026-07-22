%Function four_cams2 runs collinearity equations for 4-camera setup

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
   
            
function F=function_single_cam_cc_kc(xx,chi,XX)

num_grids = size(XX, 2);

% --- 1. UNPACK PARAMETERS ---
cx = chi(1,1);
cy = chi(1,1);         % Enforcing square pixels
cc = chi(1:2,2);       % Principal point
K  = [chi(1:2,3); 0];  % Radial distortion [k1; k2; k3]
P  = [0; 0];           % Tangential distortion [p1; p2]

om_grids = chi(:, 4 : num_grids+3);
T_grids  = chi(:, num_grids+4 : end);

% --- 2. CALCULATE 3D POINTS IN CAMERA FRAME ---
X = NaN(3, size(xx,2));
ii = 1;
for i = 1:num_grids
    X_temp = XX{i};
    cols = size(X_temp, 2);
    
    % OPTIMIZATION: Implicit expansion natively handles translation vector
    X(:, ii:ii+cols-1) = rodrigues(om_grids(:,i)) * X_temp + T_grids(:,i);
    ii = ii + cols;
end

x_all = xx(1:2:end, :);
y_all = xx(2:2:end, :);
n = size(xx, 2);

% --- 3. FORWARD DISTORTION MODEL ---
% Since the camera is at the origin, the 3D points X are 
% already in the correct coordinate frame for projection.

% Pinhole projection onto Z=1 plane
xn = X(1,:) ./ X(3,:);
yn = X(2,:) ./ X(3,:);

% Radius from the perfect pinhole projection
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
proj_pixel_x = cx * xd + cc(1);
proj_pixel_y = cy * yd + cc(2);

% --- 4. CALCULATE RESIDUALS ---
% True Pixel Residual (Model Projection - Raw Observation)
F = [
    proj_pixel_x - x_all;
    proj_pixel_y - y_all
];

Fnan = isnan(F);
F(Fnan) = 0;

% % Old version with incorrect, inverse distoration model
% om_grids=chi(:,4:size(XX,2)+3);
% T_grids=chi(:,size(XX,2)+4:end);
% 
% X=[];
% 
% for i=1:size(XX,2)
%     X_temp=XX{i};
%     X=[X rodrigues(om_grids(:,i))*X_temp+T_grids(:,i)*ones(1,size(X_temp,2))];
% end
% 
% x_all=xx(1:2:end,:);
% y_all=xx(2:2:end,:);
% 
% %Number of points
% n=size(xx,2);
% 
% cc=chi(1:2,2);
% 
% cx=chi(1,1);
% cy=chi(1,1);
% 
% dx=x_all-cc(1);
% dy=y_all-cc(2);
% 
% dx_n=dx/cx;
% dy_n=dy/cy;
% 
% r_n=sqrt(dx_n.^2+dy_n.^2);
% 
% R=rodrigues(zeros(3,1));
% T=zeros(3,n);
% 
% K=[chi(1:2,3);0];
% 
% P=[0 0];
% 
% F=[cx*(R(1,:)*(T-X))./(R(3,:)*(T-X))-(dx.*(ones(1,n)-(K(1).*r_n.^2+K(2).*r_n.^4+K(3).*r_n.^6)+2*P(1).*dx_n.*dy_n+P(2)*(r_n.^2+2*dx_n.^2)));
%    cy*(R(2,:)*(T-X))./(R(3,:)*(T-X))-(dy.*(ones(1,n)-(K(1).*r_n.^2+K(2).*r_n.^4+K(3).*r_n.^6)+2*P(2).*dx_n.*dy_n+P(1)*(r_n.^2+2*dy_n.^2)))];
% 
% Fnan=isnan(F);
% F(Fnan)=0;