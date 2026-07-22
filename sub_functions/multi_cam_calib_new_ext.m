function [om_cams_opt,T_cams_opt,om_grids_opt,T_grids_opt,residuals_ext]=multi_cam_calib_new_ext(fc,cc,kc,om_cams_2_ref,T_cams_2_ref,om_grids_2_ref,T_grids_2_ref,x_cam,XX_cam,num_cam,user_ref_cam,valid_frames)
% 
% Note this version doesn't assume the camera initally as [0,0,0] is the
% reference camera so after the bundle adjustment it corrects everything to
% that camera. This can improve performance when the reference camera has
% few shared images.
% Only does extrinsics.

% % % 
% save temp_multi_cam_calib
% return
% clear
% close all
% 
% load temp_multi_cam_calib

om_grids_2_ref=om_grids_2_ref(:,valid_frames==1);
T_grids_2_ref=T_grids_2_ref(:,valid_frames==1);

x_cam=x_cam(:,valid_frames==1);
% XX_cam=XX_cam(:,valid_frames==1);

% fc=repmat(mean(fc),2,1);
% 
% return

% for each grid point, find the number of camera views its been found in
% and remove any with no camera views (a single camera view should now be
% ok).

% New code that uses est

num_grids=size(om_grids_2_ref,2);

for i=1:num_grids
   
    xx=cell2mat(x_cam(:,i));
    X=XX_cam;
    
    xnan=sum(~isnan(xx));
    
    xx(:,xnan<4)=[];
    X(:,xnan<4)=[];
    
    x_cell{i}=xx;
    X_cell{i}=X;
    
    grid_l(i)=size(xx,2);
    
end

% remove any grids where no points have been found
x_cell(:,grid_l==0)=[];
X_cell(grid_l==0)=[];

T_grids_2_ref(:,grid_l==0)=[];
om_grids_2_ref(:,grid_l==0)=[];
grid_l(grid_l==0)=[];

% % reduce numbers, just for testing
% use_grids=1:30:300;
% grid_l=grid_l(use_grids);
% 
% T_grids_2_ref=T_grids_2_ref(:,use_grids);
% om_grids_2_ref=om_grids_2_ref(:,use_grids);
% 
% x_cell=x_cell(use_grids);
% X_cell=X_cell(use_grids);
% 
% for i=1:length(x_cell)
%     x_cell{i}=x_cell{i}(:,1:5:40);
%     X_cell{i}=X_cell{i}(:,1:5:40);
%     grid_l(i)=size(x_cell{i},2);
% 
% end

% return

xx=cell2mat(x_cell);

chi0_ext=[om_cams_2_ref T_cams_2_ref om_grids_2_ref T_grids_2_ref];

% initial_residuals = function_multi_cam_new2(xx, chi0, X_cell, fc, cc, kc, num_cam);
% return
grid_t=numel(xx);

for i=1:num_cam
    col_start=3*i-2;
    
    jcol_cam_temp=repmat([col_start:col_start+2 3*(num_cam)+col_start:3*(num_cam)+col_start+2],sum(grid_l),1);
    jcol_cam_temp(isnan(xx(2*i,:)),:)=[];
    
    jcol_cam{i,1}=repmat(jcol_cam_temp,2,1);

    jrow_cam_temp=(2*i-1:num_cam*2:sum(grid_t))';
    jrow_cam_temp(isnan(xx(2*i,:)))=[];
    
    jrow_cam{i,1}=repmat([jrow_cam_temp;jrow_cam_temp+1],1,6);
    
end

jcol_cam=cell2mat(jcol_cam);
jrow_cam=cell2mat(jrow_cam);

jcol_grid=[];
jrow_grid=[];
jrow_max=0;

for i=1:length(grid_l)
    
    jcol_grid=[jcol_grid;repmat([6*(num_cam)+3*i-2:6*(num_cam)+3*i 6*(num_cam)+3*i-2+3*length(grid_l):6*(num_cam)+3*i+3*length(grid_l)],num_cam*2*grid_l(i),1)];
    jrow_grid=[jrow_grid;repmat((1:num_cam*2*grid_l(i))',1,6)+jrow_max];
    jrow_max=max(jrow_grid(:));
end

a=ismember(jrow_grid,jrow_cam);
jcol_grid(a==0)=[];
jrow_grid(a==0)=[];

jcol=[jcol_cam(:);jcol_grid(:)];
jrow=[jrow_cam(:);jrow_grid(:)];

JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,6*(num_cam)+3*length(grid_l)+3*length(grid_l));

disp('Optimising extrinsics')

options=optimset('lsqnonlin');
options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',50,'MaxFunEvals',300*size(xx,2),'TolFun',1e-5,'TolX', 1e-5,'PrecondBandWidth', 2);
[chi_ext,resnorm_ext,residuals_ext]=lsqnonlin(@(chi_ext) function_multi_cam_ext(xx,chi_ext,X_cell,fc,cc,kc,num_cam),chi0_ext,[],[],options);

%% Extract values from Chi
om_cams_2_ref_opt=chi_ext(:,1:num_cam);
T_cams_2_ref_opt=chi_ext(:,num_cam+1:2*num_cam);

om_grids_2_ref_opt=chi_ext(:,2*num_cam+1:2*num_cam+size(X_cell,2));
T_grids_2_ref_opt=chi_ext(:,2*num_cam+size(X_cell,2)+1:end);

%% Set reference to user_ref_cam
R_ref_cam = rodrigues(om_cams_2_ref_opt(:, user_ref_cam));
T_ref_cam = T_cams_2_ref_opt(:, user_ref_cam);
   
om_cams_opt=NaN(3,num_cam);
T_cams_opt=NaN(3,num_cam);

for c=1:num_cam

    R_cams_2_ref_opt=rodrigues(om_cams_2_ref_opt(:,c));

    R_cams_opt=R_ref_cam'*R_cams_2_ref_opt; % seems correct

    om_cams_opt(:,c)=rodrigues(R_cams_opt);
    T_cams_opt(:,c)=R_ref_cam'*(T_cams_2_ref_opt(:,c)-T_ref_cam); % think is correct

end

%% Now do grids
clear om_grids T_grids
om_grids_opt=NaN(3,size(om_grids_2_ref_opt,2));
T_grids_opt=NaN(3,size(om_grids_2_ref_opt,2));

for i=1:size(om_grids_2_ref_opt,2)
    
    if ~isnan(om_grids_2_ref_opt(1,i))
        R_grids_2_ref_opt=rodrigues(om_grids_2_ref_opt(:,i));
    
        R_grids_opt=R_ref_cam'*R_grids_2_ref_opt; % seems correct
    
        om_grids_opt(:,i)=rodrigues(R_grids_opt);
        T_grids_opt(:,i)=R_ref_cam'*(T_grids_2_ref_opt(:,i)-T_ref_cam); % think is correct
    end
end

% %% --- SHIFT ENTIRE NETWORK TO GRID 1 AXIS SYSTEM ---
% ref_grid = 1;
% 
% % 1. Extract the pose of Grid 1 (currently in the Ref Camera's coordinate system)
% R_ref_grid = rodrigues(om_grids_opt(:, ref_grid));
% T_ref_grid = T_grids_opt(:, ref_grid);
% 
% % 2. Preallocate the final grid-aligned arrays
% om_cams_2_ref_grid  = NaN(3, num_cam);
% T_cams_2_ref_grid   = NaN(3, num_cam);
% om_grids_2_ref_grid = NaN(3, size(om_grids_opt, 2));
% T_grids_2_ref_grid  = NaN(3, size(om_grids_opt, 2));
% 
% % 3. Transform all Cameras into the Grid 1 frame
% for c = 1:num_cam
% 
%     R_cams_opt=rodrigues(om_cams_opt(:,c));
%     R_cams_2_ref_grid=R_ref_grid'*R_cams_opt; % seems correct
% 
%     om_cams_2_ref_grid(:,c)=rodrigues(R_cams_2_ref_grid);
%     T_cams_2_ref_grid(:,c)=R_ref_grid'*(T_cams_opt(:,c)-T_ref_grid); % think is correct
% 
% end
% 
% for i = 1:size(om_grids_opt, 2)
% 
%     if ~isnan(om_grids_opt(1,i))
% 
%         R_grids_opt=rodrigues(om_grids_opt(:,i));
%         R_cams_2_ref_grid=R_ref_grid'*R_grids_opt; % seems correct
% 
%         om_grids_2_ref_grid(:,i)=rodrigues(R_cams_2_ref_grid);
%         T_grids_2_ref_grid(:,i)=R_ref_grid'*(T_grids_opt(:,i)-T_ref_grid); % think is correct
% 
%     end
% 
% end
% 
% % P
% working_distance=mean(sqrt(sum(T_cams_2_ref_opt.^2)));
% nx=18;
% 
% Cx=[-1 1 1 -1 -1; -0.5 0.5 0.5 -0.5 -0.5]*0.1*working_distance; % makes camera shape for drawing
% Cy=[-1 -1 1 1 -1; -0.5 -0.5 0.5 0.5 -0.5]*0.1*working_distance;
% Cz=[0 0 0 0 0; -1 -1 -1 -1 -1]*0.1*working_distance;
% 
% C=[Cx(:)';Cy(:)';Cz(:)'];
% 
% figure
% 
% for c=1:num_cam
% 
%     % initial estimates
%     C2=rodrigues(om_cams_2_ref(:,c))*C+repmat(T_cams_2_ref(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=[0.8 0 0];
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     % optimised estimates
%     C2=rodrigues(om_cams_2_ref_opt(:,c))*C+repmat(T_cams_2_ref_opt(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     % col=rand(1,3);
% 
%     col=[0 0.8 0];
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('Z','X','Y');
%     set(gca,'ydir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% 
% % for i=1:size(om_grids_2_ref_opt,2)
% % 
% %     R_grids_opt=rodrigues(om_grids_2_ref_opt(:,i));
% % 
% %     XX2=R_grids_opt*XX_cam{i}+repmat(T_grids_2_ref_opt(:,i),1,size(XX_cam{i},2));
% % 
% %     X2=reshape(XX2(1,:),nx,[]);
% %     Y2=reshape(XX2(2,:),nx,[]);
% %     Z2=reshape(XX2(3,:),nx,[]);
% % 
% %     col=rand(1,3);
% % 
% %     p1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o','color',col); hold on
% %     p2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.','color',col); hold on
% % %             
% % end
% 
% title('Aligned to initial optimal reference camera')
% 
% return
% 
% figure
% 
% for c=1:num_cam
%     C2=rodrigues(om_cams_opt(:,c))*C+repmat(T_cams_opt(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=rand(1,3);
% 
%     h=plot3(Cz2,Cx2,Cy2,Cz2',Cx2',Cy2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cz2(:)),mean(Cx2(:)),mean(Cy2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cz2([1 2 4 3]),Cx2([1 2 4 3]),Cy2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cz2(:,1:4),2),mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('Z','X','Y');
%     set(gca,'ydir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% 
% for i=1:size(om_grids_2_ref_opt,2)
% 
%     R_grids_opt=rodrigues(om_grids_opt(:,i));
% 
%     XX2=R_grids_opt*XX_cam{i}+repmat(T_grids_opt(:,i),1,size(XX_cam{i},2));
% 
%     X2=reshape(XX2(1,:),nx,[]);
%     Y2=reshape(XX2(2,:),nx,[]);
%     Z2=reshape(XX2(3,:),nx,[]);
% 
%     col=rand(1,3);
% 
%     p1=plot3(Z2([1:3 nx+1]),X2([1:3 nx+1]),Y2([1:3 nx+1]),'o','color',col); hold on
%     p2=plot3(Z2([4:nx nx+2:end]),X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),'.','color',col); hold on
% %             
% end
% title('Aligned to reference camera')
% 
% figure
% for c=1:num_cam
%     C2=rodrigues(om_cams_2_ref_grid(:,c))*C+repmat(T_cams_2_ref_grid(:,c),1,size(C,2));
% 
%     Cx2=reshape(C2(1,:),2,[]);
%     Cy2=reshape(C2(2,:),2,[]);
%     Cz2=reshape(C2(3,:),2,[]);
% 
%     col=rand(1,3);
% 
%     h=plot3(Cx2,Cy2,Cz2,Cx2',Cy2',Cz2'); hold on; axis image
%     set(h,'color',col,'linewidth',2);
%     h2=text(mean(Cx2(:)),mean(Cy2(:)),mean(Cz2(:))+0.15*working_distance,num2str(c));
%     set(h2,'color',col)
%     patch(Cx2([1 2 4 3]),Cy2([1 2 4 3]),Cz2([1 2 4 3]),col); axis image
% 
%     h=plot3(mean(Cx2(:,1:4),2),mean(Cy2(:,1:4),2),mean(Cz2(:,1:4),2)); hold on; axis image
%     set(h,'color',col,'linewidth',2);
% 
%     labels('X','Y','Z');
%     set(gca,'xdir','reverse','zdir','reverse')
% 
%     axis image
%     labels
% end
% 
% for i=1%:size(om_grids_2_ref_opt,2)
% 
%     R_grids_opt=rodrigues(om_grids_2_ref_grid(:,i));
% 
%     XX2=R_grids_opt*XX_cam{i}+repmat(T_grids_2_ref_grid(:,i),1,size(XX_cam{i},2));
% 
%     X2=reshape(XX2(1,:),nx,[]);
%     Y2=reshape(XX2(2,:),nx,[]);
%     Z2=reshape(XX2(3,:),nx,[]);
% 
%     col=rand(1,3);
% 
%     p1=plot3(X2([1:3 nx+1]),Y2([1:3 nx+1]),Z2([1:3 nx+1]),'o','color',col); hold on
%     p2=plot3(X2([4:nx nx+2:end]),Y2([4:nx nx+2:end]),Z2([4:nx nx+2:end]),'.','color',col); hold on
% %             
% end
% title('Aligned to first grid')
