function [fc_opt, cc_opt, kc_opt, T_opt, om_opt, x_cam, removed, residual]=single_cam_calib_cc_kc(fc,cc,kc,T,om,x_cam,XX_cam)

% save temp_single_cam

% clear 
% close all force
% load temp_single_cam

% clear disp
% disp='off';

% use=1;
% 
% x_cam=x_cam(use);
% XX_cam=XX_cam(use);
% 
% T=T(:,use);
% om=om(:,use);

% then with distortion
use_grid=zeros(size(x_cam));

k=1;

om_use=[];
T_use=[];

for i=1:length(x_cam)
    
    xx=x_cam{i};
%     xx=xx(:,1:10);
    
    if ~isempty(xx) && any(~isnan(xx(1,:)))

        [Ay, Ax]=find(isnan(xx));
        xx(:,Ax)=[];
        
        X=XX_cam;
%         X=X(:,1:10);
        X(:,Ax)=[];
        
        x_cell{k}=xx;
        X_cell{k}=X;
        
        grid_l(k)=size(X,2);

        om_use(:,k)=om(:,i);
        T_use(:,k)=T(:,i);

        use_grid(i)=1;
        
        k=k+1;
    end
    
end

use_grid_ind=find(use_grid);

changed=0;
removed=zeros(size(x_cam));

for runs=1:2

    if runs==1 || (runs==2 && changed==1)

        xx=cell2mat(x_cell);
        grid_t=numel(xx);
        
        %BUNDLE ADJUSTMENT ROUTINE
        %Specifies structure of the Jacobian matrix in the sparse matrix JStr
        
        clear j*
        chi0=[[fc;0] [cc;0] [kc(1:2);0] om_use T_use];

        % jcol_cam=repmat([1 2 4 5 7 8 1 2 4 5 7 8],size(xx,2),1);
        jcol_cam=repmat([1 4 5 7 8 1 4 5 7 8],size(xx,2),1);
        
        jrow_cam=[repmat(1:2:2*size(xx,2),1,0.5*size(jcol_cam,2)) repmat(2:2:2*size(xx,2),1,0.5*size(jcol_cam,2))]';
        jcol_cam=jcol_cam(:);
        
        % return
        % 2:2:2*size(xx,2)
        for i=1:length(grid_l)
            jcol_om_grid{i,1}=repmat(9+(1:3)+3*(i-1),2*grid_l(i),1);
            jrow_om_grid{i,1}=repmat((1:2*grid_l(i))+2*sum(grid_l(1:i-1)),3,1)';
        
            jcol_T_grid{i,1}=jcol_om_grid{i,1}+3*length(grid_l);
            jcol_T_grid{i,1}(2:2:size(jcol_T_grid{i,1},1),1)=NaN;
            jcol_T_grid{i,1}(1:2:size(jcol_T_grid{i,1},1),2)=NaN;
        
            jrow_T_grid{i,1}=jrow_om_grid{i,1};
            jrow_T_grid{i,1}(2:2:size(jcol_T_grid{i,1},1),1)=NaN;
            jrow_T_grid{i,1}(1:2:size(jcol_T_grid{i,1},1),2)=NaN;
        
        end
        
        jcol_om_grid=cell2mat(jcol_om_grid);
        jrow_om_grid=cell2mat(jrow_om_grid);
        
        jcol_T_grid=cell2mat(jcol_T_grid);
        jrow_T_grid=cell2mat(jrow_T_grid);
        
        jcol_T_grid(isnan(jcol_T_grid))=[];
        jrow_T_grid(isnan(jrow_T_grid))=[];
        
        jcol=[jcol_cam(:);jcol_om_grid(:);jcol_T_grid(:)];
        jrow=[jrow_cam(:);jrow_om_grid(:);jrow_T_grid(:)];
        
        JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,9+6*length(grid_l));
        
        % %Nonlinear least squares optimization routine for bundle adjustment
        options=optimset('lsqnonlin');
        options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'TolFun',1e-6,'TolX', 1e-6,'PrecondBandWidth', 2);
        [chi,resnorm,residual,exitflag,output]=lsqnonlin(@(chi) function_single_cam_cc_kc(xx,chi,X_cell),chi0,[],[],options);
        
        % options=optimset('lsqnonlin');
        % options=optimset(options,'Display','off','MaxIter',150,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam);%,'OutputFcn', @outfun_single_cam);
        % [chi2,resnorm2,residual2,exitflag,output,lambda,jacob]=lsqnonlin(@(chi) function_single_cam_cc_dist(xx,chi,X_cell),chi0,[],[],options);
        % % 
        % spy(JStr,'r',25); axis equal
        % hold on
        % figure;
        % spy(jacob,'b',11); axis equal
        % return

        if runs==1

            x_cell_init_length=length(x_cell);

            for i=length(x_cell):-1:1
        
                if ~isempty(x_cell{i})
            
                    xx_temp=x_cell{i};
            
                    F=function_single_cam_cc_kc(xx_temp,chi(:,[1:3 i+3 x_cell_init_length+i+3]),X_cell(i));
        
                    RR=hypot(F(1,:),F(2,:));
                    
                    MSE=mean(RR);
        
                    if MSE>3
        
                        changed=1;
        
                        x_cell(i)=[];
                        X_cell(i)=[];

                        om_use(:,i)=[];
                        T_use(:,i)=[];

                        grid_l(i)=[];

                        x_cam{use_grid_ind(i)}=NaN(size(x_cam{i}));
                        removed(use_grid_ind(i))=1;
                        use_grid(use_grid_ind(i))=[];

                        use_grid_ind(i)=[];
        
                    end
        
                    
                end
        
            end
        end
    end

end
    
fc_opt=chi([1 1],1);
cc_opt=chi([1 2],2);

kc_opt=[chi(1:2,3);0;0;0];

om_opt=NaN(3,length(x_cam));
T_opt=NaN(3,length(x_cam));

om_opt(:,use_grid_ind)=chi(:,4:size(X_cell,2)+3);
T_opt(:,use_grid_ind)=chi(:,size(X_cell,2)+4:end);
