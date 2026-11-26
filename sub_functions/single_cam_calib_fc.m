function [fc_opt, T_opt, om_opt, x_cam, removed, residual]=single_cam_calib_fc(fc,cc,kc,T,om,x_cam,XX_cam)
% save temp_single_cam
% moo=1
% return
% clear 
% close all force
% load temp_single_cam

% clear disp
% disp='off';

% use=[1:3];
% x_cam=x_cam(use);
% XX_cam=XX_cam(use);
% 
% x_cam{1}(:,3:5)=NaN;
% x_cam{2}(:,8)=NaN;
% 
% T=T(:,use);
% om=om(:,use);

use_grid=zeros(size(x_cam));

k=1;

om_use=[];
T_use=[];

for i=1:length(x_cam)
    
    xx=x_cam{i};
    
    if ~isempty(xx) && any(~isnan(xx(1,:)))

        [Ay, Ax]=find(isnan(xx));
        xx(:,Ax)=[];

        X=XX_cam;
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

%BUNDLE ADJUSTMENT ROUTINE
%Specifies structure of the Jacobian matrix in the sparse matrix JStr

changed=0;
removed=zeros(size(x_cam));

for runs=1:2

    if runs==1 || (runs==2 && changed==1)
    
        xx=cell2mat(x_cell);
        grid_t=numel(xx);
    
        clear j*
     
        chi0=[[fc(1);0;0] om_use T_use];
        jcol_cam=repmat([1 1],size(xx,2),1);
    
        jcol_cam=jcol_cam(:);
        jrow_cam=[1:2:2*size(xx,2) 2:2:2*size(xx,2)]';
    
        for i=1:length(grid_l)
            jcol_om_grid{i,1}=repmat(3+(1:3)+3*(i-1),2*grid_l(i),1);
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
    
        JStr=sparse(jrow,jcol,ones(size(jrow)),grid_t,3+6*length(grid_l));
    
        % %Nonlinear least squares optimization routine for bundle adjustment
        options=optimset('lsqnonlin');
        
        options=optimset(options,'JacobPattern',JStr,'Display','off','MaxIter',150,'TolFun',1e-4,'TolX', 1e-4,'MaxFunEvals',200*length(xx),'PrecondBandWidth', 2);
        [chi,resnorm,residual,exitflag,output]=lsqnonlin(@(chi) function_single_cam_fc(xx,chi,X_cell,cc,kc),chi0,[],[],options);
    
        % options=optimset('lsqnonlin');
        % options=optimset(options,'Display','off','MaxIter',150,'TolFun',1e-5,'TolX', 1e-5,'MaxFunEvals',200*length(xx),'OutputFcn', @outfun_single_cam,'PrecondBandWidth', 2);
        % [chi2,resnorm2,residual,exitflag,output,lambda,jacob]=lsqnonlin(@(chi) function_single_cam_no_cc(xx,chi,X_cell,cc,kc),chi0,[],[],options);
        % return
        % disp(resnorm)
                
        if runs==1

            x_cell_init_length=length(x_cell);

             for i=length(x_cell):-1:1
        
                if ~isempty(x_cell{i})
            
                    xx_temp=x_cell{i};
            
                    F=function_single_cam_fc(xx_temp,chi(:,[1 i+1 x_cell_init_length+i+1]),X_cell(i),cc,kc);
        
                    RR=hypot(F(1,:),F(2,:));
                    
                    MSE=mean(RR)
                 
                    if MSE>5
        
                        changed=1;
                               
                        x_cell(i)=[];
                        X_cell(i)=[];

                        om_use(:,i)=[];
                        T_use(:,i)=[];

                        grid_l(i)=[];

                        use_grid_ind(i);

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

om_opt=NaN(3,length(x_cam));
T_opt=NaN(3,length(x_cam));

om_opt(:,use_grid_ind)=chi(:,2:size(X_cell,2)+1);
T_opt(:,use_grid_ind)=chi(:,size(X_cell,2)+2:end);