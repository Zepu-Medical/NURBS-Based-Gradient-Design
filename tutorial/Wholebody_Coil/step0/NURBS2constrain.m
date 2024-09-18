function [id, coeff_wire,Mat_bsp_bas, Mat_dbsp_bas,gauss_weight, l_inteval, r_inteval, ind_eq,ind_lneq,id_wire,ind_op,ind_nop] = NURBS2constrain_v01(NURBS_curve, Wire_Id_ctrl,Type)
% Copyright (c) 2024 Zepu Medical. 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at

%    http://www.apache.org/licenses/LICENSE-2.0

% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% 
l_inteval=[];
r_inteval=[];
coeff_wire=[];
for c1=1:5
ind_eq{c1,1}=[];
end
for c1=2
ind_lneq{c1,1}=[];
ind_lneq{c1,2}=[];
end
id=[];
l0=0;
l1=0;
ind_op=[];
ind_nop=[];
id2=[1;0];
id_wire=cell(0,1);

for c1=1:length(NURBS_curve)
    if strcmpi(Type,'Remove')
        ind_opt{c1,1}=~any(NURBS_curve(c1).WireId(:)==Wire_Id_ctrl{c1}(:).',2);
    else
        ind_opt{c1,1}=any(NURBS_curve(c1).WireId(:)==Wire_Id_ctrl{c1}(:).',2);
    end
    if any(ind_opt{c1,1})
    Mat_bsp_bas{c1,1}=matlab.internal.math.blkdiag(NURBS_curve(c1).Bsp_Bas{ind_opt{c1}});
    Mat_dbsp_bas{c1,1}=matlab.internal.math.blkdiag(NURBS_curve(c1).dBsp_Bas{ind_opt{c1}});
    ncp_surf(c1)=size(Mat_bsp_bas{c1,1},2);
    cell_coeff0{c1,1}=cell2mat(NURBS_curve(c1).ContrPoint(ind_opt{c1}));
    cell_coeff0{c1,2}=cell2mat(NURBS_curve(c1).Weight(ind_opt{c1}));
    gauss_weight{c1,1}=cell2mat(NURBS_curve(c1).GaussWeight(ind_opt{c1}));
    end
end

% coeff_all=cell2mat(cell_coeff0);


for c1=1:length(NURBS_curve)
    for c2=1:length(NURBS_curve(c1).Knot)
        n=NURBS_curve(c1).Knot{c2}(end);
        l1=l1+2*n+2; 

        if strcmpi(Type,'Remove')
            cond_1=all(NURBS_curve(c1).WireId(c2)~=Wire_Id_ctrl{c1}(:).');
        else
            cond_1=any(NURBS_curve(c1).WireId(c2)==Wire_Id_ctrl{c1}(:).');
        end

        if cond_1

        ind_op=[ind_op;(l1-3*n:l1).'];

        id(end+1)=length(NURBS_curve(c1).Weight{c2});

        
        if c1==2||(NURBS_curve(c1).WireId(c2,1)==[0])
        id_wire{end+1,1}=sum(id(1:end-1))+(1:id(end)).';  
        
        end
        X=NURBS_curve(c1).ContrPoint{c2}(:,1);
        Y=NURBS_curve(c1).ContrPoint{c2}(:,2);
        W=NURBS_curve(c1).Weight{c2};

        coeff_wire=[coeff_wire;X,Y,W];


        ind_eq{1,1}=[ind_eq{1,1},l0+3*(1:n-1)]; % 
        %
            if NURBS_curve(c1).ContrPoint{c2}(1,1)<1e-4
            ind_eq{2,1}=[ind_eq{2,1},l0+1];  %  
            end

            if NURBS_curve(c1).ContrPoint{c2}(end,1)<1e-4
            ind_eq{2,1}=[ind_eq{2,1},l0+2*n+1]; 
            end
        if c1>=1&c1<=length(NURBS_curve)-1
        if any(NURBS_curve(c1).WireId(c2,1)==NURBS_curve(c1+1).WireId.')
            if NURBS_curve(c1).ContrPoint{c2}(1,2)>(1-1e-4)
            ind_eq{3,1}=[ind_eq{3,1},l0+1];  %  
            elseif NURBS_curve(c1).ContrPoint{c2}(end,2)>(1-1e-4)
            ind_eq{3,1}=[ind_eq{3,1},l0+2*n+2];  %  
            end
        end
        end
        if c1>=2&c1<=length(NURBS_curve)
        if any(NURBS_curve(c1).WireId(c2,1)==NURBS_curve(c1-1).WireId.')

            if NURBS_curve(c1).ContrPoint{c2}(1,2)<4e-3
            ind_eq{4,1}=[ind_eq{4,1},l0+1];  %  
            elseif NURBS_curve(c1).ContrPoint{c2}(end,2)<4e-3
            ind_eq{4,1}=[ind_eq{4,1},l0+2*n+2];  %  
            end
            
        end
        end

       

            l0=l0+2*n+2;
        if c1==1||c1==2  %
            temp_lhx=X-0.02; %0.1
            temp_rhx=X+0.02;
            temp_lhy=Y-0.05; %0.04
            temp_rhy=Y+0.05;

             if c1==1 && any(NURBS_curve(c1).WireId(c2,1)==[])
                temp_rhy(Y>0.0&Y<0.12)=temp_rhy(Y>0.0&Y<0.12)+0.015;
             end

%              if (NURBS_curve(c1).WireId(c2,1)==[26])
%                  
%                  temp_lhy(X>=0&X<1.1&Y>0.9)=0.99;
%              end
             if (NURBS_curve(c1).WireId(c2,1)==[])
                 temp_rhy(X>=0&X<0.98&temp_rhy>0.988)=0.988;
             end
             

            if c1==1
                temp_lhy(temp_lhy<0.007)=0.007;
                temp_rhy(temp_rhy>0.993)=0.993;
            temp_lhx(temp_lhx<0)=0;
            temp_rhx(temp_rhx>pi/2-0.01)=pi/2-0.01;
            elseif c1==2
                temp_lhy(temp_lhy<0.005)=0.005;
                temp_rhy(temp_rhy>0.995)=0.995;
            temp_lhx(temp_lhx<0)=0;
            temp_rhx(temp_rhx>pi/2-0.01)=pi/2-0.01;
            end


        end



            if NURBS_curve(c1).ContrPoint{c2}(1,2)>(1-1e-4)
            temp_lhy(1)=1;
            temp_rhy(1)=1;
            elseif NURBS_curve(c1).ContrPoint{c2}(1,2)<1e-4
                temp_lhy(1)=0;
                temp_rhy(1)=0;
            elseif NURBS_curve(c1).ContrPoint{c2}(1,1)<(1e-4)
            temp_lhx(1)=0;
            temp_rhx(1)=0;    
            end



            if NURBS_curve(c1).ContrPoint{c2}(end,2)>(1-1e-4)
            temp_lhy(end)=1;
            temp_rhy(end)=1;
            elseif NURBS_curve(c1).ContrPoint{c2}(end,2)<(1e-4)
            temp_lhy(end)=0;
            temp_rhy(end)=0;
            elseif NURBS_curve(c1).ContrPoint{c2}(end,1)<(1e-4)
            temp_lhx(end)=0;
            temp_rhx(end)=0;
            end
        




        temp_lhw=W-0.1; %0.06
        temp_rhw=W+0.1;
        temp_lhw(temp_lhw<-0.4)=-0.4;
        temp_rhw(temp_rhw>0.5)=0.5;

        % temp_rhw(1:3:end)=0;
        % temp_lhw(1:3:end)=0;
        temp_rhw(1)=0;
        temp_lhw(1)=0;
        temp_rhw(end)=0;
        temp_lhw(end)=0;

        l_inteval=[l_inteval;temp_lhx,temp_lhy,temp_lhw];
        r_inteval=[r_inteval;temp_rhx,temp_rhy,temp_rhw];
        else
            ind_nop=[ind_nop;(l1-3*n:l1).'];
        end

    end
    
end


end