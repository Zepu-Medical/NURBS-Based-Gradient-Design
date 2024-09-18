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


%%
  load .\Coil2D_ini.mat
 


sz=[(size_c(:,1)+size_c(:,3))/2,abs(size_c(:,2)-size_c(:,4))];  
alpha_r={0.01*ones(length(Coil2D_ini(1).Wire2D),1),0.01*ones(length(Coil2D_ini(2).Wire2D),1),0.01*ones(length(Coil2D_ini(3).Wire2D),1)};
alpha_r{2}(end)=0.01;
f_eval=[];
for c0=1:length(Coil2D_ini)
    Wire2D_manu=Coil2D_ini(c0).Wire2D;
    count_wire=1;

    l1=length(Coil2D_ini(c0).Wire2D);
            tp_Xi=cell(l1,1);
            tp_Knot=cell(l1,1);
            tp_ContrPoint=cell(l1,1);
            tp_Weight=cell(l1,1);
            tp_WireId=zeros(l1,1);
    parfor c1=1:length(Wire2D_manu)   % parfor
        wire_fit_tmp=Wire2D_manu{c1}.*sz(c0,:);
        
        seg_tmp=wire_fit_tmp(2:end,:)-wire_fit_tmp(1:end-1,:);
        cum_length=[0;cumsum(vecnorm(seg_tmp,2,2))];
        knot_idx=Coil2D_ini(c0).seg_idx{c1,1};
        
        Npara=size(wire_fit_tmp,1);
        xi0=zeros(Npara,1);

        cp_ini=[1];
        Knots=[0,0,0:0.5:length(knot_idx)-1,length(knot_idx)-1,length(knot_idx)-1];

        Ncp=2*length(knot_idx);
        for c2=1:length(knot_idx)-1
        idx_tmp=knot_idx(c2):knot_idx(c2+1);
        xi0(idx_tmp)=(c2-1)+(cum_length(idx_tmp)-cum_length(knot_idx(c2)))/(cum_length(knot_idx(c2+1))-cum_length(knot_idx(c2)));
        cp_ini_tmp=linspace(knot_idx(c2),knot_idx(c2+1),4);
        cp_ini=[cp_ini,round(cp_ini_tmp(2:3))];
        end
        cp_ini(1,end+1)=knot_idx(end);

        coeffs=[wire_fit_tmp(cp_ini,:),zeros(Ncp,1)];
        
        coeff_ini=[xi0;coeffs(:)];

        A=spdiags(ones(Npara-1,1),0,Npara-1,length(coeff_ini))+spdiags(-ones(Npara-1,1),1,Npara-1,length(coeff_ini));
        b=sparse(Npara-1,1);
        
        n=Knots(end);
        ind_eq1=3*(1:n-1); % 
        ind_eq2=[1,3*n];  % 
        id_tmp=[];
        if abs(wire_fit_tmp(1,1))>1e-2
            id_tmp=[id_tmp,1];
        end

        if abs(wire_fit_tmp(end,1))>1e-2
            id_tmp=[id_tmp,2];
        end
        ind_eq2(id_tmp)=[];

        ind_eq=cell(0,1);
        ind_eq{1,1}=ind_eq1;
        ind_eq{2,1}=ind_eq2;
        lb=-inf*ones(length(xi0),1);
        rb=inf*ones(length(xi0),1);

        lb(knot_idx,1)=[0,-0+(1:length(knot_idx)-1)];
        rb(knot_idx,1)=[(0:length(knot_idx)-2)+0,length(knot_idx)-1];

         w_lb=-0.2*ones(Ncp,1);
         w_rb=0.3*ones(Ncp,1);
        w_lb(1)=0;
        w_rb(1)=0;
        w_lb(end)=0;
        w_rb(end)=0;
         lb=[lb;-inf*ones(Ncp,1);-inf*ones(Ncp,1);w_lb];
         rb=[rb;inf*(pi/2+0.01)*sz(c0,1)*ones(Ncp,1);inf*sz(c0,2)*ones(Ncp,1)+0.03;w_rb];

        condition=2;
                    xi_c0=linspace(0,1,3);
            xi_c0(end)=[];
            xi_c=[reshape(xi_c0(:)+(0:Knots(end)-1),[],1);Knots(end)];
        xi_span=Func_ParaSpan(xi_c,Knots);
    
        xi2cp=xi_span-flip(0:2)+1; %
        [bsp_bas,dbsp_bas,d2bsp_bas]= Func_Basisfun(xi_span, xi_c, 2, Knots);
        l_xi=size(xi2cp,1);
        ConsPara=struct();
        ConsPara.N_dC=sparse(xi2cp,repmat(1:l_xi,[1,2+1]),bsp_bas{1,2+1}(:),Ncp,l_xi);
        ConsPara.dN_dC=sparse(xi2cp,repmat(1:l_xi,[1,2+1]),dbsp_bas{1,2+1}(:),Ncp,l_xi);
        ConsPara.d2N_dC=sparse(xi2cp,repmat(1:l_xi,[1,2+1]),d2bsp_bas{1,2+1}(:),Ncp,l_xi);

        ConsPara.curvature=[-600*ones(l_xi,1),600*ones(l_xi,1)];
        ConsPara.curvature=[];
        FitPara=struct();
        FitPara.Knots=Knots;
        FitPara.WirePoint=wire_fit_tmp;
        FitPara.condition=condition;
        FitPara.alpha=alpha_r{c0}(c1);
        options=optimoptions('fmincon','MaxFunctionEvaluations',4e5,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true,'UseParallel',true,'Display','none',...
            'HessianFcn',@(coeff,lambda)Hessian_fit(coeff,lambda,FitPara,ConsPara)); 
        if isempty(ConsPara.curvature)
        [coeff_fit]=fmincon(@(coeff)Obj_NURBSfitting(coeff,FitPara),coeff_ini,A,b,[],[],lb,rb,[],options);
        else
        [coeff_fit]=fmincon(@(coeff)Obj_NURBSfitting(coeff,FitPara),coeff_ini,A,b,[],[],lb,rb,@(coeff)Constrain_NURBSfitting(coeff,FitPara,ConsPara),options);
        end

        xi=coeff_fit(1:Npara,1);
        coeffs=reshape(coeff_fit(Npara+1:end,1),[Ncp,3]);
            tp_Xi{c1,1}=xi;
            tp_Knot{c1,1}=Knots(:);
            tp_ContrPoint{c1,1}=coeffs(:,1:2)./sz(c0,:);
            tp_Weight{c1,1}=coeffs(:,3);
            tp_WireId(c1,1)=Coil2D_ini(c0).WireId(c1);
    end
            NURBS_curve_fit0(c0,1).Xi=tp_Xi;
            NURBS_curve_fit0(c0,1).Knot=tp_Knot;
            NURBS_curve_fit0(c0,1).ContrPoint=tp_ContrPoint;
            NURBS_curve_fit0(c0,1).Weight=tp_Weight;
            NURBS_curve_fit0(c0,1).WireId=tp_WireId;
end






%%  

NURBS_curve_modified=NURBS_curve_fit0;
%  


WirePara.size_coil=size_c;  %


WirePara.fold_condition='8-fold';
WirePara.DSVr=[0.3;0.4;0.5];

WirePara.Shieldr=[0.450,0.450,0.7,-0.7];


for c1=1:length(NURBS_curve_modified)
    for c2=1:length(NURBS_curve_modified(c1).ContrPoint)
        for id=[1,size(NURBS_curve_modified(c1).ContrPoint{c2},1)]
        p1=NURBS_curve_modified(c1).ContrPoint{c2}(id,:);
        distance=abs([p1(1),p1(2)-0,p1(2)-1]);
        [~,index]=sort(distance,"ascend");
        if index(1)==1

            NURBS_curve_modified(c1).ContrPoint{c2}(id,1)=0;
            if id==1
                NURBS_curve_modified(c1).ContrPoint{c2}(1:2,2)=sum(NURBS_curve_modified(c1).ContrPoint{c2}(1:2,2))/2;
            else
                NURBS_curve_modified(c1).ContrPoint{c2}(end-1:end,2)=sum(NURBS_curve_modified(c1).ContrPoint{c2}(end-1:end,2))/2;
            end
            NURBS_curve_modified(c1).ContrPoint{c2}(id,1)=0;
            NURBS_curve_modified(c1).ContrPoint{c2}(id,1)=0;
        end
        if index(1)==2
            NURBS_curve_modified(c1).ContrPoint{c2}(id,2)=0;
        end
        if index(1)==3
            NURBS_curve_modified(c1).ContrPoint{c2}(id,2)=1;
        end

        end
    end
end

for c1=1:length(NURBS_curve_modified)
    for c2=1:length(NURBS_curve_modified(c1).ContrPoint)
        for id=[1,size(NURBS_curve_modified(c1).ContrPoint{c2},1)]
        p1=NURBS_curve_modified(c1).ContrPoint{c2}(id,:);
            if p1(1)==0||p1(2)==0||p1(2)==1

            else
                error();
            end

        end
    end
end




for c0=1:length(NURBS_curve_modified)

    for c1=1:length(NURBS_curve_modified(c0).ContrPoint)
    NURBS_curve_modified(c0).ContrPoint{c1,1}=abs(NURBS_curve_modified(c0).ContrPoint{c1,1});

    ind_tmp=NURBS_curve_modified(c0).ContrPoint{c1,1}(:,2)>1;
    NURBS_curve_modified(c0).ContrPoint{c1,1}(ind_tmp,2)=2-NURBS_curve_modified(c0).ContrPoint{c1,1}(ind_tmp,2);
    end
end



figure
for c1=1:length(Coil2D_ini)
subplot(1,length(Coil2D_ini),c1)
     Plot_Wire2D(Coil2D_ini(c1));
    NURBS_plot(NURBS_curve_modified(c1),{'-'},'control point','on');
    xlim([0,pi/2]);
end




NURBS_curve_ini=NURBS_curve_modified;

    
save('NURBS_curve_ini.mat','WirePara','NURBS_curve_ini');




