function [cneq,ceq,gcneq,gceq]=Constrain_fieldv03(coeff,WirePara,OptPara,ConsPara)
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

num_out=nargout;

N=OptPara.Mat_bsp_bas_d;

dN=OptPara.Mat_dbsp_bas_d;

gw=OptPara.gauss_weight_d;
n_total_wires=length(N);

for cl0=1:n_total_wires
    N{cl0,1}=N{cl0,1};
    dN{cl0,1}=dN{cl0,1};
    gw{cl0,1}=gw{cl0,1}.'/2;  
end

size_coil=WirePara.size_coil;

target_point=ConsPara.Field;
target_point_G=ConsPara.Gradient;

cneq=[];
ceq=[];
gceq=[];
gcneq=[];


if strcmpi(WirePara.fold_condition,'8-fold')
    n_rep=8;
elseif strcmpi(WirePara.fold_condition,'4-fold')
    n_rep=4;
end




%% 不等式约束

if ~isempty(target_point)|~isempty(target_point_G)

    B=zeros(size(target_point,1)/n_rep,1);
    A=zeros(2*size(target_point,1)/n_rep,1);
    G=zeros(size(target_point_G,1)/n_rep+1,1);


    if strcmpi(ConsPara.FieldComponent,'Br')
        StrayFieldType=["Bx","By"];
    elseif strcmpi(ConsPara.FieldComponent,'Bz')
        StrayFieldType="Bz";
    end
    if strcmpi(ConsPara.FieldVP,'on')
        StrayFieldType=[StrayFieldType,"Ax","Ay","Az"];
    end
    if ~isfield(ConsPara,'TargetFieldType')
        ConsPara.TargetFieldType='Gradient';
    end

    if strcmpi(ConsPara.TargetFieldType,'Gradient') % 'Gradient' or 'Bz'
        TargetFieldType="Gradient";
    elseif strcmpi(ConsPara.TargetFieldType,'Bz')
        TargetFieldType="Bz";
    end

    % 二维平面的参数曲线
    if strcmpi(WirePara.fold_condition,'8-fold')
        n_rep=8;
    elseif strcmpi(WirePara.fold_condition,'4-fold')
        n_rep=4;
    end
    nf=size(target_point,1)/n_rep;
    ng=size(target_point_G,1)/n_rep;
    cell_coeff=mat2cell(coeff,OptPara.Cp_Num,[1,1,1]); %曲线控制点及权重
    phi=atan2(target_point(1:nf,2),target_point(1:nf,1));
    proj1=cos(phi);
    proj2=sin(phi);
    gBc=cell(n_total_wires,3);
    gGc=cell(n_total_wires,3);
    gAc=cell(n_total_wires,3);
    gTorquec=cell(n_total_wires,3);
    Torquec=cell(n_total_wires,1);
        if num_out>2
            NURBSGeneOrder=1;
        else
            NURBSGeneOrder=0;
        end
    parfor c0=1:n_total_wires  %parfor

        if OptPara.Cp_Num(c0)>0

            pu=cell_coeff{c0,1};
            pv=cell_coeff{c0,2};
            w=exp(cell_coeff{c0,3});
            N_tmp=N{c0};
            dN_tmp=dN{c0};

            size1=size(N_tmp);
            [C2,C3]=NURBSgenerator({pu,pv,w},{N_tmp,dN_tmp,[]},size_coil(c0,:),NURBSGeneOrder,1);
            C3.gw=gw{c0};

            SF_val=NURBSfield(target_point,C3,StrayFieldType,[0,1]);
            TF_val=NURBSfield(target_point_G,C3,TargetFieldType,[0,1]);

            G0_val=NURBSfield([0,0,0],C3,"Gradient",[0,1]);


            if strcmpi(ConsPara.TargetFieldType,'Gradient')
                G_t=sum(reshape(TF_val{1,1}.',[ng,n_rep]),2);
                gG_t=sum(reshape(TF_val{1,2},[3*size1(1),ng,n_rep]),3);

            elseif strcmpi(ConsPara.TargetFieldType,'Bz')
                coeffz=[1	-1	-1	1	1	-1	-1	1].';
                G_t=reshape(TF_val{1,1},[ng,n_rep])*coeffz(1:n_rep);
                gG_t=reshape(reshape(TF_val{1,2},[3*size1(1)*ng,n_rep])*coeffz(1:n_rep),[3*size1(1),ng]);
            end
                G_t=[G0_val{1,1}*n_rep;G_t];
                gG_t=[G0_val{1,2}*n_rep,gG_t];


            if strcmpi(ConsPara.FieldComponent,'Br')
                coeffx=[1	1	1	1	-1	-1	-1	-1].';
                coeffy=[1	-1	1	-1	-1	1	-1	1].';
                Bx_t=reshape(SF_val{1,1},[nf,n_rep])*coeffx(1:n_rep,1);
                By_t=reshape(SF_val{2,1},[nf,n_rep])*coeffy(1:n_rep,1);
                gBx_t=reshape(reshape(SF_val{1,2},[3*size1(1)*nf,n_rep])*coeffx(1:n_rep,1),[3*size1(1),nf]);
                gBy_t=reshape(reshape(SF_val{2,2},[3*size1(1)*nf,n_rep])*coeffy(1:n_rep,1),[3*size1(1),nf]);
                B_t=Bx_t.*proj1+By_t.*proj2;
                gB_t=gBx_t.*proj1.'+gBy_t.*proj2.';
            elseif strcmpi(ConsPara.FieldComponent,'Bz')
                coeffz=[1	-1	-1	1	1	-1	-1	1].';
                B_t=reshape(SF_val{1,1},[nf,n_rep])*coeffz(1:n_rep,1);
                gB_t=reshape(reshape(SF_val{1,2},[3*size1(1)*nf,n_rep])*coeffz(1:n_rep,1),[3*size1(1),nf]);
            end
            if strcmpi(ConsPara.FieldVP,'on')
                A_val=SF_val(end-2:end,1:2);
                coeffAx=[1	-1	1	-1	1	-1	1	-1].';
                coeffAy=[1	1	1	1	1	1	1	1].';
                coeffAz=[1	1	-1	-1	-1	-1	1	1].';

                Ax_t=reshape(A_val{1,1},[nf,n_rep])*coeffAx(1:n_rep,1);
                Ay_t=reshape(A_val{2,1},[nf,n_rep])*coeffAy(1:n_rep,1);
                Az_t=reshape(A_val{3,1},[nf,n_rep])*coeffAz(1:n_rep,1);

                gAx_t=reshape(reshape(A_val{1,2},[3*size1(1)*nf,n_rep])*coeffAx(1:n_rep,1),[3*size1(1),nf]);
                gAy_t=reshape(reshape(A_val{2,2},[3*size1(1)*nf,n_rep])*coeffAy(1:n_rep,1),[3*size1(1),nf]);
                gAz_t=reshape(reshape(A_val{3,2},[3*size1(1)*nf,n_rep])*coeffAz(1:n_rep,1),[3*size1(1),nf]);

                A_t=[Ax_t.*(-proj2)+Ay_t.*proj1;Az_t];
                gA_t=[gAx_t.*(-proj2.')+gAy_t.*proj1.',gAz_t];
                A=A+A_t;
                gAc(c0,:)={gA_t(1:size1(1),:),gA_t((1:size1(1))+size1(1),:),gA_t((1:size1(1))+size1(1)*2,:)};
            end

            B=B+B_t;
            G=G+G_t;

            gBc(c0,:)={gB_t(1:size1(1),:),gB_t((1:size1(1))+size1(1),:),gB_t((1:size1(1))+size1(1)*2,:)};
            gGc(c0,:)={gG_t(1:size1(1),:),gG_t((1:size1(1))+size1(1),:),gG_t((1:size1(1))+size1(1)*2,:)};

        end


        if strcmpi(WirePara.fold_condition,'4-fold')
           Torque_t=975*sum(C3.dCy.*C3.gw.*(C3.Cz-0));  %
           % Torque_t=975*sum(C3.dCy.*C3.gw);  %
            Torquec{c0,1}=Torque_t;
        end


        if strcmpi(WirePara.fold_condition,'4-fold')
            gTorque_t=975*sum(C3.gdCy.*(C3.gw.*(C3.Cz-0))+(C3.dCy.*C3.gw).*C3.gCz,2);
            gTorquec(c0,:)={gTorque_t(1:size1(1),1),gTorque_t((1:size1(1))+size1(1),1),gTorque_t((1:size1(1))+size1(1)*2,1)};
        end
    end



    %          out=B;
    %          dout=dB_dpu_a;


    %%
    if ~strcmpi(ConsPara.FieldType,'Uniform')
        B=(B.'*ConsPara.transformer).';
    end
    if strcmpi(ConsPara.FieldVP,'on')
        A=(A.'*ConsPara.transformerA).';
    end

    cneq=[cneq;B-ConsPara.Field_r;ConsPara.Field_l-B;   
        G-ConsPara.Gradient_r;ConsPara.Gradient_l-G];
    if strcmpi(WirePara.fold_condition,'4-fold')   
        Torque=sum(cell2mat(Torquec));
        cneq=[cneq;Torque-ConsPara.Torque_r;ConsPara.Torque_l-Torque]; 
    end
    if num_out>2   
        gB=cell2mat(gBc(:));
        gG=cell2mat(gGc(:));
    if strcmpi(ConsPara.FieldType,'Uniform')

    else 
        gB=(gB*ConsPara.transformer);
    end
        gcneq=[gcneq,gB,-gB,gG,-gG];
    if strcmpi(WirePara.fold_condition,'4-fold')
        gTorque=cell2mat(gTorquec(:));
        gcneq=[gcneq,gTorque,-gTorque];
    end
    end
    if strcmpi(ConsPara.FieldVP,'on')
        cneq=[cneq;A-ConsPara.A_r;ConsPara.A_l-A];
        if num_out>2     
        gA=cell2mat(gAc(:));
        gA=(gA*ConsPara.transformerA);
        gcneq=[gcneq,gA,-gA];
        end
    end



end


end



