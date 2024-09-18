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


load('.\NURBS_curve_ini.mat','NURBS_curve_ini','WirePara');



%% 
FormPara.XiForm='gaussian';
NURBS_curve_ini = NURBS_formal(NURBS_curve_ini,FormPara);

FormPara1.XiForm='uniform';
NURBS_curve_dc = NURBS_formal(NURBS_curve_ini,FormPara1);


Wire_Id_ctrl=cell(1, length(NURBS_curve_ini));

c0=0;
    WirePara_d=WirePara;
    for c1=1:length(NURBS_curve_ini)
        for c2=1:length(NURBS_curve_ini(c1).Xi)
            c0=c0+1;

    Mat_bsp_bas_d{c0,1}=NURBS_curve_ini(c1).Bsp_Bas{c2};
    Mat_dbsp_bas_d{c0,1}=NURBS_curve_ini(c1).dBsp_Bas{c2};
    WirePara_d.size_coil(c0,:)=WirePara.size_coil(c1,:);
    gauss_weight_d{c0,1}=(NURBS_curve_ini(c1).GaussWeight{c2});
    cell_coeff0{c0,1}=NURBS_curve_ini(c1).ContrPoint{c2};
    cell_coeff0{c0,2}=NURBS_curve_ini(c1).Weight{c2};
    ConsPara.N_dC{c0,1}=NURBS_curve_dc(c1).Bsp_Bas{c2};
    ConsPara.dN_dC{c0,1}=NURBS_curve_dc(c1).dBsp_Bas{c2};
    ConsPara.d2N_dC{c0,1}=NURBS_curve_dc(c1).d2Bsp_Bas{c2};
        end
    end


coeff_all=cell2mat(cell_coeff0);
coeff_0=coeff_all;
%%


%% 设置目标场点部分
d_r_range=[0.008,0.004,0.00,-0.002,-0.004];

global  glo_coeff count_iter 
glo_coeff=cell(0,1);
count_iter=0;
 % clear global variable


ConsPara.TargetFieldType='Bz';  %Bz or Gradient

g0=482;
Bs=12;
if strcmpi(ConsPara.TargetFieldType,'Gradient') 
       target_para=[
           0.0,0.0,0.17,2,2;
           0,0,0.1,2,2;
        0.180,0.180,0,2,2;
        0.09,0.09,0.13,2,2
        0.180,0.180,0.07,2,2;
      %   0.24,0.24,0.048+0.006,2,2;
        ];

d_r=0.02;  %0.022 0.016
Gbias_r=d_r+[0.005,0.005,0.002,0.002,0.005];
Gbase_r=[0.987,1.001,0.98,1.00,0.975];
elseif strcmpi(ConsPara.TargetFieldType,'Bz') 

R_DSV=0.25; 
 
theta=(5.6:10.55:174.4).'/180*pi;
target_para=[R_DSV*sin(theta),R_DSV*sin(theta),R_DSV*cos(theta),2*ones(length(theta),1),2*ones(length(theta),1)
        ];
Gbias_r=0.0646*ones(1,size(target_para,1));   %0.049
Gbase_r=1*ones(1,size(target_para,1));
end

dfluc=[0.1,0.1,0.1]; % for field constrain


torque_range0=[-40,40;];
    torque_range=torque_range0+[10,-10]*4+[-0.3,0.3];

NURBS_curve1 = coeff2NURBS(coeff_0, NURBS_curve_ini);

field_point={};

    for c1=1:size(target_para,1)
    if target_para(c1,1)>0 & target_para(c1,2)>0
    TargetPara.Region=target_para(c1,1:3);
    TargetPara.SurfaceShape='circle';
    TargetPara.SegNum=[32,64];
    TargetPara.DeformPara=[2,2];
    TargetPara.MirrorSymmetry=1;
    TargetPara.RotateSymmetry='4-fold';

        field_point{end+1,1}=Func_FieldRegion(TargetPara);
    else
        field_point{end+1,1}=target_para(c1,1:3);
    end
    end


    ConsPara.Gradient=UnfoldPoint(cell2mat(field_point), WirePara);

    
    neq_G={g0*(1-0.004);g0*(1+0.004)};

    if strcmpi(ConsPara.TargetFieldType,'Gradient')
    Gbias=Gbias_r;
    Gbase=Gbase_r;

    
    for c1=1:length(field_point)
        neq_G{1,1}=[neq_G{1,1};g0*ones(size(field_point{c1},1),1)*(Gbase(c1)-Gbias(c1))];
        neq_G{2,1}=[neq_G{2,1};g0*ones(size(field_point{c1},1),1)*(Gbase(c1)+Gbias(c1))];
    end

    elseif strcmpi(ConsPara.TargetFieldType,'Bz')
    Gbias=Gbias_r;
    Gbase=Gbase_r;

    for c1=1:length(field_point)
        neq_G{1,1}=[neq_G{1,1};g0*field_point{c1}(:,1)*Gbase(c1)-R_DSV*g0*Gbias(c1)];
        neq_G{2,1}=[neq_G{2,1};g0*field_point{c1}(:,1)*Gbase(c1)+R_DSV*g0*Gbias(c1)];
    end
    end
    
    ConsPara.Gradient_l=neq_G{1};
    ConsPara.Gradient_r=neq_G{2};
    ConsPara.Torque_l=torque_range(1,1);
    ConsPara.Torque_r=torque_range(1,2);

    ConsPara.FieldComponent='Bz'; % Bz, Br, 'A'
    ConsPara.FieldType='Uniform'; %'Uniform'
    ConsPara.FieldVP='on'; %

        ConsPara.field_size=[9,32];
        n_sum_tmp=2;     

if  strcmpi(ConsPara.FieldVP,'on')
        m_range=1:2:11;
        r_n=0.450+(0.5:60)*1e-4; 
        l_z=1.4;
        phi_t=(0:1)/4*pi;
        r_t=0.25;
        z_t=[0.25,0];
        rang.z=linspace(0,1,ConsPara.field_size(2)*2-1);
        rang.n=1:20;
        rang.phi=(0.5:ConsPara.field_size(1))/ConsPara.field_size(1)*pi/2;

        U_eddy=real(GeneUeddy(rang, phi_t, z_t,r_t, m_range, r_n,l_z)).';
        U_eddy_phi=U_eddy(1:end/2,:);
        U_eddy_z=U_eddy(end/2+1:end,:);
        id_tmp=(1:ConsPara.field_size(1)).'+ConsPara.field_size(1)*(ConsPara.field_size(2)*2-2:-1:ConsPara.field_size(2));
        if strcmpi(WirePara.fold_condition,'8-fold')
        U_eddy_phi=U_eddy_phi(1:ConsPara.field_size(1)*(ConsPara.field_size(2)),:)+[U_eddy_phi(id_tmp(:),:);zeros(ConsPara.field_size(1),size(U_eddy,2))];
        U_eddy_z=-U_eddy_z(1:ConsPara.field_size(1)*(ConsPara.field_size(2)),:)+[U_eddy_z(id_tmp(:),:);zeros(ConsPara.field_size(1),size(U_eddy,2))];
        end
        ConsPara.transformerA=[U_eddy_phi;U_eddy_z];
        ConsPara.A_l=-0.002/4*repmat(g0*r_t*cos(phi_t*0).',[length(z_t),1]);
        ConsPara.A_r=abs(ConsPara.A_l);
    end
    [PPhi,ZZ]=ndgrid(rang.phi,linspace(l_z/2,0.0,ConsPara.field_size(2)));%32

    shield_point=[0.450*cos(PPhi(:)),0.450*sin(PPhi(:)),ZZ(:)];
    ConsPara.Field=UnfoldPoint(shield_point, WirePara);


  Cell_wire1 = nurbs_fold2all(NURBS_curve1,WirePara.size_coil,WirePara.fold_condition);  

    Cell_shield{1,2}=shield_point; %
    Cell_shield{1,1}=[]; %
    Cell_shield{1,3}='Cartesian';
    fieldcomponent=[1,2,3];
   [Cell_shield] = Func_MagneticField(Cell_shield,Cell_wire1,fieldcomponent);
   if strcmpi(ConsPara.FieldComponent,'Bz')
    B=1e7*Cell_shield{1}(:,3);
   elseif strcmpi(ConsPara.FieldComponent,'Br')
    B=1e7*(Cell_shield{1}(:,1).*cos(PPhi(:))+Cell_shield{1}(:,2).*sin(PPhi(:)));
   end
   Bz_cell=mat2cell(abs(B).',1,ConsPara.field_size(1)*ones(ConsPara.field_size(2),1));

    if strcmpi(ConsPara.FieldType,'Uniform')
        B_s=2*abs(Cell_shield{1}(:,3)*1e7);
        B_s(B_s<0.5)=0.5;
        B_s=cellfun(@(x)max(abs(x)),Bz_cell);
        B_s=(kron(B_s,ones(1,ConsPara.field_size(1)))).';
        ConsPara.Field_l=-0.8*abs(B_s);
        ConsPara.Field_r=0.8*abs(B_s);
        ConsPara.Field_l=-Bs(1)*ones(size(shield_point,1),1);
        ConsPara.Field_r=Bs(1)*ones(size(shield_point,1),1);
        id_tmp=ZZ>0.58;
        ConsPara.Field_r(ZZ>0.58)=ConsPara.Field_r(ZZ>0.58)*100;
        ConsPara.Field_l(ZZ>0.58)=ConsPara.Field_l(ZZ>0.58)*100;
    end



if strcmpi(WirePara.fold_condition,'8-fold')
    ConsPara.Num=[[1,1]*length(ConsPara.Field_l),[1,1]*length(ConsPara.Gradient_l),[1,1]*length(ConsPara.A_l)];
elseif strcmpi(WirePara.fold_condition,'4-fold')
    ConsPara.Num=[[1,1]*length(ConsPara.Field_l),[1,1]*length(ConsPara.Gradient_l),[1,1],[1,1]*length(ConsPara.A_l)];
end

    [id, coeff_wire,Mat_bsp_bas, Mat_dbsp_bas,gauss_weight, l_inteval, r_inteval, ind_eq,id2,~,~] = NURBS2constrain(NURBS_curve1, Wire_Id_ctrl,'remove');
    ind_eq{1,1}=[];



    %%

    A=[];
    b=[];
    tic
    
    
    %
    Aeq_cell=cell(1,1);
    beq_cell=cell(1,1);
    l_tmp=length(ind_eq{2,1});
    id1_tmp=[1:l_tmp,1:l_tmp];
    id2_tmp=length(coeff_wire)+[ind_eq{2,1},ind_eq{2,1}+1];
    val_tmp=[ones(l_tmp,1),-ones(l_tmp,1)];
    Aeq_cell{1,1}=sparse(id1_tmp,id2_tmp,val_tmp,l_tmp,3*length(coeff_wire));
    beq_cell{1,1}=zeros(l_tmp,1);

    l_tmp=length(ind_eq{3,1});
    id1_tmp=[1:l_tmp,1:l_tmp];
    id2_tmp=[ind_eq{3,1},ind_eq{4,1}];
    val_tmp=[ones(l_tmp,1),-ones(l_tmp,1)];



    Aeq=cell2mat(Aeq_cell);
    beq=cell2mat(beq_cell);



%%
    n_divide=10;
    OptPara = gene_OptPara(NURBS_curve1,WirePara,'Sigmoid of Distance',n_divide);  

    OptPara.gauss_weight_d=gauss_weight_d;
    OptPara.Mat_bsp_bas_d=Mat_bsp_bas_d;
    OptPara.Mat_dbsp_bas_d=Mat_dbsp_bas_d;
    OptPara.coeff_wire0=coeff_wire;
    
     options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',20000,'MaxIterations',2000,'StepTolerance',1e-10 ...
        ,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'UseParallel',true, ...
         'HessianFcn',@(coeff,lambda)HessianFun(coeff,lambda,WirePara_d,OptPara,ConsPara));
%% 
[f,gf]=Objective_phase2(coeff_wire,OptPara,WirePara,1);

    [coeff_opt,object_value]=fmincon(@(x)Objective_phase2(x,OptPara,WirePara_d,1),coeff_wire,A,b,Aeq,beq,l_inteval,r_inteval,...
            @(x1)Constrain_fieldv03(x1,WirePara_d,OptPara,ConsPara),options);
    toc




%%

coeff_opt=glo_coeff{53};% 340
NURBS_curve_opt= coeff2NURBS(coeff_opt, NURBS_curve_ini);
    save('coeff_iter.mat','coeff_opt','glo_coeff','NURBS_curve_ini','WirePara','NURBS_curve_opt');


Cell_wire_opt = nurbs_fold2all(NURBS_curve_opt,WirePara.size_coil,WirePara.fold_condition);


Cell_wire_eddy= nurbs_fold2all(NURBS_curve_opt,WirePara.size_coil,WirePara.fold_condition,5);

    DSVPara.CoilUnit=1;
    DSVPara.Size=[0.4;0.45;0.5];
    ShieldPara.Region=[0.45,0.7,0.45,-0.7];
    ShieldPara.SurfaceShape='Cylindrical';
spec_opt=Coil_performance(Cell_wire_opt,DSVPara,ShieldPara);

if 1
    near_wire(1).Wire_arrange=num2cell([1:18;[2:18,-1];0:17].');
    near_wire(2).Wire_arrange=num2cell([19:28;[20:28,-1];[0,19:27]].');

    FormPara.XiForm='uniform';

    NURBS_curve_opt=find_mindistance(NURBS_curve_opt,WirePara,near_wire,10,'off');
    distance_opt=cell(0,0);
    for cl1=1:length(NURBS_curve_opt)
        for cl2=1:length(NURBS_curve_opt(cl1).dis)
        dis_tmp=min(NURBS_curve_opt(cl1).dis{cl2,6});
        distance_opt{cl1,1}(cl2,1)=dis_tmp;
        end
    end
end



    Cell_wire_ini = nurbs_fold2all(NURBS_curve_ini,WirePara.size_coil,WirePara.fold_condition,10);

if 1

    for c1=1:size(Cell_wire_opt,1)
    Cell_wire_opt{c1,2}(end,:)=[];
    end
    for c1=1:size(Cell_wire_ini,1)
    Cell_wire_ini{c1,2}(end,:)=[];
    end
    WirePara.size_sect=[4,4;4,8]*1e-3;
    FastHenry_import(Cell_wire_ini,'standard_60cm_ini',WirePara);
    WirePara.size_sect=[4,4;4,8]*1e-3;
    FastHenry_import(Cell_wire_opt,'standard_60cm_step0',WirePara);
end

NURBS_surf=NURBS_curve_opt;



if 1
    tic
    size_c=WirePara.size_coil;
        WirePara1=WirePara;
    WirePara1.size_coil(:,[2,4])=WirePara1.size_coil(:,[2,4])*1.3;
    NURBS_surf=NURBS_curve_opt;
    NURBS_surf(1).KnotChange([2])=1;

    NURBS_surf=nurbs_equal_seg(NURBS_surf,WirePara1);

figure
    for c1=1:length(NURBS_surf)
         subplot(1,length(NURBS_surf),c1)
         xlim([0,pi/2]);
        NURBS_plot(NURBS_curve_ini(c1),{'--'},'no control point');
        Coil2D_ini(c1)=NURBS_plot(NURBS_surf(c1),{'-'},'control point');
    end

    save('Coil2D_final.mat','Coil2D_ini','size_c');
else
    figure
    for c1=1:length(NURBS_surf)
         subplot(1,length(NURBS_surf),c1)
         xlim([0,pi/2]);
        NURBS_plot(NURBS_curve_ini(c1),{'--'},'no control point');
        NURBS_plot(NURBS_surf(c1),{'-'},'control point');
    end
end


