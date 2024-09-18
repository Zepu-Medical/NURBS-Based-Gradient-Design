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

load Results.mat    
Cell_wire_ini = nurbs_fold2all(NURBS_curve_ini,WirePara.size_coil,WirePara.fold_condition);
Cell_wire_step0 = nurbs_fold2all(NURBS_curve_step0,WirePara.size_coil,WirePara.fold_condition);
Cell_wire_step1 = nurbs_fold2all(NURBS_curve_step1,WirePara.size_coil,WirePara.fold_condition);
Cell_wire_step2 = nurbs_fold2all(NURBS_curve_step2,WirePara.size_coil,WirePara.fold_condition);
Cell_wire_step3 = nurbs_fold2all(NURBS_curve_step3,WirePara.size_coil,WirePara.fold_condition);
Cell_wire_step4 = nurbs_fold2all(NURBS_curve_step4,WirePara.size_coil,WirePara.fold_condition);
%% gradient efficiency and linearity
DSVPara.CoilUnit=1;
    DSVPara.Size=[0.26;0.3];
    ShieldPara.Region=[0.345,0.55,0.345,-0.55];
    ShieldPara.SurfaceShape='Cylindrical';
Coil_performance(Cell_wire_ini,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step0,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step1,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step2,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step3,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step4,DSVPara,ShieldPara);

%% 

figure 
subplot(3,2,1)
NURBS_plot_2D(NURBS_curve_ini,Coil2D_ini,WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(3,2,2)
NURBS_plot_2D(NURBS_curve_step0,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(3,2,3)
NURBS_plot_2D(NURBS_curve_step1,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(3,2,4)
NURBS_plot_2D(NURBS_curve_step2,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(3,2,5)
NURBS_plot_2D(NURBS_curve_step3,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(3,2,6)
NURBS_plot_2D(NURBS_curve_step4,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);

%% eddy field

Cell_wire_eddy=Cell_wire_ini;
eddy_current
Cell_wire_eddy=Cell_wire_step0;
eddy_current
Cell_wire_eddy=Cell_wire_step1;
eddy_current
Cell_wire_eddy=Cell_wire_step2;
eddy_current
Cell_wire_eddy=Cell_wire_step3;
eddy_current
Cell_wire_eddy=Cell_wire_step4;
eddy_current
%% distance 


    NURBS_curve_ini(4).WireId(7:9,1)=NURBS_curve_ini(4).WireId(6)+(1:3).';
    NURBS_curve_step0(4).WireId(7:9,1)=NURBS_curve_step0(4).WireId(6)+(1:3).';
    NURBS_curve_step1(4).WireId(7:9,1)=NURBS_curve_step1(4).WireId(6)+(1:3).';
    NURBS_curve_step2(4).WireId(7:9,1)=NURBS_curve_step2(4).WireId(6)+(1:3).';
    NURBS_curve_step3(4).WireId(7:9,1)=NURBS_curve_step2(4).WireId(6)+(1:3).';
    NURBS_curve_step4(4).WireId(7:9,1)=NURBS_curve_step2(4).WireId(6)+(1:3).';
    for cl1=1:length(NURBS_curve_ini)
        Id_tmp=NURBS_curve_ini(cl1).WireId;
        near_wire(cl1).Wire_arrange=num2cell([Id_tmp,[Id_tmp(2:end);-1],[0;Id_tmp(1:end-1)]]);
    end

    near_wire(3).Wire_arrange{3,2}=[4,34];
    near_wire(3).Wire_arrange{16,2}=[-1];
    near_wire(3).Wire_arrange{17,3}=[3,4];
    near_wire(4).Wire_arrange{1,3}=36;
    near_wire(4).Wire_arrange{5,2}=-1;
    near_wire(4).Wire_arrange{7,2}=-1;
    near_wire(4).Wire_arrange{7,3}=36;
    near_wire(4).Wire_arrange{6,2}=-1;
    near_wire(4).Wire_arrange{6,3}=[1,37];


    NURBS_curve_ini=find_mindistance(NURBS_curve_ini,WirePara,near_wire,20,'off');
    NURBS_curve_step0=find_mindistance(NURBS_curve_step0,WirePara,near_wire,20,'off');
    NURBS_curve_step1=find_mindistance(NURBS_curve_step1,WirePara,near_wire,20,'off');
    NURBS_curve_step2=find_mindistance(NURBS_curve_step2,WirePara,near_wire,20,'off');
    NURBS_curve_step3=find_mindistance(NURBS_curve_step3,WirePara,near_wire,20,'off');
 NURBS_curve_step4=find_mindistance(NURBS_curve_step4,WirePara,near_wire,20,'off');


    for cl1=1:length(NURBS_curve_ini)
        for cl2=1:length(NURBS_curve_ini(cl1).dis)

        distance_ini{cl1,1}(cl2,1)=min(NURBS_curve_ini(cl1).dis{cl2,6});
        distance_step0{cl1,1}(cl2,1)=min(NURBS_curve_step0(cl1).dis{cl2,6});
        distance_step1{cl1,1}(cl2,1)=min(NURBS_curve_step1(cl1).dis{cl2,6});
        distance_step2{cl1,1}(cl2,1)=min(NURBS_curve_step2(cl1).dis{cl2,6});
        distance_step3{cl1,1}(cl2,1)=min(NURBS_curve_step3(cl1).dis{cl2,6});
        distance_step4{cl1,1}(cl2,1)=min(NURBS_curve_step4(cl1).dis{cl2,6});
        end
    end

distance_all=[min(distance_ini{1,1}),min(distance_ini{3,1});
    min(distance_step0{1,1}),min(distance_step0{3,1});
    min(distance_step1{1,1}),min(distance_step1{3,1});
    min(distance_step2{1,1}),min(distance_step2{3,1});
    min(distance_step3{1,1}),min(distance_step3{3,1});  
    min(distance_step4{1,1}),min(distance_step4{3,1});
    ]
%%
[~,torque(1,:)] = Func_force(Cell_wire_ini,'Field_strength',[0,0,1]);
[~,torque(2,:)] = Func_force(Cell_wire_step0,'Field_strength',[0,0,1]);
[~,torque(3,:)] = Func_force(Cell_wire_step1,'Field_strength',[0,0,1]);
[~,torque(4,:)] = Func_force(Cell_wire_step2,'Field_strength',[0,0,1]);
[~,torque(5,:)] = Func_force(Cell_wire_step3,'Field_strength',[0,0,1]);
[~,torque(6,:)] = Func_force(Cell_wire_step4,'Field_strength',[0,0,1]);
%% fasthenry input

    for c1=1:size(Cell_wire_ini,1)
    Cell_wire_ini{c1,2}(end,:)=[];
    Cell_wire_step0{c1,2}(end,:)=[];
    Cell_wire_step1{c1,2}(end,:)=[];
    Cell_wire_step2{c1,2}(end,:)=[];
    Cell_wire_step3{c1,2}(end,:)=[];
    Cell_wire_step4{c1,2}(end,:)=[];
    end
    WirePara.size_sect=[5,5;5,5;5,10;5,5]*1e-3;
    FastHenry_import(Cell_wire_ini,'head_ini',WirePara);
    FastHenry_import(Cell_wire_step0,'head_step0',WirePara);
    FastHenry_import(Cell_wire_step1,'head_step1',WirePara);
    FastHenry_import(Cell_wire_step2,'head_step2',WirePara);
    FastHenry_import(Cell_wire_step3,'head_step3',WirePara);
    FastHenry_import(Cell_wire_step4,'head_step4',WirePara);
