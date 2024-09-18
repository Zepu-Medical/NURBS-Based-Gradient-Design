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
%% gradient efficiency and linearity
DSVPara.CoilUnit=1;
    DSVPara.Size=[0.45;0.5];
    ShieldPara.Region=[0.45,0.7,0.45,-0.7];
    ShieldPara.SurfaceShape='Cylindrical';
Coil_performance(Cell_wire_ini,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step0,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step1,DSVPara,ShieldPara);
Coil_performance(Cell_wire_step2,DSVPara,ShieldPara);


%% 

figure 
subplot(4,1,1)
NURBS_plot_2D(NURBS_curve_ini,Coil2D_ini,WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(4,1,2)
NURBS_plot_2D(NURBS_curve_step0,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(4,1,3)
NURBS_plot_2D(NURBS_curve_step1,[],WirePara,{'-','.'},'control point',0,1);  % 
xlabel('Integrated z-distance /m')
ylabel('$\varphi/\pi/2$',Interpreter='latex')
set(gca,'Linewidth',1,'FontName','Times New Roman','FontSize',12);
subplot(4,1,4)
NURBS_plot_2D(NURBS_curve_step2,[],WirePara,{'-','.'},'control point',0,1);  % 
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

%% distance 

    near_wire(1).Wire_arrange=num2cell([1:19;[2:19,-1];0:18].');
    near_wire(2).Wire_arrange=num2cell([1:9;[2:9,-1];0:8].');
    near_wire(3).Wire_arrange=num2cell([[1:9,20:22];[2:9,20:22,-1];[0,1:9,20:21]].');

    NURBS_curve_ini=find_mindistance(NURBS_curve_ini,WirePara,near_wire,20,'off');
    NURBS_curve_step0=find_mindistance(NURBS_curve_step0,WirePara,near_wire,20,'off');
    NURBS_curve_step1=find_mindistance(NURBS_curve_step1,WirePara,near_wire,20,'off');
    NURBS_curve_step2=find_mindistance(NURBS_curve_step2,WirePara,near_wire,20,'off');

    for cl1=1:length(NURBS_curve_ini)
        for cl2=1:length(NURBS_curve_ini(cl1).dis)

        distance_ini{cl1,1}(cl2,1)=min(NURBS_curve_ini(cl1).dis{cl2,6});
        distance_step0{cl1,1}(cl2,1)=min(NURBS_curve_step0(cl1).dis{cl2,6});
        distance_step1{cl1,1}(cl2,1)=min(NURBS_curve_step1(cl1).dis{cl2,6});
        distance_step2{cl1,1}(cl2,1)=min(NURBS_curve_step2(cl1).dis{cl2,6});
        end
    end
distance_all=[min(distance_ini{1,1}),min(distance_ini{3,1});
    min(distance_step0{1,1}),min(distance_step0{3,1});
    min(distance_step1{1,1}),min(distance_step1{3,1});
    min(distance_step2{1,1}),min(distance_step2{3,1});]

%% fasthenry input

    for c1=1:size(Cell_wire_ini,1)
    Cell_wire_ini{c1,2}(end,:)=[];
    Cell_wire_step0{c1,2}(end,:)=[];
    Cell_wire_step1{c1,2}(end,:)=[];
    Cell_wire_step2{c1,2}(end,:)=[];
    end
    WirePara.size_sect=[4,8;4,8;4,16]*1e-3;
    FastHenry_import(Cell_wire_ini,'sc_60cm_ini',WirePara);
    FastHenry_import(Cell_wire_step0,'sc_60cm_step0',WirePara);
    FastHenry_import(Cell_wire_step1,'sc_60cm_step1',WirePara);
    FastHenry_import(Cell_wire_step2,'sc_60cm_step2',WirePara);

