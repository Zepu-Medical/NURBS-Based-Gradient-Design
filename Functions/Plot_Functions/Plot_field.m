function [G0,harmonicorder,sinuorder,Cell_magneticfield]=Plot_field(Cell_wire,PlotPara)
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
if isfield(PlotPara,'DeltaX')
deltax=PlotPara.DeltaX;
else
    deltax=0.01;
end
if isfield(PlotPara,'DSVSize')
    DSVsize=PlotPara.DSVSize;
else
    error('未指定DSV尺寸，请在参数中定义DSV直径：DSVSize=[X,Y,Z]');
end
[harmonicorder,sinuorder]=Func_Wire2HarmonicOrder(Cell_wire,0.2);
if isfield(PlotPara,'CoilUnit')
    coilunit=PlotPara.CoilUnit;
    table_unit=[0,2,3,1];
    if coilunit<1
    ind==1;
    else
    ind_table=[4,2,3];
    ind=ind_table(coilunit);
    end
else
    [~,ind]=max(abs(harmonicorder(1:4,3)));
    table_unit=[0,2,3,1];
    coilunit=table_unit(ind);
end
G0=(-1)^(ind+1)*harmonicorder(ind,3)*sqrt((2+1)*2/4/pi*(factorial(0)./factorial(2)));

if isfield(PlotPara,'PlotRange')
    DSV_r=cell(0,3);
    for c1=1:size(PlotPara.PlotRange,1)
            plot_region{1,1}=PlotPara.PlotRange{c1,1}(1):deltax:PlotPara.PlotRange{c1,1}(end);
            plot_region{2,1}=PlotPara.PlotRange{c1,2}(1):deltax:PlotPara.PlotRange{c1,2}(end);
            plot_region{3,1}=PlotPara.PlotRange{c1,3}(1):deltax:PlotPara.PlotRange{c1,3}(end);
            tmp_DSV={plot_region{1},plot_region{2},plot_region{3}};
    DSV_r=[DSV_r;tmp_DSV];
    
    end
else
    plot_region{1,1}=round(-DSVsize(1)*150)/100:deltax:round(DSVsize(1)*150)/100;
    plot_region{2,1}=round(-DSVsize(2)*150)/100:deltax:round(DSVsize(2)*150)/100;
    plot_region{3,1}=round(-DSVsize(3)*150)/100:deltax:round(DSVsize(3)*150)/100;
    DSV_r={plot_region{1},plot_region{2},plot_region{3};
    plot_region{1},plot_region{2},plot_region{3};};
    if coilunit==1
        DSV_r{1,3}=0;
        DSV_r{2,2}=0;
    elseif coilunit==2
        DSV_r{1,3}=0;
        DSV_r{2,1}=0;      
    elseif coilunit==3
        DSV_r{1,2}=0;
        DSV_r(2,:)=[];
    end
end


% DSV_max=0.4;

lab={'x/m','y/m','z/m'};
labf={'Gx','Gy','Gz',};
labr={'x','y','z'};

for c1=1:size(DSV_r,1)
[XX,YY,ZZ]=ndgrid(DSV_r{c1,1},DSV_r{c1,2},DSV_r{c1,3});
Cell_magneticfield{c1,2}=[XX(1:end).',YY(1:end).',ZZ(1:end).'];
Cell_magneticfield{c1,3}='Cartesian'; %'cylindrical','spherical'
end
Cell_magneticfield{c1+1,2}=[0,0,0];
Cell_magneticfield{c1+1,3}='Cartesian'; %'cylindrical','spherical'
[Cell_magneticfield] = Func_MagneticField(Cell_magneticfield,Cell_wire,[3]);
[Cell_gradientfield] = Func_GradientField(Cell_magneticfield,Cell_wire,coilunit);


if coilunit==1
    ind=4;
elseif coilunit==2
    ind=2;
elseif coilunit==3
    ind=3;
end


for c0=1:size(DSV_r,1)
DSV_in=DSV_r(c0,:);

c2=1;
for c1=1:length(DSV_in)
    l(c1)=length(DSV_in{c1});
    
end
    ind_p=find(l==1);
    if length(ind_p)~=1
        error('磁场区域维度错误！');
    end
    ind_rem=1:3;
    ind_rem(ind_p)=[];


Bz=reshape(Cell_magneticfield{c0,1}(:,3),[length(DSV_in{ind_rem(1)}),length(DSV_in{ind_rem(2)})]); 

Gfield=reshape(Cell_gradientfield{c0,1}(:,coilunit),[length(DSV_in{ind_rem(1)}),length(DSV_in{ind_rem(2)})]); 

[UU,VV]=ndgrid(DSV_in{ind_rem(1)},DSV_in{ind_rem(2)});

%Bz(abs(Bz)>1e-5)=0;

figure('Position',[100,100,1400,600]);
subplot(1,3,1)
colormap(jet)
contour_field=[min(Bz(:)),G0*[-0.55:0.05:0.55]];
wf_field=contourf(UU,VV,Bz,contour_field);%shading interp
clim(sort([contour_field(2),contour_field(end)]));

xlabel(lab{ind_rem(1)});
ylabel(lab{ind_rem(2)});
title('B_z');

subplot(1,3,2)
colormap(jet)
contour_gra=[min(Gfield(:))/G0,(1+[-1.55:0.05:-0.05,-0.03,-0.01,0.01,0.03,0.05:0.05:0.55])];

contourf(UU,VV,Gfield/G0,contour_gra,'--','ShowText','on');

clim(sort([contour_gra(2),contour_gra(end)]));
xlabel(lab{ind_rem(1)});
ylabel(lab{ind_rem(2)});

title([labf{coilunit},'=',num2str(G0),', ',labr{ind_p},'=',num2str(DSV_in{ind_p}),'m']);

Bz0=Cell_magneticfield{end,1}(1,3);
subplot(1,3,3)
if ind_rem(1)==coilunit
ideal_field=G0*UU;
elseif ind_rem(2)==coilunit
ideal_field=G0*VV;
end

contour_relative=[-16,-12,-8,-4,-1,1,4,8,12,16];
wf_gradient=contour(UU,VV,(Bz-ideal_field)/(G0*DSVsize(coilunit)/2)*100,contour_relative,'--','ShowText','on');
clim([-0.16,0.16]);
xlabel(lab{ind_rem(1)});
ylabel(lab{ind_rem(2)});
title('Relative Linearity Deviation(%)');
hold on
if DSV_in{ind_p}<DSVsize(ind_p)/2
    phi_tmp=linspace(0,2*pi,65);
    factor_tmp=sqrt(1-4*DSV_in{ind_p}^2/DSVsize(ind_p)^2);
plot(factor_tmp*DSVsize(ind_rem(1))/2*cos(phi_tmp),factor_tmp*DSVsize(ind_rem(2))/2*sin(phi_tmp),'--');
end

end