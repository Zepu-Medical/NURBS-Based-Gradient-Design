function B_n=Coil_performance(Cell_wire,DSVPara,ShieldPara)
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


if isfield(DSVPara,'Size')
    DSVs=DSVPara.Size;
    if ~isempty(DSVs)
        if size(DSVs,1)==1&&size(DSVs,2)>3
            DSVs=DSVs.';
        end

        DSVr=[];
        if size(DSVs,2)==1
            DSVr=repmat(DSVs,[1,3]);
        elseif size(DSVs,2)==2
            DSVr(:,1)=DSVs(:,1);
            DSVr(:,2)=DSVs(:,1);
            DSVr(:,3)=DSVs(:,2);
        elseif size(DSVs,2)==3
            DSVr=DSVs;
        end
    else
    
        error('no DSV size！');
    end
else
    error('no DSV size！');
end
[harmonicorder,sinuorder]=Func_Wire2HarmonicOrder(Cell_wire,0.2);
if isfield(DSVPara,'CoilUnit')
    coilunit=DSVPara.CoilUnit;
    if coilunit>=1
    ind_unit=[4,2,3];
    ind=ind_unit(coilunit);
    elseif coilunit==0
    ind=1;
    end
else
    [~,ind]=max(abs(harmonicorder(1:4,3)));
    table_unit=[0,2,3,1];
    coilunit=table_unit(ind);
end
G0=(-1)^(ind+1)*harmonicorder(ind,3)*sqrt((2+1)*2/4/pi*(factorial(0)./factorial(2)));
if ~isfield(ShieldPara,'MirrorSymmetry')
    ShieldPara.MirrorSymmetry=0;
end

if ~isfield(ShieldPara,'RotateSymmetry')
    ShieldPara.RotateSymmetry='none';
end
if ~isfield(ShieldPara,'Region')
    ShieldPara.Region=max(DSVr)*[1.365,-1.778,1.365,1.778];
end
ShieldPara.SegNum=[32,64];
%% Linearity

Phi = (0:7.5:352.5)*pi/180;

% theta = (5.6:10.55:174.4)/180*pi;
theta = (0:5:180)/180*pi;






for c1=1:size(DSVr,1)
    DSV_linearity=DSVr(c1,:);
% DSV Points definition
NA = length(theta); 
Np = length(Phi); 
DSVposition(:,1) =DSV_linearity(1)/2*kron(sin(theta),cos(Phi)).';  %%
DSVposition(:,2) =DSV_linearity(2)/2*kron(sin(theta),sin(Phi)).';
DSVposition(:,3) =DSV_linearity(3)/2*kron(cos(theta),ones(size(Phi))).';

Cell_magneticfield_linearity=cell(1,3);
Cell_magneticfield_linearity{1,1}=[];
Cell_magneticfield_linearity{1,2}=DSVposition;
Cell_magneticfield_linearity{1,3}='Cartesian';
Cell_magneticfield_linearity=Func_MagneticField(Cell_magneticfield_linearity,Cell_wire,3);

Bz=Cell_magneticfield_linearity{1}(:,3);

GoM = DSVposition(:,coilunit)\Bz;
if exist('G0','var')
    GoM =sign(GoM)*abs(G0);
elseif c1==1
    GoM = DSVposition(:,coilunit)\Bz; %% 
end
  
TargetField = G0*DSVposition(:,coilunit);
% Gloabl Linearity
disp(['DSV ',num2str(DSV_linearity(1)),'m*',num2str(DSV_linearity(2)),'m*',num2str(DSV_linearity(3)),'m'])
disp(['Gradient efficiency ',num2str(G0*1e6),' microT/m/A'])

site_positive=find(DSVposition(:,coilunit)>0);
Bzp=Bz(site_positive);
TFp=TargetField(site_positive);
Positionp=DSVposition(site_positive,:);

MAXDBz=max((Bzp - TFp).*sign(TFp))/max(abs(TFp));
MINDBz=min((Bzp - TFp).*sign(TFp))/max(abs(TFp));

disp(['Maximum linearity Deviation = ',num2str(MAXDBz*100),',',num2str(MINDBz*100),' %']);

nlinear=find(abs(DSVposition(:,coilunit))<1e-5);
linearity=(Bz - TargetField)./TargetField;
linearity(nlinear)=[];
disp(['Maximum local Deviation = ',num2str(max(linearity)*100),',',num2str(min(linearity)*100),' %']);
end

%% Shielding 

    Cell_shield{1,2}=Func_FieldRegion(ShieldPara); 
    Cell_shield{1,1}=0*Cell_shield{1,2}(:,2); 
    Cell_shield{1,3}='Cartesian';
    fieldcomponent=[1,2,3];
   [Cell_shield] = Func_MagneticField(Cell_shield,Cell_wire,fieldcomponent);

   B_s=Cell_shield{1};
       phi=atan2(Cell_shield{1,2}(:,2),Cell_shield{1,2}(:,1));
       Br=reshape(sum(Cell_shield{1,1}(:,1:2).*[cos(phi),sin(phi)],2),[ShieldPara.SegNum]);
       B_n=fftshift(fft(Br,[],1),1);

   ratio=sum(sqrt(sum(Cell_shield{1,1}.^2,2)))/size(Cell_shield{1,2},1)/abs(GoM)*100;
   disp(['Shield Region ',num2str(ShieldPara.Region),'m']);
   disp(['Shield ratio ',num2str(ratio*100),'%']);
   ratio_max=max(sqrt(sum(Cell_shield{1,1}.^2,2)))/abs(GoM)*100;
   disp(['Max ratio ',num2str(ratio_max*100),'%']);
end

function potion=Func_maxlinear(Bz,DSVPnt,index,coilunit)
GoM=DSVPnt(index,coilunit)\Bz(index);
potion=max(abs(Bz(index)-GoM.*DSVPnt(index,coilunit)))./max(abs( GoM.*DSVPnt(index,coilunit)));
end