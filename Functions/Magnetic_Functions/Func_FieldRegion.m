function R = Func_FieldRegion(Para)
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

if isfield(Para,'SegNum')
    Num_tu=Para.SegNum(1);
    Num_tv=Para.SegNum(2);
else
    Num_tu=32;
    Num_tv=32;
end

if isfield(Para,'Region')
    if length(Para.Region)==4
    region=[Para.Region(1),Para.Region(1),Para.Region(2),Para.Region(3),Para.Region(3),Para.Region(4)];
    elseif length(Para.Region)==2
    region=[Para.Region(1),Para.Region(1),Para.Region(2),Para.Region(1),Para.Region(1),-Para.Region(2)];
    elseif length(Para.Region)==3
        region=[Para.Region(1),Para.Region(2),Para.Region(3),Para.Region(1),Para.Region(2),-Para.Region(3)];
    elseif length(Para.Region)==6
        region=Para.Region;
    end
else
    error("缺少成像区域尺寸参数!");
end

if isfield(Para,'DeformPara')
    p=Para.DeformPara;
    if length(p)==1
        p(2)=2;
    elseif length(p)==2
    
    else
        error('区域变形参数错误！');
    end
    if any(p<=0)
    error('power p must be positive!');
    end
else
    p=[2,2];
end
if isfield(Para,'MirrorSymmetry')
    mirror_symmetry=Para.MirrorSymmetry;
else
    mirror_symmetry=false;
end
if isfield(Para,'RotateSymmetry')
    rotate_symmetry=Para.RotateSymmetry;
else
    rotate_symmetry='none';
end
if ~isfield(Para,'SurfaceShape')
    error('请指定区域形状');
end
               
     if strcmpi(Para.SurfaceShape,'spherical')
            if mirror_symmetry
                v=(1:(Num_tv/2))*pi/(Num_tv);
            else
                v=(1:(Num_tv-1))*pi/(Num_tv);
            end

            if strcmpi(rotate_symmetry,'4-fold') % cylindrical symmetry
                u=(0:(Num_tu/4))*2*pi/(Num_tu);
            elseif strcmpi(rotate_symmetry,'2-fold')
                u=(0:(Num_tu/2))*2*pi/(Num_tu);
            else
                u=(0:Num_tu-1)*pi*2/Num_tu;
            end

        [uu,vv]=ndgrid(u,v);
        
        R(:,1)=(region(1)+region(4))/2*sign(sin(vv(:)).*cos(uu(:))).*abs(sin(vv(:))).^(2/p(2)).*abs(cos(uu(:))).^(2/p(1));
        R(:,2)=(region(2)+region(5))/2*sign(sin(vv(:)).*sin(uu(:))).*abs(sin(vv(:))).^(2/p(2)).*abs(sin(uu(:))).^(2/p(1));
        R(:,3)=(region(3)-region(6))/2*sign(cos(vv(:))).*abs(cos(vv(:))).^(2/p(2))+(region(6)+region(3))/2;
        
        if mirror_symmetry
                R=[R;0,0,region(3)/2];
        else
             R=[R;0,0,region(3)/2;0,0,region(6)/2];
        end
       
    elseif strcmpi(Para.SurfaceShape, 'cylindrical')

            if strcmpi(rotate_symmetry,'4-fold') % cylindrical symmetry
                u=(0:(Num_tu/4))*2*pi/(Num_tu);
            elseif strcmpi(rotate_symmetry,'2-fold')
                u=(0:(Num_tu/2))*2*pi/(Num_tu);
            else
                u=(0:Num_tu-1)*pi*2/Num_tu;
            end

            if mirror_symmetry
                v=linspace(0,1,Num_tv);
                v=v(1:end/2);
            else
                v=linspace(0,1,Num_tv);
            end

            z=v*(region(6)-region(3))+region(3);
            r_y=v*(region(5)-region(2))+region(2);
            r_x=v*(region(4)-region(1))+region(1);

        [uu,zz]=ndgrid(u,z);
        [~,rr_x]=ndgrid(u,r_x);
        [~,rr_y]=ndgrid(u,r_y);
        R(:,1)=rr_x(:).*sign(cos(uu(:))).*abs(cos(uu(:))).^(2/p(1));
        R(:,2)=rr_y(:).*sign(sin(uu(:))).*abs(sin(uu(:))).^(2/p(1));
        R(:,3)=zz(:);
    elseif strcmpi(Para.SurfaceShape, 'plane')
        R = gene_plane(rotate_symmetry, Num_tu, region(1:3));

    elseif strcmpi(Para.SurfaceShape, 'biplane')
        R = gene_plane(rotate_symmetry, Num_tu, region(1:3));
        R=[R;gene_plane(rotate_symmetry, Num_tu, region(4:6));];
    elseif strcmpi(Para.SurfaceShape, 'line')
    
    elseif strcmpi(Para.SurfaceShape, 'circle')
    
            if strcmpi(rotate_symmetry,'4-fold') % cylindrical symmetry
                u=(0:(Num_tu/4))*2*pi/(Num_tu);
            elseif strcmpi(rotate_symmetry,'2-fold')
                u=(0:(Num_tu/2))*2*pi/(Num_tu);
            else
                u=(0:Num_tu-1)*pi*2/Num_tu;
            end
        [uu,vv]=ndgrid(u,region(3));

        R(:,1)=region(1)*sign(cos(uu(:))).*abs(cos(uu(:))).^(2/p(1));
        R(:,2)=region(2)*sign(sin(uu(:))).*abs(sin(uu(:))).^(2/p(1));
        R(:,3)=vv(:);
     else
         error('未指定成像区域形状!');
     end
end


function R = gene_plane(rotate_symmetry, Num_tu, region)
if strcmpi(rotate_symmetry,'4-fold') % cylindrical symmetry
    u1=linspace(0,1,Num_tu/2+1);
    u2=linspace(0,1,Num_tu/2+1);
elseif strcmpi(rotate_symmetry,'2-fold')
    u1=linspace(-1,1,Num_tu+1);
    u2=linspace(0,1,Num_tu/2+1);
else
    u1=linspace(-1,1,Num_tu+1);
    u2=linspace(-1,1,Num_tu+1);
end
[uu,vv,x3]=ndgrid(u1*region(1),u2*region(2),region(3));
R(:,1)=uu(:);
R(:,2)=vv(:);
R(:,3)=x3(:);
end