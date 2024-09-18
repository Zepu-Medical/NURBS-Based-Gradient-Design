function Cell_wire = wire_fold2all(Coil2D_ini,size_coil,FoldCondition)
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

Cell_wire_fold=cell(0,1);
Wire_Id=zeros(0,1);
Cell_wire_id_temp=cell(0,1);
surf_size=cell(0,1);
surf_id=cell(0,1);
for c1=1:length(Coil2D_ini)

    for c2=1:length(Coil2D_ini(c1).Wire2D)


    nurbs_fold=Coil2D_ini(c1).Wire2D{c2,1};
   
    temp_wire=[size_coil(c1,1)*(1-nurbs_fold(:,2))+size_coil(c1,3)*nurbs_fold(:,2),...
        nurbs_fold.*[1,size_coil(c1,4)-size_coil(c1,2)]+[0,size_coil(c1,2)]];  %。
    
    surf_size{end+1,1}=repmat(size_coil(c1,:),[size(temp_wire,1)-1,1]); %
    surf_id{end+1,1}=c1*ones(size(temp_wire,1)-1,1);
    Cell_wire_fold{end+1,1}=temp_wire;
    Wire_Id(end+1,1)=Coil2D_ini(c1).WireId(c2);

    end
end


for c1=min(Wire_Id):max(Wire_Id)
    ind_tmp=find(Wire_Id==c1);
    temp_wire=[];
    temp_surf_size=[];
    temp_surf_id=[];
    if length(ind_tmp)==1
        Cell_wire_temp{c1,1}=Cell_wire_fold{ind_tmp};
        Cell_wire_surf_temp{c1,1}=surf_size{ind_tmp};
        Cell_wire_id_temp{c1,1}=surf_id{ind_tmp};
    elseif length(ind_tmp)>1
        c4=0;
        while(~isempty(ind_tmp))
            if isempty(temp_wire)
            temp_wire=Cell_wire_fold{ind_tmp(1),1};
            temp_surf_size=surf_size{ind_tmp(1),1};
            temp_surf_id=surf_id{ind_tmp(1),1};
            ind_tmp(1)=[];
            end
            l0=length(ind_tmp);
            d0=1e-2;
            while(length(ind_tmp)>=l0)
                for c2=1:length(ind_tmp)
                    v1=temp_wire(1,:);
                    v2=Cell_wire_fold{ind_tmp(c2),1}(end,:);
                    
                    if sqrt(v1(1)^2+v2(1)^2-2*v1(1)*v2(1)*cos(v1(2)-v2(2))+(v1(3)-v2(3))^2)<d0

                        temp_wire=[Cell_wire_fold{ind_tmp(c2),1};temp_wire(2:end,:)];
                        temp_surf_size=[temp_surf_size;surf_size{ind_tmp(c2),1}];
                        temp_surf_id=[surf_id{ind_tmp(c2),1};temp_surf_id];
                        ind_tmp(c2)=[];
                        break;
                    end
                end
            
                for c2=1:length(ind_tmp)

                    v1=temp_wire(end,:);
                    v2=Cell_wire_fold{ind_tmp(c2),1}(1,:);
                    if sqrt(v1(1)^2+v2(1)^2-2*v1(1)*v2(1)*cos(v1(2)-v2(2))+(v1(3)-v2(3))^2)<d0

                        temp_wire=[temp_wire;Cell_wire_fold{ind_tmp(c2),1}(2:end,:)];
                        temp_surf_size=[temp_surf_size;surf_size{ind_tmp(c2),1}];
                        temp_surf_id=[temp_surf_id;surf_id{ind_tmp(c2),1}];
                        ind_tmp(c2)=[];
                        break;
                    end
                end
                d0=d0*2;
            end   
            c4=c4+1; 
            if c4>4
                error('~');
            end
        end
        Cell_wire_temp{c1,1}=temp_wire;
        Cell_wire_surf_temp{c1,1}=temp_surf_size;
        Cell_wire_id_temp{c1,1}=temp_surf_id;
    end

end


Cell_wire=cell(0,3);
for c1=1:size(Cell_wire_temp,1)

    
        nurbs_fold=Cell_wire_temp{c1,1};
        nurbs_id=Cell_wire_id_temp{c1,1};
        nurbs_size=Cell_wire_surf_temp{c1,1};
        nurbs_fold=[nurbs_fold;flip(nurbs_fold(1:end-1,:),1).*[1,-1,1];];
        nurbs_id=[nurbs_id;flip(nurbs_id(1:end,:),1)];
        nurbs_size=[nurbs_size;flip(nurbs_size(:,:),1)];
        Cell_wire_temp1{1,2}=nurbs_fold;
        Cell_wire_temp1{2,2}=[0,pi,0]+[1,-1,1].*nurbs_fold;
        Cell_wire_temp1{1,4}=nurbs_id;
        Cell_wire_temp1{2,4}=nurbs_id;
        Cell_wire_temp1{1,5}=nurbs_size;
        Cell_wire_temp1{2,5}=nurbs_size;
if strcmpi(FoldCondition,'4-fold')

elseif strcmpi(FoldCondition,'8-fold')
    
    for c3=1:2
        Cell_wire_temp1{c3+2,2}=Cell_wire_temp1{c3,2}.*[1,1,-1];
        Cell_wire_temp1{c3+2,4}=Cell_wire_temp1{c3,4};
        Cell_wire_temp1{c3+2,5}=Cell_wire_temp1{c3,5}.*[1,-1,1,-1];
    end
end


    for c3=1:size(Cell_wire_temp1,1)
    Cell_wire_temp1{c3,1}=1;    
    Cell_wire_temp1{c3,3}='Cylindrical';
    end
    Cell_wire=[Cell_wire;Cell_wire_temp1];
end

end