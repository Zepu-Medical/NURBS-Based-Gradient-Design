function NURBS_surf=nurbs_equal_seg(NURBS_surf,WirePara)
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

for c_surf=1:length(NURBS_surf)

    temp_wires=NURBS_surf(c_surf);
    sz=[WirePara.size_coil(c_surf,1)/2+WirePara.size_coil(c_surf,3)/2,abs(WirePara.size_coil(c_surf,2)-WirePara.size_coil(c_surf,4))];
    for c_wire=1:length(temp_wires.Knot)
    
    coeff=[temp_wires.ContrPoint{c_wire},temp_wires.Weight{c_wire}].*[sz,1];
    
    Knot=temp_wires.Knot{c_wire};
    xi_seg=0:0.5:Knot(end);
    total_length=0;
    d_length=zeros(length(xi_seg),1);
    
   % parfor c1=1:length(xi_seg)-1
        d_length(2:end,1)=nurbs_length(xi_seg(2:end),xi_seg(1:end-1),Knot,coeff);
   % end
    
    total_length=cumsum(d_length);
    l_xi=temp_wires.Knot{c_wire}(end)*16;

    if isfield(temp_wires,'KnotChange')
        if length(temp_wires.KnotChange)>=c_wire
            l_xi=l_xi+temp_wires.KnotChange(c_wire)*16;
        end
    end
    seg_length=(1:l_xi-1)/l_xi*total_length(end);
        dl=seg_length-total_length(1:end-1);
        dr=total_length(2:end)-seg_length;
        ind1=dl>0&dr>=0;
    xi_0=zeros(length(seg_length),1);

    for c1=1:length(xi_seg)-1
        id=find(ind1(c1,:));
        xi_0_tmp=linspace(xi_seg(c1),xi_seg(c1+1),length(id)+2);
        xi_0(id)=xi_0_tmp(1,2:end-1);
    end

    xi_0=[0;xi_0];

    
    parfor c1=1:length(xi_0)-1
        d_length(c1,1)=nurbs_length(xi_0(c1+1),xi_0(c1),Knot,coeff);
    end
    
    xi_0(1)=[];

    total_length_xi_0=cumsum(d_length);
    seg_length=seg_length(:)-total_length_xi_0;
    x_val=zeros(length(seg_length),1);

    %for c1=1:l_xi-1


         options=optimoptions('lsqnonlin','Display','none','SpecifyObjectiveGradient',true);%,HessianFcn=@(xi,lambda)Hessen_nurbs_length(xi,lambda,xi_0,Knot,coeff,seg_length));

        x_val=lsqnonlin(@(xi_i)nurbs_seg_obj(xi_i,xi_0,Knot,coeff,seg_length),xi_0,zeros(l_xi-1,1),max(Knot)*ones(l_xi-1,1),[],[],[],[],[],options);

    x_val=[0;x_val;Knot(end)];

    d_length=nurbs_length(x_val(2:end),x_val(1:end-1),Knot,coeff);
    
    temp_wires.Xi{c_wire}=x_val;
    
    end

    NURBS_surf(c_surf)=temp_wires;
    toc
end
        if isfield(NURBS_surf,'knot_change')
            NURBS_surf=rmfield(NURBS_surf,'KnotChange');
        end