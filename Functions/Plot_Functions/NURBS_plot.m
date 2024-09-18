function [Coil2D,kappa_i] = NURBS_plot(NURBS_surf,Wiretype,type,Idtype)
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
dn=[16,16,16,16];
kappa_i={};
for c_surf=1:length(NURBS_surf)

    temp_wires=NURBS_surf(c_surf);
    for c_wire=1:length(temp_wires.Knot)
       if isempty(temp_wires.Xi{c_wire})
        % xi=linspace(0,temp_wires.Knot{c_wire}(end),32*temp_wires.Knot{c_wire}(end)+1);
        % xi=[0.05:temp_wires.Knot{c_wire}(end)-0.05].';

        dxi=0.5*log(linspace(0.5,1,dn(c_surf)/2+1))/log(0.5);
        dxi=[flip(dxi),1-dxi(1,2:end)]';
        xi=dxi(1:end-1,1)+[0:temp_wires.Knot{c_wire}(end)-1];
        xi=[reshape(xi,1,[]),temp_wires.Knot{c_wire}(end)];

            else
                 xi=temp_wires.Xi{c_wire};
            end
        
        [nurbs, dnurbs, d2nurbs] = Func_NURBSCurve(xi, temp_wires.Knot{c_wire},[temp_wires.ContrPoint{c_wire},temp_wires.Weight{c_wire}], 2,'Exponential');

        % dxi=0.001;
        % xi_m=xi-dxi;
        % xi_p=xi+dxi;
        % [~, dnurbs_m,~] = Func_NURBSCurve(xi_m, temp_wires.Knot{c_wire},[temp_wires.ContrPoint{c_wire},temp_wires.Weight{c_wire}], 2,'Exponential');
        % [~, dnurbs_p,~] = Func_NURBSCurve(xi_p, temp_wires.Knot{c_wire},[temp_wires.ContrPoint{c_wire},temp_wires.Weight{c_wire}], 2,'Exponential');
        % dnurbs_diff=(dnurbs_p-dnurbs_m)/2/dxi;
        if c_wire==2
            1;
        end
        
        dnurbs=dnurbs.*[0.4,0.8];
        d2nurbs=d2nurbs.*[0.4,0.8];
        kappa=abs(dnurbs(:,1).*d2nurbs(:,2)-dnurbs(:,2).*d2nurbs(:,1))./vecnorm(dnurbs,2,2).^3;
        
        kappa_i{end+1,1}=kappa;
        kappa_i{end,2}=max(kappa);
        %    nurbs=dnurbs;
        Coil2D(c_surf).Wire2D{c_wire,1}=nurbs;
        Coil2D(c_surf).WireId(c_wire,1)=NURBS_surf(c_surf).WireId(c_wire,1);
        Coil2D(c_surf).seg_idx{c_wire,1}=1:dn(c_surf):size(nurbs,1);
        if size(nurbs,1)-Coil2D(c_surf).seg_idx{c_wire,1}(end)>dn(c_surf)/2
            Coil2D(c_surf).seg_idx{c_wire,1}(end+1)=size(nurbs,1);
        elseif size(nurbs,1)-Coil2D(c_surf).seg_idx{c_wire,1}(end)>=0
            Coil2D(c_surf).seg_idx{c_wire,1}(end)=size(nurbs,1);
        end
        hold on
        if exist('Wiretype')

            plot(nurbs(:,1),nurbs(:,2),Wiretype{1},'LineWidth',1);  %原始
        else
            plot(nurbs(:,1),nurbs(:,2),'LineWidth',1);  %原始
        end
        if exist('type','var')
            if strcmpi(type,'control point')
                if length(Wiretype)>=2
                    scatter(temp_wires.ContrPoint{c_wire}(:,1),temp_wires.ContrPoint{c_wire}(:,2),Wiretype{2});
                else
                    scatter(temp_wires.ContrPoint{c_wire}(:,1),temp_wires.ContrPoint{c_wire}(:,2),'*');
                end
            end
        end
        if exist("Idtype")
            if Idtype
                text(nurbs(3,1),nurbs(3,2),num2str(NURBS_surf(c_surf).WireId(c_wire)));

            end
        end

    end
end

