function [Coil2D] = NURBS_plot_2D(NURBS_surf,Coil2D_ini,WirePara,Wiretype,type,Idtype,color_ind)
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
z0=0;
plot(z0*ones(100,1),linspace(0,1,100),'-.');
load("ColorSet.mat",'rgbset');
for c_surf=1:length(NURBS_surf)
    plot(0*ones(100,1),linspace(0,1,100),'-.',LineWidth=1,Color=[0,0,0]);
    d0=vecnorm(WirePara.size_coil(c_surf,3:4)-WirePara.size_coil(c_surf,1:2));
    temp_wires=NURBS_surf(c_surf);
    for c_wire=1:length(temp_wires.Knot)
        %if ~exist('temp_wires.Xi','var')
        xi=linspace(0,temp_wires.Knot{c_wire}(end),16*temp_wires.Knot{c_wire}(end)+1);
        %     else
        %         xi=temp_wires.Xi{c_wire};
        %     end

        [nurbs, dnurbs] = Func_NURBSCurve(xi, temp_wires.Knot{c_wire},[temp_wires.ContrPoint{c_wire},temp_wires.Weight{c_wire}], 2,'Exponential');
        %    nurbs=dnurbs;

        Coil2D(c_surf).Wire2D{c_wire,1}=nurbs;
        Coil2D(c_surf).WireId(c_wire,1)=NURBS_surf(c_surf).WireId(c_wire,1);
        Coil2D(c_surf).seg_idx{c_wire,1}=1:dn(c_surf):16*temp_wires.Knot{c_wire}(end)+1;
        hold on

            if exist('Wiretype')
                plot(nurbs(:,2)*d0+z0,nurbs(:,1)/pi*2,Wiretype{1},'LineWidth',1.5,'Color',rgbset(color_ind,:));  %原始
                if ~isempty(Coil2D_ini)
                    nurbs_ini=Coil2D_ini(c_surf).Wire2D{c_wire};
                    plot(nurbs_ini(:,2)*d0+z0,nurbs_ini(:,1)/pi*2,'--','LineWidth',1,'Color',rgbset(2,:));  %原始
                end
            else
                plot(nurbs(:,2)*d0+z0,nurbs(:,1)/pi*2,'LineWidth',1);  %原始
            end

        if exist('type','var')
            if strcmpi(type,'control point')
                if 1

                if length(Wiretype)>=2
                    scatter(temp_wires.ContrPoint{c_wire}(:,2)*d0+z0,temp_wires.ContrPoint{c_wire}(:,1)/pi*2,Wiretype{2},'LineWidth',0.5,'MarkerEdgeColor',rgbset(2,:));
                else
                    scatter(temp_wires.ContrPoint{c_wire}(:,2)*d0+z0,temp_wires.ContrPoint{c_wire}(:,1)/pi*2,'*');
                end

                end
            end
        end

        if exist("Idtype")
            if Idtype
                text(nurbs(3,1),nurbs(3,2),num2str(NURBS_surf(c_surf).WireId(c_wire)));
            end
        end

    end
        yticks(linspace(0,1,6));
                    z0=z0+d0;
        plot(z0*ones(100,1),linspace(0,1,100),'-.',LineWidth=1,Color=[0,0,0]);
        if c_surf<length(NURBS_surf)
        d=vecnorm(WirePara.size_coil(c_surf,3:4)-WirePara.size_coil(c_surf+1,1:2));
        if d<1e-6            
        else

            z0=0;
            figure
        end
        end

end

