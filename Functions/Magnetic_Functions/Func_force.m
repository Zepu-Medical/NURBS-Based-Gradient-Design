function [force,torque] = Func_force(Cell_wire,MagnetField_Type,Magnet_info,Center)
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
 
if ~exist('Center','var')
    Center=[0,0,0];
end


Vari_wiresegmentposition=[];
Vari_wiresegment=[];

tic;
for Count1=1:size(Cell_wire,1)
    %
    if size(Cell_wire{Count1,2},1)>1
        switch Cell_wire{Count1,3}
            case 'Cylindrical'
                Temp_wireposition=[Cell_wire{Count1,2}(:,1).*cos(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,1).*sin(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,3)];
            case 'Spherical' %
                Temp_wireposition=repmat(Cell_wire{Count1,2}(:,2),[1,3]).*[sin(Cell_wire{Count1,2}(:,3)).*cos(Cell_wire{Count1,2}(:,2))...
                    ,sin(Cell_wire{Count1,2}(:,3)).*sin(Cell_wire{Count1,2}(:,2)),cos(Cell_wire{Count1,2}(:,3))];
            otherwise
                Temp_wireposition=Cell_wire{Count1,2};
        end


        Temp_wiresegmentposition=(Temp_wireposition(2:end,:)+Temp_wireposition(1:end-1,:))/2;%
        Temp_wiresegment=Cell_wire{Count1,1}*(Temp_wireposition(2:end,:)-Temp_wireposition(1:end-1,:));%
        Vari_wiresegmentposition=[Vari_wiresegmentposition;Temp_wiresegmentposition];%
        Vari_wiresegment=[Vari_wiresegment;Temp_wiresegment];%
    end
end
    Num_seg=size(Vari_wiresegment,1);
%% 
if strcmp(MagnetField_Type,'Magnet_Size')
    B0=Func_Magnet_field(Vari_wiresegmentposition,Magnet_info);
elseif strcmp(MagnetField_Type,'SH_coeffs')
    B0=Func_HarmonicOrder2Field(Vari_wiresegmentposition,Magnet_info);
elseif strcmp(MagnetField_Type,'Field_strength')
    B0=repmat(Magnet_info(:).',[Num_seg,1]);
end

F_seg=cross(Vari_wiresegment,B0,2);
Torque_seg=cross(Vari_wiresegmentposition-Center,F_seg,2);
force=sum(F_seg,1);
torque=sum(Torque_seg,1);
disp(['Force is ','F(N/A)=',num2str(force)]);
disp(['Torque is ','M(Nm/A)=',num2str(torque)]);
end



