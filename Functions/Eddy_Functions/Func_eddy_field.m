function Bz = Func_eddy_field(l_z, r_n, J_phi,rt,zt)
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

phi_t=(0:7)/8*2*pi;
z_t=[-zt,0,zt];
[PPhi_t,ZZ_t]=ndgrid(phi_t,z_t);
P_targ=[rt*ones(24,1),PPhi_t(:),ZZ_t(:)];
P_targ=[P_targ;0,0,rt;0,0,-rt];
r_t=rt;
z=((0:0.02:1)-0.5)*l_z;
phi=(1:1:100)/100*2*pi;
[ZZ,PPhi]=ndgrid(z,phi);

P_curr=[ones(numel(ZZ),1),PPhi(:),ZZ(:)];
Bz=zeros(size(P_targ,1),1);

for c1=1:length(r_n)
    r1=r_n(c1);
    
    Dis_si=((P_targ(:,3)-P_curr(:,3).').^2+(P_targ(:,1).^2+r1^2-2*r1*P_targ(:,1).*cos(P_targ(:,2)-P_curr(:,2).'))).^(-3/2);
    
    sensitive=Dis_si.*(r1-P_targ(:,1).*cos(P_targ(:,2)-P_curr(:,2).'));
    J_phi_tmp=J_phi(:,:,c1);
    Bz=Bz+r1*(z(2)-z(1))*(phi(2)-phi(1))*(sensitive*J_phi_tmp(:));

    
end
Bz=Bz*1e-7;
end