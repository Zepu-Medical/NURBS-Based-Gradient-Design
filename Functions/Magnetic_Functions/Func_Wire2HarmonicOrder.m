function [harmonicorder,sinuorder,Cell_magneticfield]=Func_Wire2HarmonicOrder(Cell_wire,R)
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


Num_phi=24;
Num_theta=36;
Order1=10;
DSV={R,2*pi*(0:Num_phi-1)/Num_phi,pi*(1:2:2*Num_theta-1)/(2*Num_theta)}; %% 
[R_3,Phi_3,Theta_3]=ndgrid(DSV{1},DSV{2},DSV{3});

Cell_magneticfield{1,2}=[R_3(1:end).',Phi_3(1:end).',Theta_3(1:end).'];
Cell_magneticfield{1,3}='Spherical'; %'cylindrical','spherical'
fieldcomponent=[3];
[Cell_magneticfield] = Func_MagneticField(Cell_magneticfield,Cell_wire,fieldcomponent);



magnetic_field=reshape(Cell_magneticfield{1}(:,3),[length(DSV{2}),length(DSV{3})]);

harmonicorder=[];
sinuorder=cell(1,2);
for l=0:(Order1)
    Temp_harmonicorder=zeros(2*l+1,3);
    normfactor=DSV{1}^(-l)*sqrt((2*l+1)/4/pi*(factorial(l-abs(-l:l))./factorial(l+abs(-l:l)))).';
    legendre_theta=legendre(l,cos(DSV{3}));
    legendre_theta=[legendre_theta(end:-1:2,:);legendre_theta].*repmat(sin(DSV{3}),[2*l+1,1]);
    Temp_field_m_theta=fft(magnetic_field,[],1);
    Temp_harmonicorder(:,3)=2*pi/Num_phi*(DSV{3}(2)-DSV{3}(1))*normfactor.*sum(legendre_theta.*[Temp_field_m_theta(end-l+1:end,:);Temp_field_m_theta(1:l+1,:)],2);
    Temp_harmonicorder(:,3)=exp2sinu(l)*Temp_harmonicorder(:,3); 
    Temp_harmonicorder(:,1)=l;
    Temp_harmonicorder(:,2)=-l:l;
    harmonicorder=[harmonicorder;Temp_harmonicorder];
    sinuorder{1,1}=[sinuorder{1,1};Temp_harmonicorder(l+1:end,:)];
    sinuorder{1,2}=[sinuorder{1,2};[l,0,0]];
    if l>0
    sinuorder{1,2}=[sinuorder{1,2};Temp_harmonicorder(l:-1:1,:).*[ones(l,1),-ones(l,1),ones(l,1)]];
    else

    end
end



%% 


function trans_matrix=exp2sinu(l)
trans_matrix=zeros(2*l+1,2*l+1);
site1=[2*l+1:-1:l+1;1:l+1];
site1=[site1,[l+2:2*l+1;l+2:2*l+1]];
if l>0
sitei=[l:-1:1;l+2:2*l+1];
sitemi=[1:l;1:l];
trans_matrix(site1(1,:)+(site1(2,:)-1)*(2*l+1))=1/sqrt(2);
trans_matrix(sitei(1,:)+(sitei(2,:)-1)*(2*l+1))=1j/sqrt(2);  
trans_matrix(sitemi(1,:)+(sitemi(2,:)-1)*(2*l+1))=-1j/sqrt(2);
trans_matrix(l+1,l+1)=1;
else
trans_matrix=1;
end





