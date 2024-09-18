function [Muatul_Inductance, Resistance] = Mat_coeff(Order_z, r_n, l_z, m)
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


Muatul_Ind=cell(length(r_n),length(r_n));
rang_nz=1:Order_z;
M_tmp=1j*l_z^2/sqrt(2*pi^5)*(1+(-1).^(rang_nz+1))./rang_nz.^2;
M_tmp=kron(M_tmp.',M_tmp);
k=-100:0.2:100;
[Basis_sin,Basis_cos]=Func_BasisSinu(k,l_z,Order_z);



for c1=1:length(r_n)
    for c2=c1:length(r_n)
        r1=r_n(c1);
        r2=r_n(c2);
        m=abs(m);

        dBessel=(besseli(m+1,abs(r1*k),1)+besseli(m-1,abs(r1*k),1)).*(besselk(m-1,abs(r2*k),1)+besselk(m+1,abs(r2*k),1)).*exp(-(r2-r1)*abs(k))/4; %  -I'(r1 k)K'(r2 k) for k>0
        
        if m==0
            dBessel(k==0)=r1/r2/2;
            Muatul_Inductance=8*pi^2*r1*r2*(Basis_sin.*repmat(dBessel,[Order_z,1]))*(Basis_sin');

        else
            dBessel(k==0)=0;
            Muatul_Inductance=4*pi^2*r1*r2*(Basis_cos.*repmat(dBessel,[Order_z,1]))*(Basis_cos');
            Muatul_Inductance=Muatul_Inductance+4*pi^2*r1*r2*m/2*r1^(m-1)/r2^(m+1)*M_tmp;  
        end

        Muatul_Ind{c1,c2}=1e-7*Muatul_Inductance*abs(k(2)-k(1));
    end
end
[NN,MM]=ndgrid(1:Order_z,m);
for c1=1:length(r_n)

    Resistance{c1,1}=1/(r_n(2)-r_n(1))/37.7e6*pi*l_z*r_n(c1)/2.*((MM(:).*l_z(1)./NN(:)/pi/r_n(c1)).^2+(1+(MM(:)==0))); %电阻

end

for c1=1:length(r_n)
    for c2=c1+1:length(r_n)
        Muatul_Ind{c2,c1}=Muatul_Ind{c1,c2};
    end
end
Muatul_Inductance=cell2mat(Muatul_Ind);
end