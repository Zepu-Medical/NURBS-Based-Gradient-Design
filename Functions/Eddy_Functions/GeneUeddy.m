function U_eddy = GeneUeddy(rang, phi_t, z_t,r_t, m_range, r_n,l_z)
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

[NN,ZZ]=ndgrid(rang.n,rang.z);
nzpi=NN.*ZZ*pi;

Cosnz=cos(nzpi);
Cosnzn=cos(nzpi)./NN;
Sinnz=sin(nzpi);
Sinnzn=sin(nzpi)./NN;

U_eddy=zeros(length(phi_t)*length(z_t),2*length(rang.z)*length(rang.phi));

for m=m_range
    
    mphi=m*rang.phi;
    if m<0
        mphi(m<=0)=mphi(m<=0)-pi/2;
    end
    RRho=r_n(fix(end/2));
    Cosmphi=cos(mphi);
        SinmphiM=sin(mphi)*abs(m);
        if m==0
            Proj_phi=kron(Sinnz,RRho*l_z(1)*(rang.phi(2)-rang.phi(1))*(rang.z(2)-rang.z(1))*Cosmphi);
        else
            Proj_phi=kron(Cosnz,RRho*l_z(1)*(rang.phi(2)-rang.phi(1))*(rang.z(2)-rang.z(1))*Cosmphi);
        end
        Proj_z=kron(Sinnzn,(rang.phi(2)-rang.phi(1))*(rang.z(2)-rang.z(1))*l_z(1)^2/pi*SinmphiM);
        Proj_mat=[Proj_phi,Proj_z];



    [Muatul_Inductance, Resistance] = Mat_coeff(rang.n(end),r_n,l_z,m);
    Muatul_Inductance=real(Muatul_Inductance);
    Muatul_Inductance=(Muatul_Inductance+Muatul_Inductance.')/2;
    Mat_Mutual_inv=Muatul_Inductance^(-1);
    U0=Mat_Mutual_inv.*cell2mat(Resistance(:,1)).';   
    
    [E_vector,E_value]=eig(U0);
    e_d=real(diag(E_value));
    
    rise_time=1e-3;
    e_d0=e_d;
    
    t_measure=0.005;
    U_t=1/rise_time./e_d.*(1-exp(-rise_time*e_d)).*exp(kron(e_d0,-t_measure));
    
    k=-100:0.2:100;
    [Basis_sin,Basis_cos]=Func_BasisSinu(k,l_z,rang.n(end));
    
    
    K_fft=exp(1j*k.*(z_t(:)+l_z/2));
    mphi_t=m*phi_t;
    if m<0
        mphi_t=mphi_t-pi/2;
    end
    Kern=zeros(length(r_n),length(k));
    ind_tmp=k==0;
    for cl1=1:length(r_n)
        r_screen=0.45+0.0001*(cl1-0.5);
        Kern_tmp=r_screen*abs(k).*besseli(m,abs(k*r_t),1).*(besselk(m-1,abs(r_screen*k),1)+besselk(m+1,abs(r_screen*k),1)).*exp((r_t-r_screen)*abs(k));
        Kern_tmp(ind_tmp)=(r_t/r_screen)^abs(m);
        Kern(cl1,:)=Kern_tmp;
    end
    if m~=0
        F_cos=abs(k(2)-k(1))*reshape(matkron(Basis_cos,K_fft)*Kern.',[rang.n(end),length(z_t),length(r_n)]);
        F_cos=reshape(permute(F_cos,[1,3,2]),[rang.n(end)*length(r_n),length(z_t)]).';
        F_mat=sqrt(2*pi)*1e-7*F_cos;
    else
        F_sin=abs(k(2)-k(1))*reshape(matkron(Basis_sin,K_fft)*Kern.',[rang.n(end),length(z_t),length(r_n)]);   
        F_sin=reshape(permute(F_sin,[1,3,2]),[rang.n(end)*length(r_n),length(z_t)]).';
        F_mat=sqrt(2*pi)*1e-7*F_sin;
    end
    
    
    U_eddy_tmp=F_mat*E_vector*(U_t(:).*E_vector.')*Mat_Mutual_inv;
    
    U_eddy_tmp=sum(reshape(U_eddy_tmp,[size(U_eddy_tmp,1),rang.n(end),length(r_n)]),3);
    
    U_eddy_tmp=kron(U_eddy_tmp,cos(m*phi_t.'));
    
    U_eddy_tmp=U_eddy_tmp*Proj_mat;
    
    U_eddy=U_eddy+U_eddy_tmp;
end
end