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

Order_z=20;
Order_phi=20;
deltah=0.0001;   %
r_n=0.450+(0.5:60)*deltah; 
l_z=1.4;


%% 
rang.phi=(1:40)/40*2*pi;
rang.z=(0:50)/50;

%
rang.m=1:Order_phi;
rang.n=1:Order_z;

[MM,PPhi]=ndgrid(rang.m,rang.phi);
mphi=MM.*PPhi;
mphi(MM<=0)=mphi(MM<=0)-pi/2;

[NN,ZZ]=ndgrid(rang.n,rang.z);

nzpi=NN.*ZZ*pi;



for c1=1:size(Cell_wire_eddy,1)
    Temp_wireposition=Cell_wire_eddy{c1,2};

    [Temp_wireposition(:,1),Temp_wireposition(:,2)]=pol2cart(Temp_wireposition(:,2),Temp_wireposition(:,1));
    if c1==1
        segment_position=[];
        wire_segment=[];
    end
    segment_position=[segment_position;(Temp_wireposition(2:end,:)+Temp_wireposition(1:end-1,:))/2];
    wire_segment=[wire_segment;Cell_wire_eddy{c1,1}*(Temp_wireposition(2:end,:)-Temp_wireposition(1:end-1,:))];
end

[segment_position(:,2),segment_position(:,1)]=cart2pol(segment_position(:,1),segment_position(:,2));
Num_seg=[0,5000:5000:size(segment_position,1),size(segment_position,1)*ones(mod(size(segment_position,1),5000)~=0)];

r_n=0.450+(0.5:60)*0; %

[Z_fp,phi_fp]=ndgrid(l_z(1)*(rang.z-0.5),rang.phi);

rang.m=1:Order_phi;
rang.n=1:Order_z;

[MM,PPhi]=ndgrid(rang.m,rang.phi);
mphi=MM.*PPhi;
mphi(MM<=0)=mphi(MM<=0)-pi/2;

[NN,ZZ]=ndgrid(rang.n,rang.z);
nzpi=NN.*ZZ*pi;
Driving_c=cell(length(r_n),1);

Cosmphi=cos(mphi.');
SinmphiM=sin(mphi.').*abs(MM.');

Cosnz=cos(nzpi);
Sinnz=sin(nzpi);
Sinnzn=sin(nzpi)./NN;
r_n=0.450+(0.5:60)*deltah; %

parfor c2=1:length(r_n)
    RRho=r_n(c2);
    Driving0=zeros(length(rang.n),length(rang.m));
    for c1=1:length(Num_seg)-1
        sub1=Num_seg(c1)+1:Num_seg(c1+1);

        r_w=segment_position(sub1,1);
        phi_w=segment_position(sub1,2);
        z_w=segment_position(sub1,3);
        dl=wire_segment(sub1,:);


        Dis_si=((z_w-Z_fp(:).').^2+(r_w.^2+RRho^2-2*r_w*RRho.*cos(phi_w-phi_fp(:).'))).^(-1/2);
        A_phi=(-dl(:,1).*sin(phi_fp(:)).'+dl(:,2).*cos(phi_fp(:)).').*Dis_si;    %
        A_phi=reshape(RRho*l_z(1)*sum(A_phi,1),size(Z_fp));

        A_z=dl(:,3).*Dis_si;    %
        A_z=reshape(l_z(1)*sum(A_z,1),size(Z_fp));


       Drive_tmp1=Cosnz*(A_phi*Cosmphi);
        if any(rang.m==0)
       Drive_tmp1(:,rang.m==0)=sum(Sinnz*A_phi,2);
        end

       Drive_tmp2=l_z(1)/pi*Sinnzn*(A_z*SinmphiM);

       Drive_tmp=Drive_tmp1+Drive_tmp2;
       Driving0=Driving0+Drive_tmp;

    end
    Driving_c{c2,1}=1e-7*(rang.phi(2)-rang.phi(1))*(rang.z(2)-rang.z(1))*Driving0;
end


%% 

    J_phi=cell(length(r_n),1);
    J_z=cell(length(r_n),1);
    Stream=cell(length(r_n),1);
Bz=zeros(26,200);
r_t=0.25;
% r_screen=0.453;
z_t=[0.25,0.125,0,-0.125,-0.25];
phi_t=(0:8)/16*pi;%[0,pi/4,pi/2];
Bz_time=zeros(length(z_t)*length(phi_t),200);
[harmonicorder,sinuorder]=Func_Wire2HarmonicOrder(Cell_wire_eddy,0.2);
    [~,ind]=max(abs(harmonicorder(1:4,3)));
m_range={-5:5,-11:-1,0,1:11};

for m=m_range{ind}
    Driving=cell2mat(Driving_c);
    Driving_m=Driving(:,rang.m==m);

    
    [Muatul_Inductance, Resistance] = Mat_coeff(Order_z,r_n,l_z,m);
    
    Muatul_Inductance=real(Muatul_Inductance);
    Muatul_Inductance=(Muatul_Inductance+Muatul_Inductance.')/2;
    Mat_Mutual_inv=Muatul_Inductance^(-1);
    U0=Mat_Mutual_inv.*cell2mat(Resistance(:,1)).';   
      [E_vector,E_value]=eig(U0);
    e_d=real(diag(E_value));



    Driving_e=(E_vector.')*Mat_Mutual_inv*Driving_m;

    t=(0:0.005:0.25);
    rise_time=1e-3;
    ind=find(e_d>2/rise_time);
    J0=1/rise_time*Driving_e./e_d.*(1-exp(-rise_time*e_d));
%    J0(ind)=0;
    e_d0=e_d;
 %   e_d0(ind)=0;
    J_t=real(J0).*exp(kron(e_d0,-t));


    % 
    J_real_t=E_vector*J_t;


k=-100:0.2:100;
[Basis_sin,Basis_cos]=Func_BasisSinu(k,l_z,Order_z);


K_fft=exp(1j*k.*(z_t(:)+l_z/2));

Kern=zeros(length(r_n),length(k));
ind_tmp=find(k==0);
parfor cl1=1:length(r_n)
    r_screen=0.45+0.0001*(cl1-0.5);
    Kern_tmp=r_screen*abs(k).*besseli(m,abs(k*r_t),1).*(besselk(m-1,abs(r_screen*k),1)+besselk(m+1,abs(r_screen*k),1)).*exp((r_t-r_screen)*abs(k));
    Kern_tmp(ind_tmp)=(r_t/r_screen)^abs(m);
    Kern(cl1,:)=Kern_tmp;
end

    F_sin=abs(k(2)-k(1))*reshape(matkron(Basis_sin,K_fft)*Kern.',[Order_z,length(z_t),length(r_n)]);   
    F_cos=abs(k(2)-k(1))*reshape(matkron(Basis_cos,K_fft)*Kern.',[Order_z,length(z_t),length(r_n)]);
    F_sin=reshape(permute(F_sin,[1,3,2]),[Order_z*length(r_n),length(z_t)]);
    F_cos=reshape(permute(F_cos,[1,3,2]),[Order_z*length(r_n),length(z_t)]);
for n_time=1:length(t)
    J_n=J_real_t(:,n_time); %

if m~=0
    Bz_t=sqrt(2*pi)*1e-7*J_n.'*F_cos;
else
    Bz_t=sqrt(2*pi)*1e-7*J_n.'*F_sin;
end
    Bz_w=kron(Bz_t(:).',cos(m*phi_t));

    Bz_time(:,n_time) =Bz_time(:,n_time)+Bz_w(:);
end
end
%% 



   
%% 绘图



G=48.0e-6;
display(num2str(100*max(abs(Bz_time(:,2)/G/r_t)),2))

% figure
% plot(0.005*(0:125),100*Bz_time(1,1:126)/(G*r_t))
% hold on
% plot(0.005*(0:125),100*Bz_time(9,1:126)/(G*r_t))
% plot(0.005*(0:125),100*Bz_time(17,1:126)/(G*r_t))
% plot(0.005*(0:125),100*Bz_time(25,1:126)/(G*r_t))
% plot(0.005*(0:125),100*Bz_time(33,1:126)/(G*r_t))
% 
% xlabel('t/s')
% ylabel('x relative field/%')
% 
% figure
% G=48.8e-6;
% 
% plot(0.005*(0:125),100*Bz_time(2,1:126)/(G*r_t))
% hold on
% plot(0.005*(0:125),100*Bz_time(10,1:126)/(G*r_t))
% plot(0.005*(0:125),100*Bz_time(18,1:126)/(G*r_t))
% 
% xlabel('t/s')
% ylabel('x relative field/%')

