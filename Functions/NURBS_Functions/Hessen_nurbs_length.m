function hessen=Hessen_nurbs_length(xi_e,lambda,xi_s,Knot,coeff,val0)
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

xi_e=xi_e(:);

dim=size(coeff,2)-1;
Npara=length(xi_e);
Nknot=length(Knot);
p=2;
Ncp=Nknot-p-1; %
s1=Func_ParaSpan(xi_e,Knot); % 

[bsp_bas,dbsp_bas,d2bsp_bas]= Func_Basisfun(s1, xi_e, p, Knot); 
xi2cpn=s1-flip(0:p)+1;  
coeff(:,end)=exp(coeff(:,end)); 



l_xi=size(xi_e,1);
ncp=size(coeff,1);
N=sparse(xi2cpn(:),repmat(1:l_xi,[1,p+1]),bsp_bas{1,p+1}(:),ncp,l_xi);
dN=sparse(xi2cpn(:),repmat(1:l_xi,[1,p+1]),dbsp_bas{1,p+1}(:),ncp,l_xi);
d2N=sparse(xi2cpn(:),repmat(1:l_xi,[1,p+1]),d2bsp_bas{1,p+1}(:),ncp,l_xi);
pu=coeff(:,1);
pv=coeff(:,2);
w=coeff(:,3);

m_i=N.*w;
nu_i=N.*(pu.*w);
nv_i=N.*(pv.*w);
nu=sum(nu_i,1);
nv=sum(nv_i,1);
m=sum(m_i,1);
Cu=nu./m; % 曲线u分量
Cv=nv./m; % 曲线v分量


    dm_i=dN.*w;
    dnu_i=dN.*(pu.*w);
    dnv_i=dN.*(pv.*w);

    dnu=sum(dnu_i,1);
    dnv=sum(dnv_i,1);
    dm=sum(dm_i,1);
    dCu=(dnu-dm.*Cu)./m; %曲线u分量相对参数的导数 C'_u
    dCv=(dnv-dm.*Cv)./m; %曲线v分量相对参数的导数 C'_v


    d2m_i=d2N.*w;
    d2nu_i=d2N.*(pu.*w);
    d2nv_i=d2N.*(pv.*w);
    d2nu=sum(d2nu_i,1);
    d2nv=sum(d2nv_i,1);
    d2m=sum(d2m_i,1);

    d2Cu=(d2nu-d2m.*Cu-2*dm.*dCu)./m;
    d2Cv=(d2nv-d2m.*Cv-2*dm.*dCv)./m;
    dCurve=sqrt(dCu.^2+dCv.^2);
hessen_tmp=(dCu.*d2Cu+dCv.*d2Cv)./dCurve;

parfor c1=1:length(xi_e)
    dx(c1,1)=integral(@(xi)Func_dcurve(xi,Knot,coeff),xi_s(c1),xi_e(c1),ArrayValued=true)-val0(c1);
end

    hessen_tmp=hessen_tmp.*dx+dCurve.^2;


hessen=spdiags(hessen_tmp(:),0,length(Cu),length(Cu));
    


end