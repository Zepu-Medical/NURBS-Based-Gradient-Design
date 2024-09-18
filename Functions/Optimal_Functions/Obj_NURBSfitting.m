function [f,gf]=Obj_NURBSfitting(coeff,FitPara)
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

Knots=FitPara.Knots;
wire_point=FitPara.WirePoint;
condition=FitPara.condition;

l_xi=length(wire_point);
xi=coeff(1:l_xi,1);
p=2;
ncp=length(Knots)-(p+1);
coeffs=reshape(coeff(l_xi+1:end,1),[ncp,3]);


%%
s1=Func_ParaSpan(xi,Knots); %
[bsp_bas,dbsp_bas]= Func_Basisfun(s1, xi, p, Knots); %
xi2cpn=s1-flip(0:p)+1;  %

pu=coeffs(:,1);
pv=coeffs(:,2);
w=exp(coeffs(:,3));
N=sparse(xi2cpn(:),repmat(1:l_xi,[1,2+1]),bsp_bas{(p+1)}(:),ncp,l_xi);
dN=sparse(xi2cpn(:),repmat(1:l_xi,[1,2+1]),dbsp_bas{(p+1)}(:),ncp,l_xi);


C2=NURBSgenerator({pu,pv,w},{N,dN,[]},[1,1],2,1);

DeltaCu=C2.Cu-wire_point(:,1).';
DeltaCv=C2.Cv-wire_point(:,2).';
if isfield(FitPara,'alpha')
    alpha=FitPara.alpha;
else
    alpha=0.01;
end

if condition==1
    f=sum(DeltaCu.^2+DeltaCv.^2,'all')/2; %+3*sum((coeffs(2:end,1:2)-coeffs(1:end-1,1:2)).^2,'all'); %加上控制点坐标权重，防止偏离过大。
elseif condition==2
    f=sum(DeltaCu.^2+DeltaCv.^2,'all')/2+alpha/2*sum((coeffs(2:end,1:2)-coeffs(1:end-1,1:2)).^2,'all'); %加上控制点坐标权重，防止偏离过大。
end



    gf_xi=C2.dCu.*DeltaCu+C2.dCv.*DeltaCv;
    gf_coeff=C2.gCu*DeltaCu.'+C2.gCv*DeltaCv.';

    if condition==2
        dCPu=alpha*(pu(2:end,1)-pu(1:end-1,1));
        dCPv=alpha*(pv(2:end,1)-pv(1:end-1,1));
        gf_u=zeros([ncp,1]);
        gf_v=zeros([ncp,1]);
        gf_u(1:end-1,1)=gf_u(1:end-1,1)+dCPu;
        gf_u(2:end,1)=gf_u(2:end,1)-dCPu;
        gf_v(1:end-1,1)=gf_v(1:end-1,1)+dCPv;
        gf_v(2:end,1)=gf_v(2:end,1)-dCPv;
        gf_coeff=gf_coeff+[gf_u;gf_v;sparse(ncp,1)];
    end

    
    gf=[gf_xi.';gf_coeff];


end