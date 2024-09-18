function [nurbs, dnurbs,d2nurbs] = Func_NURBSCurve(xi, Knot,coeff, p,WeightType)
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

xi=xi(:);

dim=size(coeff,2)-1;
Npara=length(xi);
Nknot=length(Knot);
Ncp=Nknot-p-1; %
s1=Func_ParaSpan(xi,Knot); % 
% ctrl_p_xi=s1-flip(0:p)+1; %
[bsp_bas,dbsp_bas,d2bsp_bas]= Func_Basisfun(s1, xi, p, Knot); %
xi2cpn=s1-flip(0:p)+1;  % 
if exist('WeightType','var')
    if strcmpi(WeightType,'Exponential')
    coeff(:,dim+1)=exp(coeff(:,dim+1)); %
    end
end


l_xi=size(xi,1);
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
nurbs=full([Cu;Cv].');




if nargout>1


    dm_i=dN.*w;
    dnu_i=dN.*(pu.*w);
    dnv_i=dN.*(pv.*w);

    dnu=sum(dnu_i,1);
    dnv=sum(dnv_i,1);
    dm=sum(dm_i,1);
    dCu=(dnu-dm.*Cu)./m; 
    dCv=(dnv-dm.*Cv)./m; 
    dnurbs=full([dCu;dCv].');

end

if nargout>2 && p>=2

    d2m_i=d2N.*w;
    d2nu_i=d2N.*(pu.*w);
    d2nv_i=d2N.*(pv.*w);
    d2nu=sum(d2nu_i,1);
    d2nv=sum(d2nv_i,1);
    d2m=sum(d2m_i,1);

    d2Cu=(d2nu-d2m.*Cu-2*dm.*dCu)./m;
    d2Cv=(d2nv-d2m.*Cv-2*dm.*dCv)./m;
    d2nurbs=full([d2Cu;d2Cv].');
    
end
end