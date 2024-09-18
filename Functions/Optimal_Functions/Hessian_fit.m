function [Hessen_obj]=Hessian_fit(coeff,lambda,FitPara,ConsPara)
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

N_dC=ConsPara.N_dC;
dN_dC=ConsPara.dN_dC;
d2N_dC=ConsPara.d2N_dC;

    l_xi=numel(coeff)-3*size(N_dC,1);
    ncp=size(N_dC,1);

xi=coeff(1:l_xi,1);
p=2;
    pu=coeff(l_xi+(1:ncp),1);
    pv=coeff(l_xi+ncp+(1:ncp),1);
    w=exp(coeff(l_xi+2*ncp+(1:ncp),1));
    coeffs=[pu,pv,w];


            C2_dC=NURBSgenerator({pu,pv,w},{N_dC,dN_dC,d2N_dC},[1,1],2,2);
            sz0=size(pu,1);

            dx=C2_dC.dCu*1;
            dy=C2_dC.dCv*1;
            ndC=sqrt(dx.^2+dy.^2);
            d2x=C2_dC.d2Cu*1;
            d2y=C2_dC.d2Cv*1;
            dphi=(dx.*d2y-dy.*d2x);
            % d2nurbs=C2.

            % Curvature=(dphi./ndC.^3).';
                gdx=1*C2_dC.gdCu;
                gdy=1*C2_dC.gdCv;
                gd2x=1*C2_dC.gd2Cu;
                gd2y=1*C2_dC.gd2Cv;

                gdphi=(gdx.*d2y+dx.*gd2y-gdy.*d2x-dy.*gd2x);

                gdn2=dx.*gdx+dy.*gdy;

                hdx=1*C2_dC.hdCu;
                hdy=1*C2_dC.hdCv;
                hd2x=1*C2_dC.hd2Cu;
                hd2y=1*C2_dC.hd2Cv;

                hdphi=2*matkron(gdx,gd2y)-2*matkron(gd2x,gdy)+hdx.*d2y+dx.*hd2y-hd2x.*dy-d2x.*hdy;
                
                hCurva_tmp=hdphi./ndC.^3-6*matkron(gdn2./ndC.^5,gdphi)-3*dphi./ndC.^5.*(matkron(gdx,gdx)+matkron(gdy,gdy)+dx.*hdx+dy.*hdy)+15.*dphi./ndC.^7.*matkron(gdn2,gdn2);
if ~isempty(ConsPara.curvature)
h_curva=[hCurva_tmp,-hCurva_tmp]*lambda.ineqnonlin(:);

[id1,~,val]=find(h_curva);

[idx,idy]=ind2sub([ncp*3,ncp*3],id1);

Hessen_cons1=sparse(idx+l_xi,idy+l_xi,val,numel(coeff),numel(coeff));
else
    Hessen_cons1=sparse(numel(coeff),numel(coeff));
end
%% objective 
s1=Func_ParaSpan(xi,Knots); %
[bsp_bas,dbsp_bas,d2bsp_bas]= Func_Basisfun(s1, xi, p, Knots); %B-spline的基函数及其一阶导
xi2cpn=s1-flip(0:p)+1;  %

N=sparse(xi2cpn(:),repmat(1:l_xi,[1,2+1]),bsp_bas{(p+1)}(:),ncp,l_xi);
dN=sparse(xi2cpn(:),repmat(1:l_xi,[1,2+1]),dbsp_bas{(p+1)}(:),ncp,l_xi);
d2N=sparse(xi2cpn(:),repmat(1:l_xi,[1,2+1]),d2bsp_bas{(p+1)}(:),ncp,l_xi);

C2=NURBSgenerator({pu,pv,w},{N,dN,d2N},[1,1],2,2);

DeltaCu=C2.Cu-wire_point(:,1).';
DeltaCv=C2.Cv-wire_point(:,2).';

alpha=0.01;
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
    
    gx=[spdiags(C2.dCu(:),0,l_xi,l_xi);C2.gCu];
    gy=[spdiags(C2.dCv(:),0,l_xi,l_xi);C2.gCv];

Hessen_obj1=gx*(gx).'+gy*(gy).';

Hessen_obj2_tmp4=reshape(C2.hCu*DeltaCu.'+C2.hCv*DeltaCv.',[3*ncp,3*ncp]);
if condition==2
    mat_tmp=[ones(ncp-1,1);0;ones(ncp-1,1);0];
    Hessen_obj2_tmp4=Hessen_obj2_tmp4+spdiags(alpha*[mat_tmp,ones(2*ncp,1),mat_tmp;zeros(ncp,3)],0,3*ncp,3*ncp);
end
Hessen_obj2_tmp1=spdiags((C2.d2Cu.*DeltaCu+C2.d2Cv.*DeltaCv).',0,l_xi,l_xi);
Hessen_obj2_tmp3=C2.gdCu.*DeltaCu+C2.gdCv.*DeltaCv;
Hessen_obj2=[Hessen_obj2_tmp1,Hessen_obj2_tmp3.';Hessen_obj2_tmp3,Hessen_obj2_tmp4];




Hessen_obj=Hessen_obj1+Hessen_obj2+Hessen_cons1;




