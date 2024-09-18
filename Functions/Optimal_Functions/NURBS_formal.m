function NURBS_curve = NURBS_formal(NURBS_curve,FormPara)
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

if isfield(FormPara,'XiForm')

else
    FormPara.XiForm='gaussian';
end

if strcmpi(FormPara.XiForm,'uniform')
    if isfield(FormPara,'xi0')
        xi0=FormPara.xi0(:);
    else
        xi0=linspace(0,1,17).';
        xi0(end)=[];
    end

    Xi0=reshape(xi0+(0:99),[],1);
    gauss_order=length(xi0);
    GaussWeight0=repmat(ones(gauss_order,1),[100,1]);
elseif strcmpi(FormPara.XiForm,'gaussian')
    load('Gauss_coeff_1d.mat','gauss_1d_order16');
    gauss_order=16;
    Xi0=reshape((gauss_1d_order16(:,1)+1)/2+(0:99),[16*100,1]);
    GaussWeight0=repmat(gauss_1d_order16(:,2),[100,1]);
end

% 
for c1=1:length(NURBS_curve)
    for c2=1:length(NURBS_curve(c1).Knot)
        n=NURBS_curve(c1).Knot{c2,1}(end);
        NURBS_curve(c1).CenterNum(c2,1)=n+1;
        NURBS_curve(c1).p=2;
        if strcmpi(FormPara.XiForm,'uniform')
            NURBS_curve(c1).Xi{c2,1}=Xi0(1:n*gauss_order+1,1);
            elseif strcmpi(FormPara.XiForm,'gaussian')
            NURBS_curve(c1).Xi{c2,1}=Xi0(1:n*gauss_order,1);
        end
        NURBS_curve(c1).GaussWeight{c2,1}=GaussWeight0(1:length(NURBS_curve(c1).Xi{c2,1}),1);
        xi_span=Func_ParaSpan(NURBS_curve(c1).Xi{c2,1},NURBS_curve(c1).Knot{c2,1});
        NURBS_curve(c1).SpanIndex{c2,1}=xi_span; %
        NURBS_curve(c1).Xi2Cpn{c2,1}=xi_span-flip(0:NURBS_curve(c1).p)+1; %
        [bsp_bas,dbsp_bas,d2bsp_bas]= Func_Basisfun(xi_span, NURBS_curve(c1).Xi{c2,1}, NURBS_curve(c1).p, NURBS_curve(c1).Knot{c2,1});
        l_xi=size(NURBS_curve(c1).Xi2Cpn{c2,1},1);
        ncp=size(NURBS_curve(c1).ContrPoint{c2,1},1);
        NURBS_curve(c1).Bsp_Bas{c2,1}=sparse(NURBS_curve(c1).Xi2Cpn{c2,1}(:),repmat(1:l_xi,[1,NURBS_curve(c1).p+1]),bsp_bas{1,NURBS_curve(c1).p+1}(:),ncp,l_xi);
        NURBS_curve(c1).dBsp_Bas{c2,1}=sparse(NURBS_curve(c1).Xi2Cpn{c2,1}(:),repmat(1:l_xi,[1,NURBS_curve(c1).p+1]),dbsp_bas{1,NURBS_curve(c1).p+1}(:),ncp,l_xi);
        NURBS_curve(c1).d2Bsp_Bas{c2,1}=sparse(NURBS_curve(c1).Xi2Cpn{c2,1}(:),repmat(1:l_xi,[1,NURBS_curve(c1).p+1]),d2bsp_bas{1,NURBS_curve(c1).p+1}(:),ncp,l_xi);

    end
end
end