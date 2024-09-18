function [NURBS_curve]=find_mindistance(NURBS_curve,WirePara,near_wire,n_divide,check_sign,divide_type)
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


for c1=1:length(NURBS_curve)
    sc=WirePara.size_coil(c1,:);
    wire_id=NURBS_curve(c1).WireId;
    sz=[min(sc(1),sc(3)),sqrt(abs(sc(2)-sc(4)).^2+abs(sc(1)-sc(3)).^2)];
        tmp_cell=cell(0,5);
    for c2=1:length(near_wire(c1).Wire_arrange(:,1))

        wire_num0=find(NURBS_curve(c1).WireId==near_wire(c1).Wire_arrange{c2,1});
        if ~exist('divide_type','var')
            divide_type='uniform';
        end
        if strcmpi(divide_type,'lorentz')
        xi_single_seg=0.5*log(linspace(0.2,1,round(n_divide/2)+1))/log(0.2);
        
        xi_single_seg=[flip(xi_single_seg),1-xi_single_seg(1,2:end)]';
        xi_single_seg=[0,0.5,1]';
            elseif strcmpi(divide_type,'uniform')
        xi_single_seg=linspace(0,1,n_divide).'; %[0;0.333333;0.666666];
        end
        xi_single_seg(end)=[];
        xi0=xi_single_seg+repmat(0:NURBS_curve(c1).CenterNum(wire_num0)-2,[length(xi_single_seg),1]);
        xi0=[xi0(:).',NURBS_curve(c1).CenterNum(wire_num0)-1];
        coeffs_0=[NURBS_curve(c1).ContrPoint{wire_num0},NURBS_curve(c1).Weight{wire_num0}];
        
        wire_point=Func_NURBSCurve(xi0.', NURBS_curve(c1).Knot{wire_num0},coeffs_0,NURBS_curve(c1).p,'Exponential');
        

        for c3=2:3
                xi_fit=cell(length(near_wire(c1).Wire_arrange{c2,c3}),1);
                sign_dwire=zeros(1,length(near_wire(c1).Wire_arrange{c2,c3}));
                wire_num=cell(length(near_wire(c1).Wire_arrange{c2,c3}),1);

                distance=zeros(length(xi0),length(near_wire(c1).Wire_arrange{c2,c3}));

            for c4=1:length(near_wire(c1).Wire_arrange{c2,c3})
                


                wire_num{c4,1}=find(wire_id==near_wire(c1).Wire_arrange{c2,c3}(c4));

                

                if ~isempty(wire_num{c4,1})
                    xi_0=linspace(0,NURBS_curve(c1).CenterNum(wire_num{c4,1})-1,3*NURBS_curve(c1).CenterNum(wire_num{c4,1})).';
                    [wire,~] = Func_NURBSCurve(xi_0, NURBS_curve(c1).Knot{wire_num{c4,1}},[NURBS_curve(c1).ContrPoint{wire_num{c4,1}},NURBS_curve(c1).Weight{wire_num{c4,1}}], 2,'Exponential');
                    dr=sc(1)^2*(wire(:,1)-wire_point(:,1).').^2+abs(sc(2)-sc(4))^2*(wire(:,2)-wire_point(:,2).').^2;
                    [~,id]=min(dr,[],1);
                    xi_ini=xi_0(id);
                    options=optimoptions('lsqnonlin',FunctionTolerance=1e-9,OptimalityTolerance=1e-9,MaxFunctionEvaluations=4e4,SpecifyObjectiveGradient=true,Display='off');

                    lb=zeros(length(xi_ini),1);

                    rb=(NURBS_curve(c1).CenterNum(wire_num{c4,1})-1)*ones(length(xi_ini),1);

                    coeffs=[NURBS_curve(c1).ContrPoint{wire_num{c4,1}},NURBS_curve(c1).Weight{wire_num{c4,1}}];

                    [xi,~]=lsqnonlin(@(xi)Obj_min_distance(xi,wire_point,NURBS_curve(c1).Knot{wire_num{c4,1}},coeffs,sz),xi_ini,lb,rb,[],[],[],[],[],options);
                
                    [wire_fit, deriv_wire_fit] = Func_NURBSCurve(xi, NURBS_curve(c1).Knot{wire_num{c4,1}},[NURBS_curve(c1).ContrPoint{wire_num{c4,1}},NURBS_curve(c1).Weight{wire_num{c4,1}}], 2,'Exponential');
                    

                    deriv_wire_fit=deriv_wire_fit.*sz;

                    delta_wire=(wire_fit-wire_point).*sz;

                    distance(:,c4)=vecnorm(delta_wire,2,2);

                    temp_sign=sign(-deriv_wire_fit(:,2).*delta_wire(:,1)+deriv_wire_fit(:,1).*delta_wire(:,2));
                    temp_signp=temp_sign/temp_sign(1);

                    if ~all(temp_signp==1)
                        xi_ini(temp_signp~=1)=NURBS_curve(c1).CenterNum(wire_num{c4,1})-1-xi_ini(temp_signp~=1);

                    [xi,~]=lsqnonlin(@(xi)Obj_min_distance(xi,wire_point,NURBS_curve(c1).Knot{wire_num{c4,1}},coeffs,sz),xi_ini,lb,rb,[],[],[],[],[],options);
                
                    [wire_fit, deriv_wire_fit] = Func_NURBSCurve(xi, NURBS_curve(c1).Knot{wire_num{c4,1}},[NURBS_curve(c1).ContrPoint{wire_num{c4,1}},NURBS_curve(c1).Weight{wire_num{c4,1}}], 2,'Exponential');
                    
                    deriv_wire_fit=deriv_wire_fit.*sz;
                    delta_wire=(wire_fit-wire_point).*sz;
                    distance(:,c4)=vecnorm(delta_wire,2,2);
                    temp_sign=sign(-deriv_wire_fit(:,2).*delta_wire(:,1)+deriv_wire_fit(:,1).*delta_wire(:,2));
                    temp_signp(:,c4)=temp_sign/temp_sign(1);
                    if strcmpi(check_sign,'on')
                        if ~all(temp_signp==1)
                            error(['初态不对,电流面',num2str(c1),',导线',num2str(near_wire(c1).Wire_arrange{c2,1}),'和导线',near_wire(c1).Wire_arrange{c2,c3}]);
                        end
                    end
                    end

                    sign_dwire(1,c4)=temp_sign(1);

                    xi_fit{c4,1}=xi;
                end
            end

                [a,b]=sort(distance,2,"ascend");
                for c4=1:length(near_wire(c1).Wire_arrange{c2,c3})
                    if ~isempty(wire_num{c4,1})
                    tmp_id=find(b(:,1)==c4);

                    if ~isempty(tmp_id)

                        wire_tmp={wire_num0,xi0(tmp_id),wire_num{c4,1},xi_fit{c4,1}(tmp_id).',sign_dwire(c4),a(tmp_id,1).'};
                        tmp_cell=[tmp_cell;wire_tmp];

                    end

                end
                end
            end

        
    end
        NURBS_curve(c1).dis=tmp_cell;
        
end

end



function [f,gf]=Obj_min_distance(xi,wire_point,Knots,coeffs,sz)
p=2;
ncp=length(Knots)-(p+1);
l_xi=length(wire_point);
s1=Func_ParaSpan(xi,Knots); %
[bsp_bas,dbsp_bas]= Func_Basisfun(s1, xi, p, Knots); %B-spline的基函数及其一阶导

xi2cpn=s1(:)-flip(0:p)+1;  %

CPu=coeffs(:,1);
CPv=coeffs(:,2);
w=exp(coeffs(:,3));

bsp_bas=sparse(repmat(1:l_xi,[1,2+1]),xi2cpn(:),bsp_bas{(p+1)}(:),l_xi,ncp);
dbsp_bas=sparse(repmat(1:l_xi,[1,2+1]),xi2cpn(:),dbsp_bas{(p+1)}(:),l_xi,ncp);
m_i=bsp_bas.*w.';
nu_i=bsp_bas.*(CPu.*w).';
nv_i=bsp_bas.*(CPv.*w).';

nu=sum(nu_i,2);
nv=sum(nv_i,2);
m=sum(m_i,2);
Cu=nu./m; % 曲线u分量
Cv=nv./m; % 曲线v分量

gf=[];
DeltaCu=Cu-wire_point(:,1);
DeltaCv=Cv-wire_point(:,2);
f=[DeltaCu*sz(1);DeltaCv*sz(2)];
% f=sum(sz(1)^2*DeltaCu.^2+sz(2)^2*DeltaCv.^2,'all')/2; %+3*sum((coeffs(2:end,1:2)-coeffs(1:end-1,1:2)).^2,'all'); %加上控制点坐标权重，防止偏离过大。
if nargout>1
    mr1=1./m;
    mr2=mr1.^2;
    dm_i=dbsp_bas.*w.';
    dnu_i=dbsp_bas.*(CPu.*w).';
    dnv_i=dbsp_bas.*(CPv.*w).';
    dnu=sum(dnu_i,2);
    dnv=sum(dnv_i,2);
    dm=sum(dm_i,2);
    dCu=(dnu-dm.*Cu)./m; %曲线u分量相对参数的导数 C'_u
    dCv=(dnv-dm.*Cv)./m; %曲线v分量相对参数的导数 C'_v
    
%    gf_xi=sz(1)^2*dCu.*DeltaCu+sz(2)^2*dCv.*DeltaCv;

    gf=[diag(sz(1)*dCu);diag(sz(2)*dCv)];
end
end