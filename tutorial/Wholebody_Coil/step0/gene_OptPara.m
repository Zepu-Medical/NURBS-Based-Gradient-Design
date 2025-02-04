function OptPara = gene_OptPara(NURBS_curve,WirePara,Method,n_divide)
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


    near_wire(1).Wire_arrange=num2cell([1:18;[2:18,-1];0:17].');
    near_wire(2).Wire_arrange=num2cell([19:28;[20:28,-1];[0,19:27]].');

    NURBS_curve=find_mindistance(NURBS_curve,WirePara,near_wire,n_divide,'on');

p=NURBS_curve(1).p;
if strcmpi(Method,'Sigmoid of Distance')
    OptPara.Method='Sigmoid of Distance';
Cp_Num=zeros(0,1);
id=0;
dis_para=cell(0,5);
d0_rang=[0.0051,0.0095]*0.8; 
% d0_rang=[0.0032,0.0064,0.008,0.004];
ratio_rang=[1,0.4,0.0,0.2]; 
ratio1_rang=[2,1,2,1.5]*50;  
d0=[];
d_reg={};
ratio={};
ratio1=[];
sz=zeros(0,2);
sz1=zeros(0,2);
Knots=cell(0,1);
size_coil=WirePara.size_coil;

for c0=1:length(NURBS_curve)

    tmp_dis=NURBS_curve(c0).dis;
        tmp_sz=repmat([min(size_coil(c0,1),size_coil(c0,3)),abs(size_coil(c0,4)-size_coil(c0,2))],[size(tmp_dis,1),1]);
        tmp_sz1=repmat([min(size_coil(c0,1),size_coil(c0,3)),abs(size_coil(c0,4)-size_coil(c0,2))],[length(NURBS_curve(c0).Xi),1]);
        tmp_d_reg=cell(0,3);
    for c1=1:length(NURBS_curve(c0).Xi)
        Cp_Num(end+1,1)=length(NURBS_curve(c0).Weight{c1});
        Knots{end+1,1}=NURBS_curve(c0).Knot{c1};

        ratio1(end+1,1)=ratio1_rang(c0);

        ind1=[1,3:2:Cp_Num(end)];
        ind2=[2:2:Cp_Num(end)];
        X=tmp_sz1(c1,1)*NURBS_curve(c0).ContrPoint{c1}(:,1);
        Y=tmp_sz1(c1,2)*NURBS_curve(c0).ContrPoint{c1}(:,2);

        tmp_d_reg{c1,1}=sqrt((X(ind2)-X(ind1)).^2+(Y(ind2)-Y(ind1)).^2);
        tmp_d_reg{c1,2}=X(ind2)-X(ind1);
        tmp_d_reg{c1,3}=Y(ind2)-Y(ind1);
    end


    tmp_d0=d0_rang(c0).*ones(size(tmp_dis,1),1);
    tmp_ratio=cell(0,1);
    
    for c1=1:size(tmp_dis,1) 
        tmp_dis{c1,1}=tmp_dis{c1,1}+id;
        tmp_dis{c1,3}=tmp_dis{c1,3}+id;

        tmp_ratio{c1,1}=ratio_rang(c0).*exp(-NURBS_curve(c0).dis{c1,6}/(d0_rang(c0)*2));
    end
    

    dis_para=[dis_para;tmp_dis];
    id=id+length(NURBS_curve(c0).Xi); 
    sz=[sz;tmp_sz];
    sz1=[sz1;tmp_sz1];
    ratio=[ratio;tmp_ratio];
    d_reg=[d_reg;tmp_d_reg];
    d0=[d0;tmp_d0];
    
end

for c1=1:size(dis_para,1)
    s0=Func_ParaSpan(dis_para{c1,2}(:),Knots{dis_para{c1,1}});
    s1=Func_ParaSpan(dis_para{c1,4}(:),Knots{dis_para{c1,3}});
    [bsp_bas0,bsp_dbas0]= Func_Basisfun(s0, dis_para{c1,2}, 2, Knots{dis_para{c1,1}}); %
    [bsp_bas1,bsp_dbas1]= Func_Basisfun(s1, dis_para{c1,4}, 2, Knots{dis_para{c1,3}}); %

    Xi2Cpn0=s0-flip(0:p)+1; %
    Xi2Cpn1=s1-flip(0:p)+1; 

    ncp0=Cp_Num(dis_para{c1,1});   
    ncp1=Cp_Num(dis_para{c1,3});

        
    l0=size(Xi2Cpn0,1);
    l1=size(Xi2Cpn1,1);

    N0{c1,1}=sparse(repmat(1:l0,[1,p+1]),Xi2Cpn0(:),bsp_bas0{1,p+1}(:),l0,ncp0);
    dN0{c1,1}=sparse(repmat(1:l0,[1,p+1]),Xi2Cpn0(:),bsp_dbas0{1,p+1}(:),l0,ncp0);

    N1{c1,1}=sparse(repmat(1:l1,[1,p+1]),Xi2Cpn1(:),bsp_bas1{1,p+1}(:),l1,ncp1);
    dN1{c1,1}=sparse(repmat(1:l1,[1,p+1]),Xi2Cpn1(:),bsp_dbas1{1,p+1}(:),l1,ncp1);

    N0{c1,1}=N0{c1,1}.';
    dN0{c1,1}=dN0{c1,1}.';
    N1{c1,1}=N1{c1,1}.';
    dN1{c1,1}=dN1{c1,1}.';
end
OptPara.ratio1=ratio1;
OptPara.Cp_Num=Cp_Num;
OptPara.dis_para=dis_para;
OptPara.d0=d0/4;
OptPara.dis_0=d0;
OptPara.ratio=ratio;
OptPara.sz=sz;
OptPara.sz1=sz1;
OptPara.d_reg=d_reg;
OptPara.N0=N0;
OptPara.dN0=dN0;
OptPara.N1=N1;
OptPara.dN1=dN1;
end
end