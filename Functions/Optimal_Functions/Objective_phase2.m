function [f,gf]=Objective_phase2(coeff_wire,OptPara,WirePara,n)
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

global glo_coeff count_iter glo_f_reg

count_iter=count_iter+1;
if mod(count_iter,n)==0
    glo_coeff{end+1,1}=coeff_wire;
end

if strcmpi(OptPara.Method,'Sigmoid of Distance')
    sz=OptPara.sz;
    d0=OptPara.d0;
    ratio=OptPara.ratio;
    coeff_wire_cell=mat2cell(coeff_wire,OptPara.Cp_Num,[1,1,1]);
    N0=OptPara.N0;
    N1=OptPara.N1;
    dN0=OptPara.dN0;
    dN1=OptPara.dN1;
    dis_para=OptPara.dis_para;
    f_i=zeros(size(dis_para,1),1);
    gf0=cell(size(dis_para,1),1);
    gf1=cell(size(dis_para,1),1);

    parfor c1=1:size(dis_para,1)    %parfor

        p0u=coeff_wire_cell{dis_para{c1,1},1};
        p0v=coeff_wire_cell{dis_para{c1,1},2};
        w0=exp(coeff_wire_cell{dis_para{c1,1},3});

        p1u=coeff_wire_cell{dis_para{c1,3},1};
        p1v=coeff_wire_cell{dis_para{c1,3},2};
        w1=exp(coeff_wire_cell{dis_para{c1,3},3});


        m0_i=N0{c1,1}.*w0;
        dm0_i=dN0{c1,1}.*w0;

        n0u_i=N0{c1,1}.*(p0u.*w0);
        dn0u_i=dN0{c1,1}.*(p0u.*w0);
        n0v_i=N0{c1,1}.*(p0v.*w0);
        dn0v_i=dN0{c1,1}.*(p0v.*w0);

        n0u=sum(n0u_i,1);
        dn0u=sum(dn0u_i,1);
        n0v=sum(n0v_i,1);
        dn0v=sum(dn0v_i,1);

        m0=sum(m0_i,1);
        dm0=sum(dm0_i,1);

        C0u=n0u./m0; 
        C0v=n0v./m0; 
        dC0u=(dn0u-dm0.*C0u)./m0; 
        dC0v=(dn0v-dm0.*C0v)./m0; 
       

        m1_i=N1{c1,1}.*w1;
        dm1_i=dN1{c1,1}.*w1;

        n1u_i=N1{c1,1}.*(p1u.*w1);
        dn1u_i=dN1{c1,1}.*(p1u.*w1);
        n1v_i=N1{c1,1}.*(p1v.*w1);
        dn1v_i=dN1{c1,1}.*(p1v.*w1);

        n1u=sum(n1u_i,1);
        dn1u=sum(dn1u_i,1);
        n1v=sum(n1v_i,1);
        dn1v=sum(dn1v_i,1);

        m1=sum(m1_i,1);
        dm1=sum(dm1_i,1);

        C1u=n1u./m1; 
        C1v=n1v./m1; 
        dC1u=(dn1u-dm1.*C1u)./m1; 
        dC1v=(dn1v-dm1.*C1v)./m1; 



        N0w_m=N0{c1,1}./m0.*w0;

        C0upu=N0w_m;  
        C0vpv=N0w_m;
        C0upv=sparse(size(N0w_m,1),size(N0w_m,2));
        C0vpu=sparse(size(N0w_m,1),size(N0w_m,2));
        C0uw=(p0u-C0u).*N0w_m;
        C0vw=(p0v-C0v).*N0w_m;
        





        N1w_m=N1{c1,1}./m1.*w1;

        C1upu=N1w_m;  
        C1vpv=N1w_m;
        C1upv=sparse(size(N1w_m,1),size(N1w_m,2));
        C1vpu=sparse(size(N1w_m,1),size(N1w_m,2));
        C1uw=(p1u-C1u).*N1w_m;
        C1vw=(p1v-C1v).*N1w_m;


        dN1w_m=dm1_i./m1-dm1./m1.*N1w_m;
        dC1upu=dN1w_m; 
        dC1upv=sparse(size(dN1w_m,1),size(dN1w_m,2));
        dC1uw=-dC1u.*N1w_m+(p1u-C1u).*dN1w_m;

        dC1vpu=sparse(size(dN1w_m,1),size(dN1w_m,2));
        dC1vpv=dN1w_m;  
        dC1vw=-dC1v.*N1w_m+(p1v-C1v).*dN1w_m;

% 二阶导
        sz0=[size(N0{c1,1},1),size(N0{c1,1},2)];

        

        tmp1_2d=kronicdelta(N0w_m);   

        tmp2=matkron(N0w_m,N0w_m);    
        
        
        C0upuw=tmp1_2d-tmp2;
        C0vpvw=C0upuw;


        tmp1_2du=kronicdelta((p0u(:)-C0u).*N0w_m);  
        tmp1_2dv=kronicdelta((p0v(:)-C0v).*N0w_m);
        [pup,puq]=ndgrid(p0u,p0u);
        [pvp,pvq]=ndgrid(p0v,p0v);

        C0uww=tmp1_2du-(pup(:)+puq(:)-2*C0u).*tmp2;
        C0vww=tmp1_2dv-(pvp(:)+pvq(:)-2*C0v).*tmp2;


        sz1=[size(N1{c1,1},1),size(N1{c1,1},2)];


        tmp1_2d=kronicdelta(N1w_m);
       
        tmp2=matkron(N1w_m,N1w_m);   
        
        
        C1upuw=tmp1_2d-tmp2;
        C1vpvw=C1upuw;


        tmp1_2du=kronicdelta((p1u(:)-C1u).*N1w_m);  

        tmp1_2dv=kronicdelta((p1v(:)-C1v).*N1w_m);
        [pup,puq]=ndgrid(p1u,p1u);
        [pvp,pvq]=ndgrid(p1v,p1v);

        C1uww=tmp1_2du-(pup(:)+puq(:)-2*C1u).*tmp2;
        C1vww=tmp1_2dv-(pvp(:)+pvq(:)-2*C1v).*tmp2;


        dtmp1_2d=kronicdelta(dN1w_m);
        dtmp2=matkron(dN1w_m,N1w_m)+matkron(N1w_m,dN1w_m);
        dtmp1_2du=kronicdelta((p1u-C1u).*dN1w_m);
        dtmp1_2dv=kronicdelta((p1v-C1v).*dN1w_m);
       
        dC1upuw=dtmp1_2d-dtmp2;
        dC1vpvw=dC1upuw;

        dC1uww=dtmp1_2du-dC1u.*(C1upuw-tmp2)-(pup(:)+puq(:)-2*C1u).*dtmp2;
        dC1vww=dtmp1_2dv-dC1v.*(C1vpvw-tmp2)-(pvp(:)+pvq(:)-2*C1v).*dtmp2;


        C0={C0u,C0v};
        C1={C1u,C1v,dC1u,dC1v};

        C0dpdq=cell(2,1);
        C0dpdq{1,1}=hessen_cat([1,3;3,3],{C0upuw;C0uww});
        C0dpdq{2,1}=hessen_cat([2,3;3,3],{C0vpvw;C0vww});
        C1dpdq=cell(4,1);
        C1dpdq{1,1}=hessen_cat([1,3;3,3],{C1upuw;C1uww});
        C1dpdq{2,1}=hessen_cat([2,3;3,3],{C1vpvw;C1vww});
        C1dpdq{3,1}=hessen_cat([1,3;3,3],{dC1upuw;dC1uww});
        C1dpdq{4,1}=hessen_cat([2,3;3,3],{dC1vpvw;dC1vww});
        
        C0dp={[C0upu;C0upv;C0uw];[C0vpu;C0vpv;C0vw]};

        C1dp={[C1upu;C1upv;C1uw];[C1vpu;C1vpv;C1vw];[dC1upu;dC1upv;dC1uw];[dC1vpu;dC1vpv;dC1vw]};
        C0dpC1dq=cell(2,4);
        for c2=1:2
            for c3=1:4
                C0dpC1dq{c2,c3}=matkron(C0dp{c2},C1dp{c3});
            end
        end
        C1dpC1dq=cell(4,4);

        for c2=1:4
            for c3=1:4
                C1dpC1dq{c2,c3}=matkron(C1dp{c2},C1dp{c3});
            end
        end
        %
        V1=[-dC1v*sz(c1,2);dC1u*sz(c1,1)];
        V1n=vecnorm(V1,2,1);
        V1u=V1(1,:)./V1n;
        V1v=V1(2,:)./V1n;
        
        distance=[C1u-C0u;C1v-C0v].*sz(c1,:).';

        dr=(dis_para{c1,5}*(distance(1,:).*V1u+distance(2,:).*V1v)-OptPara.dis_0(c1))./d0(c1);

        alpha_i=exp(dr.');
        alpha_r=1./alpha_i;
        sigmoid_i2=-ratio{c1}.'/d0(c1)./((1+alpha_i).*(1+alpha_r));

        f_i(c1,1)=sum(ratio{c1}.'./(1+alpha_i));
        factor0=dis_para{c1,5}*sz(c1,1)*sz(c1,2);
        d_series=[factor0,1,1,0,0,0,0,0,0,1,0,0,0;-factor0,1,0,0,0,0,1,0,0,1,0,0,0;factor0,1,0,0,0,0,0,1,1,0,0,0,0;-factor0,1,0,1,0,0,0,0,1,0,0,0,0;];
        gd0_s=D_obj(d_series,sz(c1,:),0);
        gd1_s=D_obj(d_series,sz(c1,:),1);

        Vn={1./V1n;1./V1n.^3;1./V1n.^5};

        dis_dp0 = term2val([3*sz0(1),sz0(2)], gd0_s, Vn, C0, C1, C0dpdq, C1dpdq, C1dpC1dq, C0dpC1dq, C0dp, C1dp);

        gf0{c1,1}=dis_dp0*sigmoid_i2;

        dis_dp1 = term2val([3*sz1(1),sz0(2)], gd1_s, Vn, C0, C1, C0dpdq, C1dpdq, C1dpC1dq, C0dpC1dq, C0dp, C1dp);

        gf1{c1,1}=dis_dp1*sigmoid_i2;
        
    end

    f=sum(f_i);
    gf=cell(length(OptPara.Cp_Num),1);
    
    for c1=1:length(OptPara.Cp_Num)
        gf{c1,1}=sparse(OptPara.Cp_Num(c1),3);
    end

    for c1=1:length(gf0)
        if isempty(gf{dis_para{c1,1}})
            gf{dis_para{c1,1}}=reshape(gf0{c1},[length(gf0{c1})/3,3]);
        else
           gf{dis_para{c1,1}}=gf{dis_para{c1,1}}+reshape(gf0{c1},[length(gf0{c1})/3,3]); 
        end
        if isempty(gf{dis_para{c1,3}})
            gf{dis_para{c1,3}}=reshape(gf1{c1},[length(gf1{c1})/3,3]);
        else
            gf{dis_para{c1,3}}=gf{dis_para{c1,3}}+reshape(gf1{c1},[length(gf1{c1})/3,3]);
        end

    end

     gf=reshape(cell2mat(gf),[numel(coeff_wire),1]);


end

end


