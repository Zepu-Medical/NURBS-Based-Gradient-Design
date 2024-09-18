function [C2,C3]=NURBSgenerator(coeff,bspbas,size_coil,JacobiOrder,DerivationOrder)
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


pu=coeff{1};
pv=coeff{2};
w=coeff{3};
N=bspbas{1};
dN=bspbas{2};
sz1=size(N);
mi=N.*w;
nui=N.*(pu.*w);
nvi=N.*(pv.*w);
nu=sum(nui,1);
nv=sum(nvi,1);
m=sum(mi,1);
Cu=nu./m; % 
Cv=nv./m; % 
C2.Cu=Cu;
C2.Cv=Cv;
if DerivationOrder >=1
    dmi=dN.*w;
    dnui=dN.*(pu.*w);
    dnvi=dN.*(pv.*w);
    dnu=sum(dnui,1);
    dnv=sum(dnvi,1);
    dm=sum(dmi,1);
    dnu_m=dnu./m;
    dnv_m=dnv./m;
    dm_m=dm./m;

    dCu=dnu_m-dm_m.*Cu; 
    dCv=dnv_m-dm_m.*Cv; 

    C2.dCu=dCu;
    C2.dCv=dCv;
       if DerivationOrder >=2
            d2N=bspbas{3};
            d2mi=d2N.*w;
            d2nui=d2N.*(pu.*w);
            d2nvi=d2N.*(pv.*w);
            d2nu=sum(d2nui,1);
            d2nv=sum(d2nvi,1);
            d2m=sum(d2mi,1);
            d2mi_m=d2mi./m;
            d2nu_m=d2nu./m;
            d2nv_m=d2nv./m;
            d2m_m=d2m./m;
            C2.d2Cu=(d2nu_m-d2m_m.*Cu-2*dm_m.*dCu);
            C2.d2Cv=(d2nv_m-d2m_m.*Cv-2*dm_m.*dCv);
        end
end
if JacobiOrder>=1
    % mr1=m.^(-1);
    Nw_m=mi./m;
    dmi_m=dmi./m;
    % Cupu=Nw_m;  
    % Cvpv=Nw_m;
    Cupu=Nw_m;
    Cvpv=Nw_m;
    Cuw=(pu-Cu).*Nw_m;
    Cvw=(pv-Cv).*Nw_m;
    empmat1=sparse(size(N,1),size(N,2));
    C2.gCu=[Cupu;empmat1;Cuw];
    C2.gCv=[empmat1;Cvpv;Cvw];

    if DerivationOrder >=1
        dNw_m=(dmi_m-dm_m.*Nw_m);
        dCupu=dNw_m; %
        dCvpv=dNw_m;  %
        dCuw=-dCu.*Nw_m+(pu-Cu).*dNw_m;
        dCvw=-dCv.*Nw_m+(pv-Cv).*dNw_m;
        C2.gdCu=[dNw_m;empmat1;dCuw];
        C2.gdCv=[empmat1;dNw_m;dCvw];
        if DerivationOrder >=2
            
            d2Nw_m=d2mi_m-dm_m.*dmi_m-d2m_m.*Nw_m+dm_m.^2.*Nw_m-dm_m.*dNw_m;
            
            d2Cuw=-C2.d2Cu.*Nw_m-2*dCu.*dNw_m+(pu-Cu).*d2Nw_m;
            
            d2Cvw=-C2.d2Cv.*Nw_m-2*dCv.*dNw_m+(pv-Cv).*d2Nw_m;
            
            C2.gd2Cu=[d2Nw_m;empmat1;d2Cuw];
            
            C2.gd2Cv=[empmat1;d2Nw_m;d2Cvw];
            
        end

    end

    if JacobiOrder>=2
        sz2=[sz1(1)^2,sz1(2)];

        Nw_ms=matkron(Nw_m,Nw_m);    % 
        Cupuw=kronicdelta(Nw_m)-Nw_ms;
        Cvpvw=Cupuw;
        [pup,puq]=ndgrid(pu,pu);
        [pvp,pvq]=ndgrid(pv,pv);
        Cuww=kronicdelta((pu(:)-Cu).*Nw_m)-(pup(:)+puq(:)-2*Cu).*Nw_ms;
        Cvww=kronicdelta((pv(:)-Cv).*Nw_m)-(pvp(:)+pvq(:)-2*Cv).*Nw_ms;

        C2.hCu=hessen_reorder({[],[],Cupuw;[],[],[];[],[],Cuww});
        C2.hCv=hessen_reorder({[],[],[];[],[],Cvpvw;[],[],Cvww});
        if DerivationOrder >=1
            dNw_ms=matkron(dNw_m,Nw_m)+matkron(Nw_m,dNw_m);
            dCupuw=kronicdelta(dNw_m)-dNw_ms; 
            dCvpvw=dCupuw;
            dCuww=kronicdelta((pu-Cu).*dNw_m)-dCu.*(Cupuw-Nw_ms)-(pup(:)+puq(:)-2*Cu).*dNw_ms;
            dCvww=kronicdelta((pv-Cv).*dNw_m)-dCv.*(Cvpvw-Nw_ms)-(pvp(:)+pvq(:)-2*Cv).*dNw_ms;

            C2.hdCu=hessen_reorder({[],[],dCupuw;[],[],[];[],[],dCuww});
            C2.hdCv=hessen_reorder({[],[],[];[],[],dCvpvw;[],[],dCvww});
            if DerivationOrder >=2
                d2Nw_ms=matkron(d2Nw_m,Nw_m)+2*matkron(dNw_m,dNw_m)+matkron(Nw_m,d2Nw_m);
                d2Cupuw=kronicdelta(d2Nw_m)-d2Nw_ms;
                d2Cvpvw=d2Cupuw;
                d2Cuww=-C2.d2Cu.*(Cupuw-Nw_ms)-2*dCu.*(dCupuw-dNw_ms)+kronicdelta((pu-Cu).*d2Nw_m)-(pup(:)+puq(:)-2*Cu).*d2Nw_ms;
                d2Cvww=-C2.d2Cv.*(Cvpvw-Nw_ms)-2*dCv.*(dCvpvw-dNw_ms)+kronicdelta((pv-Cv).*d2Nw_m)-(pvp(:)+pvq(:)-2*Cv).*d2Nw_ms;

                C2.hd2Cu=hessen_reorder({[],[],d2Cupuw;[],[],[];[],[],d2Cuww});
                C2.hd2Cv=hessen_reorder({[],[],[];[],[],d2Cvpvw;[],[],d2Cvww});
            end
        end
    end
end


if nargout>1
    % 
    rp=(size_coil(1,3)-size_coil(1,1));  %
    zp=(size_coil(1,4)-size_coil(1,2));
    r=size_coil(1,1)+rp*Cv;

    xbar=cos(Cu); %
    ybar=sin(Cu);

    Cx=r.*xbar;
    Cy=r.*ybar;
    Cz=size_coil(1,2)+zp*Cv;

    C3.Cx=Cx;
    C3.Cy=Cy;
    C3.Cz=Cz;
    if DerivationOrder >=1
        dCx=-Cy.*dCu+rp*dCv.*xbar;
        dCy=Cx.*dCu+rp*dCv.*ybar;
        dCz=zp*dCv;

        C3.dCx=dCx;
        C3.dCy=dCy;
        C3.dCz=dCz;
    end

    if JacobiOrder>=1
        Cxpu=-Cy.*Cupu;    %
        Cxpv=rp*Cvpv.*xbar;
        Cxw=-Cy.*Cuw+rp*Cvw.*xbar;

        Cypu=Cx.*Cupu;
        Cypv=rp*Cvpv.*ybar;
        Cyw=Cx.*Cuw+rp*Cvw.*ybar;

        Czpu=sparse(sz1(1),sz1(2));
        Czpv=zp*Cvpv;
        Czw=zp*Cvw;

        C3.gCx=[Cxpu;Cxpv;Cxw];
        C3.gCy=[Cypu;Cypv;Cyw];
        C3.gCz=[Czpu;Czpv;Czw];
        if DerivationOrder >=1
            dCxpu=-Cx.*Cupu.*dCu-Cy.*dCupu-rp*dCv.*ybar.*Cupu;  %
            dCxpv=rp*(-Cvpv.*ybar.*dCu+dCvpv.*xbar);
            dCxw=-dCy.*Cuw-Cy.*dCuw-rp*ybar.*dCu.*Cvw+rp*xbar.*dCvw;

            dCypu=-Cy.*Cupu.*dCu+ Cx.*dCupu+rp.*dCv.*xbar.*Cupu;
            dCypv=rp*(Cvpv.*xbar.*dCu+dCvpv.*ybar);
            dCyw=dCx.*Cuw+Cx.*dCuw+rp*xbar.*dCu.*Cvw+rp*ybar.*dCvw;

            dCzpv=zp*dCvpv;
            dCzpu=sparse(size(Czpv,1),size(Czpv,2));
            dCzw=zp*dCvw;

            C3.gdCx=[dCxpu;dCxpv;dCxw];
            C3.gdCy=[dCypu;dCypv;dCyw];
            C3.gdCz=[dCzpu;dCzpv;dCzw];
        end
        if JacobiOrder>=2
            CvwCuw=matkron(Cvw,Cuw);
            CvpvCuw=matkron(Cvpv,Cuw);
            CupuCyw=matkron(Cupu,Cyw);
            CupuCxw=matkron(Cupu,Cxw);
            CuwCyw=matkron(Cuw,Cyw);
            CuwCxw=matkron(Cuw,Cxw);
            Cxpupu=-matkron(Cupu,Cypu);
            Cxpupv=-matkron(Cupu,Cypv);

            Cxpuw=-CupuCyw-Cy.*Cupuw;
            Cxpvw=-rp*ybar.*CvpvCuw+rp*xbar.*Cvpvw;
            Cxww=-CuwCyw-rp*ybar.*CvwCuw-Cy.*Cuww+rp*xbar.*Cvww;


            Cypupu=matkron(Cupu,Cxpu);
            Cypupv=matkron(Cupu,Cxpv);

            Cypuw=CupuCxw+Cx.*Cupuw;
            Cypvw=rp*xbar.*CvpvCuw+rp*ybar.*Cvpvw;
            Cyww=CuwCxw+rp*xbar.*CvwCuw+Cx.*Cuww+rp*ybar.*Cvww;


            Czpvw=zp*Cvpvw;
            Czww=zp*Cvww;

            C3.hCx=hessen_reorder({Cxpupu,Cxpupv,Cxpuw;[],[],Cxpvw;[],[],Cxww});
            C3.hCy=hessen_reorder({Cypupu,Cypupv,Cypuw;[],[],Cypvw;[],[],Cyww});
            C3.hCz=hessen_reorder({[],[],[];[],[],Czpvw;[],[],Czww});


            if DerivationOrder >=1
                dCuwCyw=matkron(dCuw,Cyw);
                CuwdCyw=matkron(Cuw,dCyw);
                dCvwCuw=matkron(dCvw,Cuw);
                CvwdCuw=matkron(Cvw,dCuw);
                CvpvdCuw=matkron(Cvpv,dCuw);
                dCvpvCuw=matkron(dCvpv,Cuw);

                dCxpupu=-matkron(Cupu,dCypu)-matkron(dCupu,Cypu);
                dCxpupv=-matkron(Cupu,dCypv)-matkron(dCupu,Cypv);
                dCxpvpv=sparse(sz2(1),sz2(2));

                dCxpuw=-matkron(Cupu,dCyw)-matkron(dCupu,Cyw)-dCy.*Cupuw-Cy.*dCupuw;
                dCxpvw=-rp*ybar.*dCvpvCuw-rp*xbar.*dCu.*CvpvCuw-rp*ybar.*CvpvdCuw-rp*ybar.*dCu.*Cvpvw+rp*xbar.*dCvpvw;
                dCxww=-CuwdCyw-dCuwCyw-dCy.*Cuww-Cy.*dCuww-rp*xbar.*dCu.*CvwCuw-rp*ybar.*dCvwCuw-rp*ybar.*CvwdCuw-rp*ybar.*dCu.*Cvww+rp*xbar.*dCvww;
                dCypupu=matkron(Cupu,dCxpu)+matkron(dCupu,Cxpu);
                dCypupv=matkron(Cupu,dCxpv)+matkron(dCupu,Cxpv);

                dCypuw=matkron(Cupu,dCxw)+matkron(dCupu,Cxw)+dCx.*Cupuw+Cx.*dCupuw;
                dCypvw=-rp*ybar.*dCu.*CvpvCuw+rp*xbar.*CvpvdCuw+rp*xbar.*dCvpvCuw+rp*xbar.*dCu.*Cvpvw+rp*ybar.*dCvpvw;
                dCyww=matkron(Cuw,dCxw)+matkron(dCuw,Cxw)+dCx.*Cuww+Cx.*dCuww...
                    -rp*ybar.*dCu.*CvwCuw+rp*xbar.*CvwdCuw+rp*xbar.*dCvwCuw+rp*xbar.*dCu.*Cvww+rp*ybar.*dCvww;

                dCzpvw=zp*dCvpvw;
                dCzww=zp*dCvww;
                C3.hdCx=hessen_reorder({dCxpupu,dCxpupv,dCxpuw;[],dCxpvpv,dCxpvw;[],[],dCxww});
                C3.hdCy=hessen_reorder({dCypupu,dCypupv,dCypuw;[],[],dCypvw;[],[],dCyww});
                C3.hdCz=hessen_reorder({[],[],[];[],[],dCzpvw;[],[],dCzww});

            end
        end

    end
end
end

