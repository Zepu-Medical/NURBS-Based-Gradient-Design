function FastHenry_import(Cell_wire,filename,WirePara)

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


date_today=datestr(today,'yy_mm_dd');
date_now=datestr(now,'HH_MM_SS');



string1={'X','Y','Z'};
%% 
for c1=1
    if ~exist('Cell_wire')
eval(['Cell_data=[Cell_',string1{c1},'P;Cell_',string1{c1},'S];']);
eval(['n1=length(Cell_',string1{c1},'P);']);
Cell_wire=cell(size(Cell_data,1),3);
for Count1=1:size(Cell_data,1)
    if Count1<n1+1
    Cell_wire{Count1,1}=1;
    else 
    Cell_wire{Count1,1}=-1;
    end
    Cell_wire{Count1,2}=Cell_data{Count1,1};
    Cell_wire{Count1,3}='Cylindrical';
end
    end

fid1 = fopen([filename,'_',string1{c1},'_',date_today,'.inp'],'wt');

fprintf(fid1, '* FastHenry input file\n\n');
fprintf(fid1, '.units mm\n\n');
fprintf(fid1, '.default sigma=58000.0 nhinc=1 nwinc=1 rh=2 rw=2\n\n');
fprintf(fid1, '* Nodes\n');

% 
for Count1=1:size(Cell_wire,1)
    seg_num=size(Cell_wire{Count1,2},1);  %导线点数目
    d=Cell_wire{Count1,2}([1:2,seg_num],:);  % 
    if strcmp(Cell_wire{Count1,3},'Cylindrical')
    d(:,1:2)=[d(:,1).*cos(d(:,2)),d(:,1).*sin(d(:,2))]; 
    end

end

c1=0;

for Count1=1:size(Cell_wire,1)
    
    U=Cell_wire{Count1,2};

    if strcmp(Cell_wire{Count1,3},'Cylindrical')
    U=[U(:,1).*cos(U(:,2)),U(:,1).*sin(U(:,2)),U(:,3)]*1000;
    else
    U=U*1000;
    end

    if Cell_wire{Count1,1}<0
        U=U(end:-1:1,:);
    end
    for c2=1:size(U,1)

    fprintf(fid1, ['NFHNode',num2str(c1),' x=',num2str(U(c2,1),'%.12f'),' y=',num2str(U(c2,2),'%.12f'),' z=',num2str(U(c2,3),'%.12f'),'\n']);
    c1=c1+1;
    end
end

%
fprintf(fid1, '\n* Segments from paths\n');
c1=0;
size_sect=WirePara.size_sect;

for Count1=1:size(Cell_wire,1)
        U=Cell_wire{Count1,2};
        
    for c2=1:size(U,1)
        d0=Cell_wire{Count1,5}(:,3:4)-Cell_wire{Count1,5}(:,1:2);
        clinic=d0(:,1)./d0(:,2);
        surf_id=Cell_wire{Count1,4}(c2,1);
        if c2~=size(U,1)
            u0=(U(c2+1,:)+U(c2,:))/2;
            u0(1,3)=-vecnorm(u0(1,:),2,2)*clinic(c2);
            w_v=cross(U(c2+1,:)-U(c2,:),u0);
            
            
            
            fprintf(fid1, ['EFHPath',num2str(c1),' NFHNode',num2str(c1),' NFHNode',num2str(c1+1),' w=',num2str(1000*size_sect(surf_id,2)),' h=',num2str(1000*size_sect(surf_id,1)), ...
                ' wx=',num2str(w_v(1)),' wy=',num2str(w_v(2)),' wz=',num2str(w_v(3)),'\n']);   
        
        end
    c1=c1+1;
    end
end

fprintf(fid1,'\n* Node shorts\n');
c1=0;
for Count1=1:size(Cell_wire,1)
    U=Cell_wire{Count1,2};
    n=size(U,1);
    if Count1<size(Cell_wire,1)
        fprintf(fid1, ['.equiv NFHNode',num2str(c1+n-1),' NFHNode',num2str(c1+n),'\n']);
    end
    c1=c1+n;
end
%%

fprintf(fid1, '\n* Port\n');
fprintf(fid1, ['.external NFHNode',num2str(0),' NFHNode',num2str(c1-1),'\n\n']);
fprintf(fid1, ['.freq fmin=1.0 fmax=1.0 ndec=1.0\n\n']);
fprintf(fid1, '.end');

fclose(fid1);
end
end