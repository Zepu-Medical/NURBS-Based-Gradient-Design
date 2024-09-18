function [Cell_magneticfield] = Func_MagneticField(Cell_magneticfield,Cell_wire,fieldcomponent)
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


%
Ctrl_maxfieldpoint=20000; %
n_wire=size(Cell_wire,1);

%
Vari_wiresegmentposition=cell(n_wire,1);
Vari_wiresegment=cell(n_wire,1);
%
tic;
for Count1=1:size(Cell_wire,1)

    if size(Cell_wire{Count1,2},1)>1
        if strcmpi(Cell_wire{Count1,3},'Cylindrical')
            Temp_wireposition=[Cell_wire{Count1,2}(:,1).*cos(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,1).*sin(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,3)];
        elseif strcmpi(Cell_wire{Count1,3},'Spherical')  %
            Temp_wireposition=repmat(Cell_wire{Count1,2}(:,2),[1,3]).*[sin(Cell_wire{Count1,2}(:,3)).*cos(Cell_wire{Count1,2}(:,2))...
                ,sin(Cell_wire{Count1,2}(:,3)).*sin(Cell_wire{Count1,2}(:,2)),cos(Cell_wire{Count1,2}(:,3))];
        elseif strcmpi(Cell_wire{Count1,3},'Cartesian')
            Temp_wireposition=Cell_wire{Count1,2};
        else
            error('wrong coordinates')
        end
        if vecnorm(Temp_wireposition(end,:)-Temp_wireposition(1,:))<vecnorm(Temp_wireposition(2,:)-Temp_wireposition(1,:))*1.5
            Temp_wireposition=[Temp_wireposition;Temp_wireposition(1,:)];
        end

        Temp_wiresegmentposition=(Temp_wireposition(2:end,:)+Temp_wireposition(1:end-1,:))/2;%
        Temp_wiresegment=Cell_wire{Count1,1}*(Temp_wireposition(2:end,:)-Temp_wireposition(1:end-1,:));%
        Vari_wiresegmentposition{Count1,1}=Temp_wiresegmentposition;%
        Vari_wiresegment{Count1,1}=Temp_wiresegment;%
    end

end
wiresegment=cell(0,1);
wiresegposition=cell(0,1);
for c1=1:length(Vari_wiresegment)
    tmp_l=size(Vari_wiresegment{c1,1},1);
    tmp_wiresegment=cell(0,1);
    tmp_wiresegposition=cell(0,1);
    if tmp_l>800*1.5
        id=round(linspace(0,tmp_l,round(tmp_l/800)+1));
        for c2=1:length(id)-1
            tmp_wiresegment{c2,1}=Vari_wiresegment{c1,1}(id(c2)+1:id(c2+1),:);
            tmp_wiresegposition{c2,1}=Vari_wiresegmentposition{c1,1}(id(c2)+1:id(c2+1),:);
        end
    else
        tmp_wiresegment=Vari_wiresegment{c1,1};
        tmp_wiresegposition=Vari_wiresegmentposition{c1,1};
    end

    wiresegment=[wiresegment;tmp_wiresegment];
    wiresegposition=[wiresegposition;tmp_wiresegposition];
end
Vari_wiresegment=wiresegment;
Vari_wiresegmentposition=wiresegposition;

%% 
for Count1=1:size(Cell_magneticfield,1)
    Cell_magneticfield{Count1,1}=[];
    Ctrl_fieldsequence=[0:Ctrl_maxfieldpoint:size(Cell_magneticfield{Count1,2},1),size(Cell_magneticfield{Count1,2},1)]; %判断点

    %
    for Count3=1:length(Ctrl_fieldsequence)-1



        if Ctrl_fieldsequence(Count3)~=Ctrl_fieldsequence(Count3+1)
            Vari_fieldposition=Cell_magneticfield{Count1,2}(Ctrl_fieldsequence(Count3)+1:Ctrl_fieldsequence(Count3+1),:); 
            %

            if strcmpi(Cell_magneticfield{Count1,3},'Cylindrical')
                Vari_fieldposition=[Vari_fieldposition(:,1).*cos(Vari_fieldposition(:,2)),Vari_fieldposition(:,1).*sin(Vari_fieldposition(:,2)),Vari_fieldposition(:,3)];
            elseif strcmpi(Cell_magneticfield{Count1,3},'Spherical')   
                Vari_fieldposition=repmat(Vari_fieldposition(:,1),[1,3]).*[sin(Vari_fieldposition(:,3)).*cos(Vari_fieldposition(:,2)),sin(Vari_fieldposition(:,3)).*sin(Vari_fieldposition(:,2)),cos(Vari_fieldposition(:,3))];
            elseif strcmpi(Cell_magneticfield{Count1,3},'Cartesian')

            else
                error('wrong coordinates')

            end
            

            field_all=zeros(size(Vari_fieldposition,1),3);
            parfor c1=1:length(Vari_wiresegment)


                Vari_r_rs={[],[],[]}; 
                for Count2=1:3
                    Vari_r_rs{Count2}=repmat(Vari_fieldposition(:,Count2),[1,size(Vari_wiresegmentposition{c1},1)])-repmat(Vari_wiresegmentposition{c1}(:,Count2),[1,size(Vari_fieldposition,1)])';
                end
                Vari_r3=(Vari_r_rs{1}.^2+Vari_r_rs{2}.^2+Vari_r_rs{3}.^2).^(3/2);
                for Count2=1:3 
                    Vari_r_rs{Count2}=Vari_r_rs{Count2}./Vari_r3;
                end
                %
                Temp_magneticfield=zeros(size(Vari_fieldposition,1),3);
                for Count2=fieldcomponent
                    i=mod(Count2,3)+1;
                    j=mod(Count2+1,3)+1;
                    Temp_magneticfield(:,Count2)=Vari_r_rs{j}*Vari_wiresegment{c1}(:,i)-Vari_r_rs{i}*Vari_wiresegment{c1}(:,j);
                end
                %
                field_all=field_all+Temp_magneticfield;
            end
        end
        Cell_magneticfield{Count1,1}(Ctrl_fieldsequence(Count3)+1:Ctrl_fieldsequence(Count3+1),:)=1e-7*field_all;
        %
    end

end

