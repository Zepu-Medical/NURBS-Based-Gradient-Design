function [Cell_gradientfield] = Func_GradientField(Cell_gradientfield,Cell_wire,fieldcomponent)
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


maxfieldpoint=20000; 
n_wire=size(Cell_wire,1);

%
Vari_wiresegmentposition=cell(n_wire,1);
Vari_wiresegment=cell(n_wire,1);
%
tic;
for Count1=1:size(Cell_wire,1)
    %
    if size(Cell_wire{Count1,2},1)>1

            if strcmpi(Cell_wire{Count1,3},'Cylindrical')
                Temp_wireposition=[Cell_wire{Count1,2}(:,1).*cos(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,1).*sin(Cell_wire{Count1,2}(:,2)),Cell_wire{Count1,2}(:,3)];
            elseif strcmpi(Cell_wire{Count1,3},'Spherical')  %
                Temp_wireposition=repmat(Cell_wire{Count1,2}(:,2),[1,3]).*[sin(Cell_wire{Count1,2}(:,3)).*cos(Cell_wire{Count1,2}(:,2))...
                    ,sin(Cell_wire{Count1,2}(:,3)).*sin(Cell_wire{Count1,2}(:,2)),cos(Cell_wire{Count1,2}(:,3))];
            elseif strcmpi(Cell_wire{Count1,3},'Cartesian')
                Temp_wireposition=Cell_wire{Count1,2};
            else
                error('导线坐标系错误');
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
for Count1=1:size(Cell_gradientfield,1)
    Ctrl_fieldsequence=[0:maxfieldpoint:size(Cell_gradientfield{Count1,2},1),size(Cell_gradientfield{Count1,2},1)]; %

    %
    for Count3=1:length(Ctrl_fieldsequence)-1
        if Count3>1
        (Count3-1)*maxfieldpoint
        toc;
        end

        if Ctrl_fieldsequence(Count3)~=Ctrl_fieldsequence(Count3+1)
        Vari_fieldposition=Cell_gradientfield{Count1,2}(Ctrl_fieldsequence(Count3)+1:Ctrl_fieldsequence(Count3+1),:); %
        %
        
             if strcmpi(Cell_gradientfield{Count1,3},'Cylindrical')
                Vari_fieldposition=[Vari_fieldposition(:,1).*cos(Vari_fieldposition(:,2)),Vari_fieldposition(:,1).*sin(Vari_fieldposition(:,2)),Vari_fieldposition(:,3)];
             elseif strcmpi(Cell_gradientfield{Count1,3},'Spherical')  %
                Vari_fieldposition=repmat(Vari_fieldposition(:,1),[1,3]).*[sin(Vari_fieldposition(:,3)).*cos(Vari_fieldposition(:,2)),sin(Vari_fieldposition(:,3)).*sin(Vari_fieldposition(:,2)),cos(Vari_fieldposition(:,3))];
             elseif strcmpi(Cell_gradientfield{Count1,3},'Cartesian')

            else
                error('场点坐标系错误');

        end


            field_all=zeros(size(Vari_fieldposition,1),3);
    parfor c1=1:length(Vari_wiresegment)
        DR={[],[],[]}; 
        for Count2=1:3
            DR{Count2}=repmat(Vari_fieldposition(:,Count2),[1,size(Vari_wiresegmentposition{c1},1)])-repmat(Vari_wiresegmentposition{c1}(:,Count2),[1,size(Vari_fieldposition,1)])';
        end
        iDR25=(DR{1}.^2+DR{2}.^2+DR{3}.^2).^(-5/4);
        for Count2=1:3 
            DR{Count2}=DR{Count2}.*iDR25;
        end
        %
        Temp_gradientfield=zeros(size(Vari_fieldposition,1),3);

        fieldcomponent_par=sort(fieldcomponent);
        if length(fieldcomponent_par)>1
            fieldcomponent_par=[fieldcomponent_par(fieldcomponent_par(2:end)-fieldcomponent_par(1:end-1)~=0),fieldcomponent_par(end)];
        end
        
        for Count2=fieldcomponent_par
            if Count2==1
                Temp_gradientfield(:,1)=-3*DR{1}.*DR{2}*Vari_wiresegment{c1}(:,1)+(2*DR{1}.^2-DR{2}.^2-DR{3}.^2)*Vari_wiresegment{c1}(:,2);
            elseif Count2==2
                Temp_gradientfield(:,2)=(DR{1}.^2-2*DR{2}.^2+DR{3}.^2)*Vari_wiresegment{c1}(:,1)+3*DR{1}.*DR{2}*Vari_wiresegment{c1}(:,2);
            elseif Count2==3
                Temp_gradientfield(:,3)=-3*DR{2}.*DR{3}*Vari_wiresegment{c1}(:,1)+3*DR{1}.*DR{3}*Vari_wiresegment{c1}(:,2);
            end
        end
        %
        field_all=field_all+Temp_gradientfield;
    end
    Cell_gradientfield{Count1,1}(Ctrl_fieldsequence(Count3)+1:Ctrl_fieldsequence(Count3+1),:)=1e-7*field_all;
    end
    end
end

