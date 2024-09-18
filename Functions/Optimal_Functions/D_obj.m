function out=D_obj(in1,sz,term)
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

    out1=cell(size(in1,1),1);
    for c1=1:size(in1,1)
        if term==0

            id_tmp=find(in1(c1,3:4)>0);
            if ~isempty(id_tmp)
                temp_out=in1(c1,:);
                temp_out(2+id_tmp)=0;
                temp_out(5)=id_tmp;
                out1{c1,1}=[out1{c1,1};temp_out];
            end
            if in1(c1,5)>0
                temp_out=in1(c1,:);
                temp_out(6)=temp_out(5);
                temp_out(5)=0;
                out1{c1,1}=[out1{c1,1};temp_out];
            end
            
        elseif term==1
            if in1(c1,2)>0
            temp_out=repmat(in1(c1,:),[2,1]);
            temp_out(:,1)=-temp_out(:,1).*temp_out(:,2).*(sz.^2).';
            temp_out(:,2)=temp_out(:,2)+2;
            temp_out(1,9)=temp_out(1,9)+1;
            temp_out(2,10)=temp_out(2,10)+1;
            if in1(c1,11)~=0
                temp_out(1,12)=3;
                temp_out(2,12)=4;
            else
                temp_out(1,11)=3;
                temp_out(2,11)=4;
            end
            out1{c1,1}=[out1{c1,1};temp_out];
            end

            id_tmp=find(in1(c1,7:10)>0);
            if ~isempty(id_tmp)
                temp_out=repmat(in1(c1,:),[length(id_tmp),1]);
                for c2=1:length(id_tmp)
                temp_out(c2,1)=temp_out(c2,1)*temp_out(c2,6+id_tmp(c2));
            if in1(c1,11)~=0
                temp_out(c2,12)=id_tmp(c2);
            else
                temp_out(c2,11)=id_tmp(c2);
            end
                temp_out(c2,6+id_tmp(c2))=temp_out(c2,6+id_tmp(c2))-1;
                end
                out1{c1,1}=[out1{c1,1};temp_out];
            end
            if in1(c1,11)>0
                temp_out=in1(c1,:);
                temp_out(1,13)=in1(c1,11);
                temp_out(1,11)=0;
            out1{c1,1}=[out1{c1,1};temp_out];
            end
            
        end
    end
    out=cell2mat(out1);


