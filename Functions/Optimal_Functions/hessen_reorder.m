function mat_out = hessen_reorder(mat)
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

id1=[];
id2=[];
val=[];
for c1=1:size(mat,1)
    for c2=c1:size(mat,2)
        if nnz(mat{c1,c2})>0
            sz1=size(mat{c1,c2});
            sz1(1)=sqrt(sz1(1));
            [id1_tmp,id2_tmp,val_tmp]=find(mat{c1,c2});
            [idx,idy]=ind2sub([sz1(1),sz1(1)],id1_tmp);

            if c1==c2
            id1=[id1;sub2ind([sz1(1)*3,sz1(1)*3],idx+(c1-1)*sz1(1),idy+(c2-1)*sz1(1))];
            id2=[id2;id2_tmp];
            val=[val;val_tmp];
            
            else
            id1=[id1;sub2ind([sz1(1)*3,sz1(1)*3],idx+(c1-1)*sz1(1),idy+(c2-1)*sz1(1))];
            id1=[id1;sub2ind([sz1(1)*3,sz1(1)*3],idx+(c2-1)*sz1(1),idy+(c1-1)*sz1(1))];
            id2=[id2;id2_tmp;id2_tmp];
            val=[val;val_tmp;val_tmp];
            end
        end
    end
end

mat_out=sparse(id1,id2,val,sz1(1).^2*9,sz1(2));
end
    
