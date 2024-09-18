function out=hessen_cat(mat_position,mat_cell)
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

    n=sqrt(size(mat_cell{1},1));
    i1=[];
    i2=[];
    val=[];
for c1=1:size(mat_position,1)
    tp=mat_position(c1,:)-[1,1];
    [i1_tmp,i2_tmp,val_tmp]=find(mat_cell{c1,1});
    [ix,iy]=ind2sub([n,n],i1_tmp);
    
    if tp(1)==tp(2)
        i1_tmp=sub2ind([3*n,3*n],tp(1)*n+ix,tp(2)*n+iy);

    else
        i1_tmp=sub2ind([3*n,3*n],[tp(1)*n+ix;tp(2)*n+iy],[tp(2)*n+iy;tp(1)*n+ix]);
        i2_tmp=[i2_tmp;i2_tmp];
        val_tmp=[val_tmp;val_tmp];
        
    end
    i1=[i1;i1_tmp];
    i2=[i2;i2_tmp];
    val=[val;val_tmp];
end
out=sparse(i1,i2,val,9*n^2,size(mat_cell{1},2));

