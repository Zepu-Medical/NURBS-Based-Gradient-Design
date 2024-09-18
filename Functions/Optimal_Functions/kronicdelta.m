function out=kronicdelta(input1)
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

sz1=size(input1);
ind=sub2ind([sz1(1),sz1(1)],1:sz1(1),1:sz1(1));

if ~issparse(input1)
    out=zeros(sz1(1)^2,sz1(2));
    for cl1=1:sz1(2)
       out(ind,cl1)=input1(:,cl1); 
    end

else
    [id1,id2,val]=find(input1);
    out=sparse(ind(id1),id2,val,sz1(1)^2,sz1(2));
end


