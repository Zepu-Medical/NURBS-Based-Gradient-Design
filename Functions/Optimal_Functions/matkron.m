function out=matkron(in1,in2)
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

if ~issparse(in1) && ~issparse(in2)
    out=zeros(size(in2,1)*size(in1,1),size(in1,2));
    for cl1=1:size(in1,2)
        out(:,cl1)=kron(in2(:,cl1),in1(:,cl1)); 
    end
else 
    A=cell(size(in1,2),1);
    for cl1=1:size(in1,2)
        A{cl1}=kron(in2(:,cl1),in1(:,cl1)); 
    end
        out=cat(2,A{:});
end

end