function u_span=Func_ParaSpan(u,knots)
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

knots=sort(knots,'ascend');
ind=find(knots==knots(end));
knots(ind)=[];
knots=knots(:);
t1=[knots(2:end);u(:)];
index=1:length(t1);    
[a,b]=sort(t1);
    index=index(b);
    n=find(index<numel(knots));
    u_span=zeros(size(u));
    for c1=1:length(n)-1
    u_span(index(n(c1)+1:n(c1+1)-1)-length(knots)+1)=c1;
    end
    if n(c1+1)~=length(b)
    u_span(index(n(end)+1:end)-length(knots)+1)=length(n);
    end

    
end