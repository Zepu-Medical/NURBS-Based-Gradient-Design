function [out] = term2val(sz0, eq, Vn, C0, C1, C0dpdq, C1dpdq, C1dpC1dq, C0dpC1dq, C0dp, C1dp)
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

out=sparse(sz0(1),sz0(2));
for c2=1:size(eq,1)
    check1=0;
    tp=eq(c2,1)*Vn{(eq(c2,2)+1)/2};
    id=find(eq(c2,3:4)>0);
    for c3=1:length(id)
        tp=tp.*C0{id(c3)}.^eq(c2,2+id(c3));

    end
    id=find(eq(c2,7:10)>0);
    for c3=1:length(id)
        tp=tp.*C1{id(c3)}.^eq(c2,6+id(c3));

    end
    
    if eq(c2,6)>0
        tp=tp.*C0dpdq{eq(c2,6)};
        check1=check1+2;
    elseif eq(c2,13)>0
        
        tp=tp.*C1dpdq{eq(c2,13)};
        check1=check1+2;

    elseif eq(c2,11)>0&&eq(c2,12)>0
        tp=tp.*C1dpC1dq{eq(c2,11),eq(c2,12)};
        check1=check1+2;

    elseif eq(c2,5)>0&&eq(c2,11)>0
        tp=tp.*C0dpC1dq{eq(c2,5),eq(c2,11)};
        check1=check1+2;

    elseif eq(c2,5)>0
        tp=tp.*C0dp{eq(c2,5)};
        check1=check1+1;

    elseif eq(c2,11)>0    
        tp=tp.*C1dp{eq(c2,11)};
        check1=check1+1;
    end
    out=out+tp;
    if check1>2
    error('!')
    end
end


end