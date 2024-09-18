function NURBS_curve = coeff2NURBS(coeff, NURBS_curve)
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

id=[];

for c1=1:length(NURBS_curve)
    for c2=1:length(NURBS_curve(c1).Knot)
        id(end+1)=length(NURBS_curve(c1).Weight{c2});
    end
end

coeff=mat2cell(coeff,id,[2,1]);
c3=1;
for c1=1:length(NURBS_curve)
    for c2=1:length(NURBS_curve(c1).Knot)

        NURBS_curve(c1).ContrPoint{c2}(:,1:2)=coeff{c3,1};
        NURBS_curve(c1).Weight{c2}=coeff{c3,2};
        c3=c3+1;
    end
end
end

