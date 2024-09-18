function point_unfold = UnfoldPoint(points_fold, WirePara)
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

    point_unfold=zeros(4*size(points_fold,1),3); 
    point_unfold(:,1)=kron([1;-1;-1;1],points_fold(:,1));
    point_unfold(:,2)=kron([1;1;-1;-1],points_fold(:,2));
    point_unfold(:,3)=kron([1;1;1;1],points_fold(:,3));
    if strcmpi(WirePara.fold_condition,'8-fold')
    point_unfold=[point_unfold;point_unfold.*[1,1,-1]];
    end
end