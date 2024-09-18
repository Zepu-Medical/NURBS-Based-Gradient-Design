function [multi]=multiplicity(in)
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

    [in,~]=sort(in(:),"ascend");
    id=find((in(2:end)-in(1:end-1))~=0);
    id=[0;id;length(in)];
        multi=id(2:end)-id(1:end-1);
end