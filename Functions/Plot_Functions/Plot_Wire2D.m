function CoilObj=Plot_Wire2D(Coil2D,WireIdcontrol)
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
for c0=1:length(Coil2D)

for c1=1:length(Coil2D(c0).Wire2D)

    Wire_temp=Coil2D(c0).Wire2D{c1};

    %subplot(1,3,c0)
    hold on
    CoilObj(c0).Wire(c1)=plot(Wire_temp(:,1),Wire_temp(:,2),'--');

    if exist('WireIdcontrol','var')
    if strcmpi(WireIdcontrol,'wireId on')
    text(Wire_temp(3,1),Wire_temp(3,2),num2str(Coil2D(c0).WireId(c1)));
    end
    end
end
end
end