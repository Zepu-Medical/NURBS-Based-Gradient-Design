function Plot_Wire3D(Cell_wire)
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
for c1=1:size(Cell_wire,1)
    hold on
    x=Cell_wire{c1,2}(:,1).*cos(Cell_wire{c1,2}(:,2));
    y=Cell_wire{c1,2}(:,1).*sin(Cell_wire{c1,2}(:,2));
    load ColorSet.mat
    if Cell_wire{c1,1}>0
    plot3(x,y,Cell_wire{c1,2}(:,3),'LineWidth',0.75,'Color',rgbset(1,:));
    elseif Cell_wire{c1,1}<0
    plot3(x,y,Cell_wire{c1,2}(:,3),'LineWidth',0.75,'Color',rgbset(2,:));
    else
        plot3(x,y,Cell_wire{c1,2}(:,3),'LineWidth',0.75,'Color',rgbset(6,:));
    end
end
xlabel('y/m')
ylabel('x/m')
zlabel('z/m')
view([-8,6])