function [Basis_sin,Basis_cos]=Func_BasisSinu(k,l_z,Order_z)
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

q=1:Order_z;
p1=q.'*pi-k*l_z;
p2=q.'*pi+k*l_z;
f1=(exp(1j*p1)-1)./p1;
f2=(exp(-1j*p2)-1)./p2;
f1(p1==0)=1j;
f2(p2==0)=-1j;
Basis_sin=-l_z/(2*sqrt(2*pi))*(f1+f2);  
Basis_cos=-1j*l_z/(2*sqrt(2*pi))*(f1-f2);

end