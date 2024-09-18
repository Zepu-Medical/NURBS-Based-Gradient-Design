function [abs_length,gf]=nurbs_length(xi_e,xi_s,Knot,coeff,val0)
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

if ~exist("val0","var")
    val0=0;
end
dx=[];
parfor c1=1:length(xi_e)
    dx(c1,1)=integral(@(xi)Func_dcurve(xi,Knot,coeff),xi_s(c1),xi_e(c1),ArrayValued=true)-val0;
end
    
    abs_length=abs(dx);

    gf=sign(dx).*Func_dcurve(xi_e,Knot,coeff);



