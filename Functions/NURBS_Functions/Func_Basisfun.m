function [bs,dbs,d2bs]= Func_Basisfun(xi_i, xi, p, Knot)
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

  xi=xi(:);

  Knot=Knot(:);

xi_i=xi_i(:);



    bs{1}=ones([length(xi),1]);
    dbs{1}=zeros([length(xi),1]);
    d2bs{1}=zeros([length(xi),1]);
    if any(xi<0) || any(xi>Knot(end))
        1;
    end
    if any(xi_i<0) || any(xi>Knot(end))
        1;
    end
for k=2:p+1
    
    try  
        d_n=reshape(Knot(xi_i+k-flip(0:k-2))-Knot(xi_i+1-flip(0:k-2)),[],length(0:k-2));
    catch
        errorvalue={xi_i,xi};
    end
    d_l=xi-reshape(Knot(xi_i+1-flip(0:k-2)),[],length(0:k-2));
    d_r=reshape(Knot(xi_i+k-flip(0:k-2)),[],length(0:k-2))-xi;
    bs{k}=zeros([length(xi),k]);
    bs{k}(:,1:end-1)=bs{k}(:,1:end-1)+(d_r)./(d_n).*bs{k-1};
    bs{k}(:,2:end)=bs{k}(:,2:end)+(d_l)./(d_n).*bs{k-1};
    
    dbs{k}=zeros([length(xi),k]);
    dbs{k}(:,1:end-1)=dbs{k}(:,1:end-1)-(k-1)*bs{k-1}./d_n;
    dbs{k}(:,2:end)=dbs{k}(:,2:end)+(k-1)*(bs{k-1})./d_n;

    d2bs{k}=zeros([length(xi),k]);

    d2bs{k}(:,1:end-1)=d2bs{k}(:,1:end-1)-(k-1)*(dbs{k-1})./d_n;
    d2bs{k}(:,2:end)=d2bs{k}(:,2:end)+(k-1)*dbs{k-1}./d_n;

end

end



