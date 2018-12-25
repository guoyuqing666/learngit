function ret = gsf(leftposFlame,rightposFlame,N,res)
%  leftposFlame is the left boundary position of flame layer
%  rightposFlame is the right boundary position of flame layer
%  N is the total points of the domain
%  res is a 1-D vector to be partly filtered
ret = zeros(N, 1);
theta=1.0;
for i=1:N
    if (i>=leftposFlame && i<=rightposFlame)
        ftotal=0.0;
        for j=i-50:i+50
            f(j)=1/(theta*(2*3.1416)^0.5)*exp(-(j-i)^2/(2*theta));
            ftotal=ftotal+f(j);
        end
        restotal=0.0;
        for j=i-50:i+50
            restotal=f(j)*res(j)+restotal;
        end
        ret(i)=restotal/ftotal;
    else
        ret(i)=res(i);    
    end
end
end

