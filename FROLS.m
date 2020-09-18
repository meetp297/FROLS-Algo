clear all
close all
clc
rho=0.05
ny=2                                         % nummber of y terms
nu=2                                         % nummber of u terms
ne=0                                         % nummber of e terms
n=ny+nu+ne
x=3
M= factorial(n+x)/(factorial(n)*factorial(x))
mu=0;
sigma=0.2;
e=sigma*randn(200,1)+mu;
u=-1+2*rand(200,1);
y=[0;0.5]

for k=3:200
    y(k)=(-0.605*y(k-1))-(0.163*y(k-2)^2)+(0.588*u(k-1))-(0.240*u(k-2))+e(k);
end
ly=200

b=2;
C=zeros(M,n);   % all terms c(35,4)
% columns of c is u(k-1),u(k-2),y(k-1),y(k-2)
% each row represent one term 
% elements of row represent the powers of u(k-1),u(k-2),y(k-1),y(k-2) respectively in each terms
for p=1:x
    if p==1
      for i_1=1:n
          C(b,i_1)=C(b,i_1)+1;  
          b=b+1;
      end
    end
    if p==2
        for i_1=1:n
            for i_2=i_1:n
                C(b,i_1)=C(b,i_1)+1;
                C(b,i_2)=C(b,i_2)+1;
                b=b+1;
            end
        end
    end
    
    if p==3
        for i_1=1:n
            for i_2=i_1:n
                for i_3=i_2:n
                    C(b,i_1)=C(b,i_1)+1;
                    C(b,i_2)=C(b,i_2)+1;
                    C(b,i_3)=C(b,i_3)+1;  
                     b=b+1;
                end
            end
        end
    end     
   
end
C 
sizeC=size(C)

% creating D

nm=max(nu,ny);
D=ones(ly-nm,M);
%size(D)
for i=nm+1:ly   %iteraing y
    for j=1:M   % iterating C
       k=1;
       for l=1:nu
       D(i-nm,j)=(D(i-nm,j))*(u(i-l)^(C(j,k)));
       k=k+1;
       end
       for m=1:ny
       D(i-nm,j)=(D(i-nm,j))*(y(i-m)^(C(j,k)));
       k=k+1;
       end
    end
end
size(D)





Y=y(3:200,:)
sig=Y'*Y;

for m=1:35
    q(:,m)=D(:,m);
    g(m,:)=(Y'*q(:,m))/(q(:,m)'*q(:,m));
    ERR(m,:)=(g(m,:))^2*q(:,m)'*q(:,m)/sig;
end
[err(1,:),l_1 ]=max(ERR);
Q(:,1)=D(:,l_1);
a11=1;



t=1;
l(1)=l_1;
for o=1:34
    t=t+1;
    ERR=zeros(34,1);
    for m=1:35
        
        s=zeros(198,1);
        A=ismember(m,l);
        if A==0
            for r=1:size(Q,2)
            s=s+(((D(:,m)'*Q(:,r))/(Q(:,r)'*Q(:,r)))*Q(:,r));
            end
            q(:,m)=D(:,m)-s;
            g(m,:)=(Y'*q(:,m))/(q(:,m)'*q(:,m));
            ERR(m,:)=(g(m,:)^2)*(q(:,m)'*q(:,m))/sig;
        else
        end
    end
        [err(t,:),l_2]=max(ERR);
        l(t)=l_2;               % l=unknown model terms
        Q(:,t)=q(:,l_2);
        ESR=1-sum(err);
        if ESR<=rho
            break
        end
    end
    

     [s1,s2]=size(l)
   for f=1:s2
       X(:,f)=D(:,l(f));
   end
   w=inv(X'*X)*X'*Y;    % w=unknown model parameters
