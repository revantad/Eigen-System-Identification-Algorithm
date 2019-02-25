function [system,A,B,C,D] = eigenSys(yMat,uMat,numMP,numH,dt)
%function inputs=(size of Input Vector, Output Vector,numMP, 
%Size of Hankel Matrix, Sample Time)
%% fOKID
p=2;
[m,N]=size(yMat);
[r,N]=size(uMat);
V=zeros(p*(r+m)+r,N);
V(1:r,1:N)=uMat;

for i=1:N
    V(r*i+1+(i-1)*m:(i+1)*r+(i-1)*m,(i+1):N)=uMat(:,1:N-i);
    V((i+1)*r+(i-1)*m+1:(i+1)*r+i*m,(i+1):N)=yMat(:,1:N-i);
end
V=V([1:p*(r+m)+r],:);
yBar=yMat*pinv(V);

%% fSYSMP

y(1:m,1:r)=yBar(1:m,1:r);
Y1=zeros(p*numMP,r);
Y2=zeros(p*numMP,m);
for k=1:p
    Ybar((k-1)*m+1:k*m,:)=yBar(1:m,(k-1)*m+(k*r)+1:k*m+(k+1)*r);
    Y1((k-1)*m+1:k*m,1:r)=Ybar((k-1)*m+1:k*m,1:r);
    Y2((k-1)*m+1:k*m,1:m)=Ybar((k-1)*m+1:k*m,1+r:r+m);
end
y2=zeros(m,r);
for k=1:numMP
    for j=1:k 
        y2=y2+Y2((j-1)*m+1:j*m,1:m)*y((k-j)*m+1:(k-j+1)*m,1:r);
    end
    y(k*m+1:(k+1)*m,1:r)=Y1((k-1)*m+1:k*m,1:r)+y2;
    y2=zeros(m,r);
end
Y=y(m+1:end,:);

%% System Identification
for k=1:numH
    for j=1:numH
        H0(k,j)=Y(j+k-1,1);
        H1(k,j)=Y(j+k,1);
    end
end
[R,E,S]=svd(H0);
p=rank(E);
En=E(1:p,1:p);
Sn=S(:,1:p);
Rn=R(:,1:p);
pHat=Rn*sqrt(En); % Solving for PHat
qHat=sqrt(En)*Sn'; % Solving for QHat
pseudo=pinv(pHat); % Pseudo Inverse of PHat
qseudo=pinv(qHat); % Pseudo Inverse of QHat
B=qHat(:,1:r);
C=pHat(1:m,:);
A=pseudo*H1*qseudo;
D=yBar(1:m,1:r);
sys=ss(A,B,C,D);

[num,den]=ss2tf(A,B,C,D);
system=tf(num,den,dt);
end