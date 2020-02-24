%Assignment 2 for ELEC4700
%Philippe Masson

%Q1 A -------------------
%clear;
close all;

Width = 50;
Length = Width * 3/2;
numX = Length;
numY = Width;
delta = 1;
V0 = 10;

GMatrix = zeros(numX*numY,numX*numY);

for i = 1:numX
    for j = 1:numY
        n = j+(i-1)*numY;
        nxm = j+(i-2)*numY;
        nxp = j + i*numY;
        
        if(i == 1) %Left Boundary
            GMatrix(n,n) = V0;
        elseif(i == numX) %Right Boundary
            GMatrix(n,n) = 1; %using 1 here since 0 is bugged for eigs
        else            
            GMatrix(n,n) = -2;
            GMatrix(nxm,n) = 1;
            GMatrix(nxp,n) = 1;
        end      
        
    end    
end

[E,D] = eigs(GMatrix,2,'sm');

jC = 1;
iC = 1;

mapMatrix = zeros(numX,numY);
for count = 1:numX*numY
        mapMatrix(iC,jC) = E(count,1);

        if (jC == numY)
            jC = 1;
            iC = iC + 1;        
        else
            jC = jC + 1;
        end
end
fig = figure(2);
plot(mapMatrix(:,1));
title('Q1.A) 2D Plot of V(x)');
xlabel('X');
ylabel('V(x)');
%surf(mapMatrix*D(1,1));
saveas(fig, ['Q1A', '.png']);


%Q1 B -------------------
clear;

Width = 50;
Length = Width * 3/2;
numX = Length;
numY = Width;
delta = 1;
V0 = 10;

GMatrix = zeros(numX*numY,numX*numY);

for i = 1:numX
    for j = 1:numY
        n = j+(i-1)*numY;
        nxm = j+(i-2)*numY;
        nxp = j + i*numY;
        nym = j-1+(i-1)*numY;
        nyp = j+1 + (i-1)*numY;
        
        if(i == 1) %Left Boundary
            GMatrix(n,n) = V0;
        elseif(i == numX) %Right Boundary
            GMatrix(n,n) = V0;
        elseif(j == 1)
            GMatrix(n,n) = 1; %using 1 here since 0 is bugged for eigs
        elseif (j == numY)
            GMatrix(n,n) = 1; %using 1 here since 0 is bugged for eigs
        else            
            GMatrix(n,n) = -4;
            GMatrix(nxm,n) = 1;
            GMatrix(nxp,n) = 1;
            GMatrix(n,nym) = 1;
            GMatrix(n,nyp) = 1;
        end      
        
    end    
end

[E,D] = eigs(GMatrix,2,'sm');
%[E,D] = eig(GMatrix);
% [Ds,Pr] = sort(diag(D));
% D = D(Pr,Pr);
% E = E(:,Pr);


jC = 1;
iC = 1;

mapMatrix = zeros(numX,numY);
for count = 1:numX*numY
        mapMatrix(iC,jC) = E(count,1);

        if (jC == numY)
            jC = 1;
            iC = iC + 1;        
        else
            jC = jC + 1;
        end
end

fig = figure(3);
surf(mapMatrix*D(1,1));
title('Q1.B) Surface Plot of V(x,y)');
xlabel('X');
ylabel('Y');
zlabel('V(x,y)');
saveas(fig, ['Q1B', '.png']);

%Q2 -------------------
clear

Width = 50;
Length = Width * 3/2;
numX = Length;
numY = Width;
V0 = 10;

GMatrix = zeros(numX*numY,numX*numY);
SigMatrix = zeros(numX,numY);

for i = 1:numX
    for j = 1:numY
        n = j+(i-1)*numY;
        nxm = j+(i-2)*numY;
        nxp = j + i*numY;
        nym = j-1+(i-1)*numY;
        nyp = j+1 + (i-1)*numY;
        
        %Setup Sigma
        Sigma = 1;
        if(i >=25 && i <= 45 && (j<15 || j>numY-15))
            Sigma = 10e-2;
        end
        SigMatrix(i,j) = Sigma;
        
        if(i == 1) %Left Boundary
            GMatrix(n,n) = V0;
        elseif(i == numX) %Right Boundary
            GMatrix(n,n) = 1;
        else            
            GMatrix(n,n) = -4*Sigma;
            GMatrix(nxm,n) = 1*Sigma;
            GMatrix(nxp,n) = 1*Sigma;
            GMatrix(n,nym) = 1*Sigma;
            GMatrix(n,nyp) = 1*Sigma;
        end

    end    
end

[V,D] = eigs(GMatrix,1,'sm');
jC = 1;
iC = 1;
mapMatrix = zeros(numX,numY);
for count = 1:numX*numY
        mapMatrix(iC,jC) = V(count,1);

        if (jC == numY)
            jC = 1;
            iC = iC + 1;        
        else
            jC = jC + 1;
        end
end

V=mapMatrix*D(1,1);
%V=mapMatrix;
I=V./SigMatrix; %get the current from V
q = 1.6e-19; %charge
E=V*q; %Energy from V


fig = figure(4);
surf(SigMatrix);
title('2 - sigma(x,y)');
xlabel('X');
ylabel('Y');
zlabel('sigma(x,y)');
saveas(fig, ['Q2Sig', '.png']);

fig = figure(5);
surf(V);
title('2 - V(x,y)');
xlabel('X');
ylabel('Y');
zlabel('V(x,y)');
saveas(fig, ['Q2V', '.png']);

fig = figure(6);
surf(I);
title('2 - J(x,y)');
xlabel('X');
ylabel('Y');
zlabel('E(x,y)');
saveas(fig, ['Q2J', '.png']);

fig = figure(7);
surf(E);
title('2 - E(x,y)');
xlabel('X');
ylabel('Y');
zlabel('J(x,y)');
saveas(fig, ['Q2E', '.png']);

