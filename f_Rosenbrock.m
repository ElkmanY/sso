function [ Y ] = f_Rosenbrock( X )
%   optimization function named Rosenbrock
%   Y=sum( 100*(X_bk-X_ft.^2).^2+(X_ft-1).^2, 2);   X is a NxD matrix
X_bk=X;
X_bk(:,1)=[];
X_ft=X;
X_ft(:,size(X,2))=[];
Y=sum( 100*(X_bk-X_ft.^2).^2+(X_ft-1).^2, 2);
end
%%  draw a figure 2-D problems
% x1=-30:0.1:30;
% x2=x1;
% [X1, X2]=meshgrid(x1,x2);
% Y=100*(X2-X1.^2).^2 + (X1-1).^2;
% mesh(X1,X2,Y);


% X=[1 2 3;3 4 5;5 6 7;7 8 9];  % test data 