function [A,b] = FEM(mesh,f)
% FEM assembles the stiffness matrix and right-hand side for Poisson
% equation with the Finite Element Method. 
%
% SYNOPSIS: [A,b] = FEM(mesh,f)
%
% INPUT: mesh:	structure with information of the mesh:
%        mesh.verts:       matrix with coordinates of the vertices
%        mesh.elems:       matrix element connectivity
%        
%        f:     right-hand side of PDE, -\nabla u = f
%
% OUTPUT: A: stiffness matrix with entries 
%            \int_D \nabla\phi_i \cdot \nabla\phi_j
%         b: rhs with entries \int_D \phi_j f
%

% AUTHOR: Juan G. Calvo and collaborators, 2024
node = mesh.verts;
elem = reshape(cell2mat(mesh.elems),3,[])';
N = size(node,1); NT = size(elem,1);
ii = zeros(9*NT,1); jj = zeros(9*NT,1); sA = zeros(9*NT,1);
ve(:,:,3) = node(elem(:,2),:)-node(elem(:,1),:);
ve(:,:,1) = node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2) = node(elem(:,1),:)-node(elem(:,3),:);
area = 0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
index = 0;
for i = 1:3
    for j=1:3
        ii(index+1:index+NT) = elem(:,i);
        jj(index+1:index+NT) = elem(:,j);
        sA(index+1:index+NT) = dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        index = index + NT;
    end
end
A= sparse(ii,jj,sA,N,N);
% right hand side with one-point quadrature rule
barycenter = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
bt = f(barycenter(:,1),barycenter(:,2))/3;
bt = bt.*repmat(area,1,3);
b  = accumarray(elem(:),bt(:),[N 1]);
end