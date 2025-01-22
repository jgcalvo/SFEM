function [A,b] = ESFEM(mesh,f)
% SFEM assembles the stiffness matrix and right-hand side for Poisson
% equation with the edge-smoothed Finite Element Method. 
%
% SYNOPSIS: [A,b] = SFEM(mesh,f)
%
% INPUT: mesh:	structure with information of the mesh:
%        mesh.verts:       matrix with coordinates of the vertices
%        mesh.elems:       matrix element connectivity
%        mesh.edges:       matrix with edge connectivity
%        mesh.edges2elems: matrix with elements shared by an edge
%        
%        f:     right-hand side of PDE, -\nabla u = f
%
% OUTPUT: A: stiffness matrix with entries for the edge-smoothed FEM
%         b: rhs with entries \int_D \phi_j f
%

% AUTHOR: Juan G. Calvo and collaborators, 2024

% read mesh vertices and number of vertices
verts = mesh.verts;
NN    = size(verts,1);
% change element connectivity to a matrix
elems = reshape(cell2mat(mesh.elems),3,size(mesh.elems,1))';
% compute gradients of basis functions for all elements
ve(:,:,3) = verts(elems(:,2),:)-verts(elems(:,1),:);
ve(:,:,1) = verts(elems(:,3),:)-verts(elems(:,2),:);
ve(:,:,2) = verts(elems(:,1),:)-verts(elems(:,3),:);
area      = 0.5*(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));
gradients = [-ve(:,2,:)./(2*area), ve(:,1,:)./(2*area)];
% allocate matrices to create stiffness matrix as A = sparse(ii,jj,ss)
% 16 dof for internal edges + 9 dof for boundary edges
numStore = 16*(size(mesh.edges,1)-numel(mesh.bdEdge))+9*numel(mesh.bdEdge);
ii = nan(numStore,1);
jj = nan(numStore,1);
ss = nan(numStore,1);
pt = 0;
% loop over all edges of the mesh
for edgeID = 1:size(mesh.edges,1)
    edge   = mesh.edges(edgeID,:);
    neighs = mesh.edges2elems{edgeID};
    if(numel(neighs)==2) % internal edge with two neighbors
        t1 = neighs(1);                          % t1 - element 1
        t2 = neighs(2);                          % t2 - element 2
        i1 = (elems(t1,:)~=edge(1))&(elems(t1,:)~=edge(2));
        i2 = (elems(t2,:)~=edge(1))&(elems(t2,:)~=edge(2));
        nn = [edge, elems(t1,i1), elems(t2,i2)]; % four vertices
        w1 = area(t1)/(area(t1)+area(t2));       % weights for 
        w2 = area(t2)/(area(t1)+area(t2));       % average gradient
        Aloc = nan(4);                           % local matrix
        for i=1:4                                % index phi_i
            if(i==1 || i==2)                     % vertex on edge
                grad1 = gradients(t1,:,elems(t1,:)==nn(i));
                grad2 = gradients(t2,:,elems(t2,:)==nn(i));
            elseif(i==3)                         % vertex on t1
                grad1 = gradients(t1,:,i1);
                grad2 = [0 0];
            else                                 % vertex on t2
                grad1 = [0 0];
                grad2 = gradients(t2,:,i2);
            end
            gradi = (w1*grad1+w2*grad2);         % weighted average phi_i

            for j=i:4                            % index phi_j
                if(j==1 || j==2)                 % vertex on edge
                    grad1 = gradients(t1,:,elems(t1,:)==nn(j));
                    grad2 = gradients(t2,:,elems(t2,:)==nn(j));
                elseif(j==3)                     % vertex on t1
                    grad1 = gradients(t1,:,elems(t1,:)==nn(j));
                    grad2 = [0 0];
                else                             % vertex on t2
                    grad1 = [0 0];
                    grad2 = gradients(t2,:,elems(t2,:)==nn(j));
                end
                gradj = (w1*grad1+w2*grad2);     % weighted average phi_i
                % entry for the sfiffness matrix
                valij = (area(t1)+area(t2))/3*(gradi*gradj');
                % store entries for local symmetric matrix
                Aloc(i,j) = valij;
                Aloc(j,i) = valij;
            end
        end
        xx = repmat(nn,4,1);
        yy = xx';
        ii(pt+1:pt+16) = xx(:);   % 
        jj(pt+1:pt+16) = yy(:);   %
        ss(pt+1:pt+16) = Aloc(:); % 
        pt = pt+16;               % update pointer 
    else                     % boundary edge with one neighbor
        t1 = neighs(1);
        nn = elems(neighs(1),:);
        Aloc = nan(3);
        for i=1:3
            gradi = gradients(t1,:,i);
            for j=1:3
                gradj = gradients(t1,:,j);
                Aloc(i,j) = area(t1)/3*(gradi*gradj');
            end
        end
        xx = repmat(nn,3,1);
        yy = xx';
        ii(pt+1:pt+9) = xx(:);    %
        jj(pt+1:pt+9) = yy(:);    %
        ss(pt+1:pt+9) = Aloc(:);  %
        pt = pt + 9;              % update pointer 
    end
end
% construct sparse matrix
A = sparse(ii,jj,ss,NN,NN);
% rhs is similar as FEM
barycntr = (verts(elems(:,1),:)+verts(elems(:,2),:)+verts(elems(:,3),:))/3;
bt = f(barycntr(:,1),barycntr(:,2))/3;
bt = bt.*repmat(area,1,3);
b = accumarray(elems(:),bt(:),[NN 1]);
end