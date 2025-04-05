function [AVEM,ASSEVEM,b,out] = SSEVEM(mesh,rhs,varargin)
if(size(varargin,2)==1)
    rho = varargin{1};
else
    rho = ones(size(mesh.elems,1),1);
end


NN = size(mesh.verts,1);
NT = size(mesh.elems,1);
%
[AVEM,b,PP,stab,area,diam] = VEM(mesh,rhs,rho);
out.proj = PP;
out.mesh = mesh;
out.polys = [0 0; 1 0; 0 1];


nmax = max(cellfun(@numel,mesh.elems));
numStore = nmax*(nmax-1)*NT;
ii = nan(numStore,1);
jj = nan(numStore,1);
ss = nan(numStore,1);
pt = 0;
%
for elemID = 1:NT
    nodes  = mesh.elems{elemID};       % nodes of element
    edges  = mesh.elem2edge{elemID};   % edges of element
    neighs = mesh.edge2elem(edges,:);  % info of neighbors
    % define local nodes (with at most three triangles)
    nodeAdd = [];
    for edg = 1:numel(edges)
        if(neighs(edg,1)~=neighs(edg,2))
            if(neighs(edg,1)~=elemID)
                nodeAdd = [nodeAdd; setdiff(mesh.elems{neighs(edg,1)},nodes)];
            else
                nodeAdd = [nodeAdd; setdiff(mesh.elems{neighs(edg,2)},nodes)];
            end
        end
    end
    nodes = unique([nodes; nodeAdd]');
    % local matrix
    Aloc   = nan(numel(nodes)); % local matrix 3x3 to 6x6
    % loop for all pair of basis functions
    for i = 1:numel(nodes) % order: nodes and external nodes as edges order
        ni = nodes(i);
        for j = i:numel(nodes)
            nj = nodes(j);
            gradEdgi = nan(numel(edges),2);
            gradEdgj = nan(numel(edges),2);
            for edg = 1:numel(edges)               
                % find triangles that share edg
                if(neighs(edg,1)~=neighs(edg,2)) % two neighborhs
                    t1 = neighs(edg,1);
                    nodest1 = mesh.elems{t1};
                    t2 = neighs(edg,2);
                    nodest2 = mesh.elems{t2};
                    % phi_i with t1 and t2
                    ind = nodest1 == ni;
                    if(sum(ind)>0)
                        gradi1 = PP{t1}(2:3,ind)'/diam(t1);
                    else
                        gradi1 = [0 0];
                    end
                    ind = nodest2 == ni;
                    if(sum(ind)>0)
                        gradi2 = PP{t2}(2:3,ind)'/diam(t2);
                    else
                        gradi2 = [0 0];
                    end
                    % phi_j with t1 and t2
                    ind = nodest1 == nj;
                    if(sum(ind)>0)
                        gradj1 = PP{t1}(2:3,ind)'/diam(t1);
                    else
                        gradj1 = [0 0];
                    end
                    ind = nodest2 == nj;
                    if(sum(ind)>0)
                        gradj2 = PP{t2}(2:3,ind)'/diam(t2);
                    else
                        gradj2 = [0 0];
                    end
                    % average
                    w1 = area(t1)/(area(t1)+area(t2));
                    w2 = area(t2)/(area(t1)+area(t2));
                    gradEdgi(edg,:) = w1*gradi1+w2*gradi2;
                    gradEdgj(edg,:) = w1*gradj1+w2*gradj2;
                else                    % one neighbor
                    t1 = neighs(edg,1);
                    nodest1 = mesh.elems{t1};
                    % phi_i with t1
                    ind = nodest1 == ni;
                    if(sum(ind)>0)
                        gradi1 = PP{t1}(2:3,ind)'/diam(t1);
                    else
                        gradi1 = [0 0];
                    end
                    % phi_j with t1
                    ind = nodest1 == nj;
                    if(sum(ind)>0)
                        gradj1 = PP{t1}(2:3,ind)'/diam(t1);
                    else
                        gradj1 = [0 0];
                    end
                    % average
                    gradEdgi(edg,:) = gradi1;
                    gradEdgj(edg,:) = gradj1;
                end
            end
            % averages for each pair of consecutives edges
            gradEdgi = (gradEdgi([2:end 1],:)+gradEdgi)/2;
            gradEdgj = (gradEdgj([2:end 1],:)+gradEdgj)/2;
            % approx for integral
            valij = area(elemID)*sum(sum(gradEdgi.*gradEdgj,2))/size(gradEdgi,1);
            Aloc(i,j) = valij;
            Aloc(j,i) = valij;
        end
    end
    % store data
    yy = repmat(nodes,numel(nodes),1);
    xx = yy';
    adv = numel(xx);
    ii(pt+1:pt+adv) = xx(:);
    jj(pt+1:pt+adv) = yy(:);
    ss(pt+1:pt+adv) = Aloc(:);
    pt = pt + adv;
end
ii = ii(1:pt); jj = jj(1:pt); ss = ss(1:pt);
ASSEVEM = sparse(ii,jj,ss,NN,NN) + stab;
ASSEVEM(mesh.bdNodes,:) = 0;
ASSEVEM(:,mesh.bdNodes) = 0;
ASSEVEM(mesh.bdNodes,mesh.bdNodes) = speye(numel(mesh.bdNodes));

end