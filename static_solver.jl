using LinearAlgebra
#= Author: Alexander Greer

FUNCTION INPUTS

connecttable describes the ELEMENTS
   Column
   1: element number
   2: node A
   3: node B
   4: cross-sectional area
   5: elastic modulus
   6: coefficient of thermal expansion
   7: temperature change

nodetable describes the NODES
   Column
   1: node number
   2: x-coordinate
   3: y-coordinate
   4: z-coodinate
   5: x-displacement (if known, NaN if unknown)
   6: y-displacement (if known, NaN if unknown)
   7: z-displacement (if known, NaN if unknown)
   8: external force in x-direction (if known, NaN if unknown)
   9: external force in y-direction (if known, NaN if unknown)
   10: external force in z-direction (if known, NaN if unknown)

dof is the number of degrees of freedom in the truss


FUNCTION OUTPUTS
u = matrix of displacements for each node (size n x dof, where n is the number of nodes and dof is the number of degress of freedom)
   For row n:
   Column 1: x-displacement of node n
   Column 2: y-displacement of node n
   Column 3: z-displacement of node n

F = matrix of external force on each node (size n x dof)
   For row n:
   Column 1: x-component of external force on node n
   Column 2: y-component of external force on node n
   Column 3: z-component of external force on node n

sigma = vector of axial stress in each element (size m x 1, where m is the number of elements)
   For row m:
   Column 1: axial stress in element m
=#

#function [u, F, sigma] = static_solver(nodetable, connecttable, dof)
function static_solver(nodetable, connecttable, dof)
    #----------- PREPROCESSING -----------
    # coordinate matrix for each node with the specified dof
    co = nodetable[:,2:dof+1];

    # external force matrix for each node with the specified dof
    F = nodetable[:,8:dof+7];

    # displacement matrix for each node with the specified dof
    u = nodetable[:,5:dof+4];

    Nel = size(connecttable,1); # number of elements
    Nnodes = size(co,1);        # number of nodes

    #-----------  Initializing Variables -----------
    # initialize internal force vector
    sigma = fill(NaN,Nel,1); # how to do this in julia

    #initialize unitvector array for elements
    e = fill(NaN,dof,1,Nel);

    #initialize R array for elements
    R = fill(NaN,dof,dof,Nel);

    #initialize K hat vector for elements
    K_hat = fill(NaN,dof,dof,Nel);

    #initialize array of augmented K matrices for elements
    K_aug = zeros(Nnodes.*dof,Nnodes.*dof,Nel);

    #initialize the global K matrix
    K = zeros(Nnodes.*dof,Nnodes.*dof);

    #initializing the augmented temperature vector
    T_aug1 = zeros(Nnodes,dof,Nel);
    T_aug2 = zeros(Nnodes,dof,Nel);

    #initializing the global temperature effect vector
    T = zeros(Nnodes,dof);

    #-----------  Form Global Stiffness Matrix -----------
    for i = 1:Nel
        #make the e vector for each element
        #subtract the second from the first node
        #indexes are the node numbers for the element
        index1 = Int64(connecttable[i,2]);
        index2 = Int64(connecttable[i,3]);
        second_node = co[index2,:];
        first_node = co[index1,:];
        e[:,:,i] = (second_node - first_node)./(norm(second_node-first_node)); # check if norm is a function

        #make the R vector for each element
        R[:,:,i] = e[:,:,i]*e[:,:,i]';

        #make the K_hat vector for each element
        E = connecttable[i,5];
        A = connecttable[i,4];
        L = norm(second_node-first_node);
        K_hat[:,:,i] = (E.*A)./(L).*R[:,:,i];

        #----------- make the augmented stiffness matrix for each element-----------

        #place the positive and negative K hat vectors for an element in an
        #augmented matrix based off of node numbers
        K_aug[((index1-1)*dof+1):((index1-1)*dof+dof),((index1-1)*dof+1):((index1-1)*dof+dof),i] = K_hat[:,:,i];
        K_aug[((index1-1)*dof+1):((index1-1)*dof+dof),((index2-1)*dof+1):((index2-1)*dof+dof),i] = -1*K_hat[:,:,i];
        K_aug[((index2-1)*dof+1):((index2-1)*dof+dof),((index1-1)*dof+1):((index1-1)*dof+dof),i] = -1*K_hat[:,:,i];
        K_aug[((index2-1)*dof+1):((index2-1)*dof+dof),((index2-1)*dof+1):((index2-1)*dof+dof),i] = K_hat[:,:,i];

        #make the global stiffness matrix, add the augmented stiffness matrix
        #for the current element to the global stiffness matrix
        K = K + K_aug[:,:,i];

        #Calculate the temperature effects for each element
        alpha = connecttable[i,6];
        deltaT = connecttable[i,7];
        #create the augmented temperature matrices
        T_aug1[index1,:,i] = alpha*deltaT*E*A*e[:,:,i]';
        T_aug2[index2,:,i] = -1*alpha*deltaT*E*A*e[:,:,i]';
        #make the global Temmperature matrix, add the augmented temperature
        #matrices
        T = T + T_aug1[:,:,i] + T_aug2[:,:,i];
    end


    #-----------  Calculate all of the output variables -----------

    #manipulate the F and umatrix so that it is a column vector
    F_column = reshape(F',:,1); # we need to see what the reshape command is going to look like in Julia
    u_column = reshape(u',:,1);
    T_column = reshape(T',:,1);
    #Create a logical matrix for values that must be taken out
    F_logical = map(!isnan,F_column); # this does not work in julia
    #creating the reduces values for u, F, A, and K
    F_reduced = F_column[F_logical];
    u_reduced = u_column[F_logical];

    T_reduced = T_column[F_logical];
    mask = convert(Matrix{Bool}, F_logical*F_logical')
    #println(mask)
    K_reduced = K[mask]; # THIS IS INEFFICIENT - LOOK TO DO THIS BETTER
    K_reduced = reshape(K_reduced,length(F_reduced),length(F_reduced))
    #solve the values for u_reduced
    u_reduced = K_reduced\(F_reduced - T_reduced);

    #Put u_reduced back into u
    u_column[F_logical] = u_reduced;
    u = reshape(u_column',dof,Nnodes)';

    #solve the values for F
    F_vert = K*u_column + T_column;
    F = reshape(F_vert',dof,Nnodes)';

    #calculate the internal stress from the displacements
    for m = 1:Nel
        E = connecttable[m,5];
        alpha = connecttable[m,6];
        deltaT = connecttable[m,7];
        second_node = Int64(connecttable[m,3]);
        first_node = Int64(connecttable[m,2]);
        L = norm(co[second_node,:] - co[first_node,:]);
        sigma[m] = (E./L).*(((u[second_node,:]-u[first_node,:])'*e[:,:,m])[1] - alpha.*deltaT.*L); # also BS solutions with the [1]
    end
    # return the displacements, element forces and stresses
    return [u, F, sigma]

end
