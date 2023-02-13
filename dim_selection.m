function [dim] = dim_selection(Mask,Dec,Dec2,D,T)
    Base = zeros(1,D);
    Ind = SOLUTION([Base;Dec.*Mask]);
    Ind2 = SOLUTION([Base;Dec2.*Mask]);
    [FrontNo,~] = NDSort(Ind.objs,inf);
    [temp_FrontNo,~] = NDSort(Ind2.objs,inf);
    
    dim = find(FrontNo<=FrontNo(1,1))-1;
    temp_dim = find(temp_FrontNo<=temp_FrontNo(1,1))-1;
    NS_num = size(dim,2);
    

    if NS_num <= T  
        dim = union(temp_dim,dim);
    else
        dim = intersect(temp_dim,dim);
    end
end
