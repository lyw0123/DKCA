function dim_saea = dim_base(Mask_dim,dim,num)  
    non_zero = zeros(1,size(Mask_dim,1));  
    for i = 1:size(Mask_dim,1)
        non_zero(i) =size(find(Mask_dim(i,:)),2);
    end
    index = find(non_zero==num);
    n = size(index,2); 
    dim_b = zeros(n,num); 
    for i = 1:n
        dim_b(i,:) = dim(find(Mask_dim(index(i),:)));
    end 
    dim_saea = unique(dim_b)'; 
end