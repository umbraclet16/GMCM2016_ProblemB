function [ non_duplicate ] = fun_delete_duplicate( origin )
non_duplicate=origin;
cont=0;
for i=2:numel(origin)
    for j=1:i-1
        if origin(i)==origin(j)
            non_duplicate(i-cont)=[];
            cont=cont+1;
            break
        end
    end
end