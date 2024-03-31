function [rank1,rank2] = rank_compute(Lambda,N)


mat1 = reshape(Lambda,N^2,N);
[~,mat1v,mat1r] = svd(mat1');
rank1 = rank(mat1);
mat2 = reshape(abs(mat1v(1,1)*mat1r(:,1)),N,N);
rank2 = rank(mat2);
end

