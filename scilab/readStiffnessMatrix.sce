// read mtx file

fd = mopen(uiputfile(["*.mtx"]),'r')

data = mgetl(fd,-1)

mat = zeros(6,6)

a = strsplit(data(8),",")
    mat(1,1) = strtod(a(1))

a = strsplit(data(9),",")
for i=1:size(a,1)
    mat(2,i) = strtod(a(i))
end

a = strsplit(data(10),",")
for i=1:size(a,1)
    mat(3,i) = strtod(a(i))
end

a = strsplit(data(11),",")
for i=1:size(a,1)
    mat(4,i) = strtod(a(i))
end

a = strsplit(data(12),",")
for i=1:size(a,1)
    mat(5,i) = strtod(a(i))
end
a = strsplit(data(13),",")
mat(5,5) = strtod(a(1))

a = strsplit(data(14),",")
for i=1:size(a,1)
    mat(6,i) = strtod(a(i))
end
a = strsplit(data(15),",")
mat(6,5) = strtod(a(1))
mat(6,6) = strtod(a(2))

mat = mat + mat' - diag(diag(mat))

mat(1:3,1:3) = mat(1:3,1:3) * 1e3
mat(4:6,4:6) = mat(4:6,4:6) * 1e-3

mat2 = mat
mat2(1,:) = mat(3,:)
mat2(3,:) = mat(1,:)
mat2(4,:) = mat(6,:)
mat2(6,:) = mat(4,:)

mat3 = mat2
mat3(:,1) = mat2(:,3)
mat3(:,3) = mat2(:,1)
mat3(:,4) = mat2(:,6)
mat3(:,6) = mat2(:,4)

mat = mat/2 + mat3/2

mclose(fd)
