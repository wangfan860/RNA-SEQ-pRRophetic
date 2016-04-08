
x= read.csv("521OV__broad.mit.edu__ht_hg-u133a__gene.quantification__Apr-07-2016.csv")

gnames <- x[,1]

mat = data.matrix(x[,-1])
vec=c(1,2,3,4,5)
mat2=10^mat
mat3=sweep(mat2, 2, vec, `/`)
mat4=log10(mat3)
