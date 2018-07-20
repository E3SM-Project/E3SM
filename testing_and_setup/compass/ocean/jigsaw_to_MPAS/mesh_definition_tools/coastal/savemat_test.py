from scipy import io

io.savemat('test.mat',mdict={'test':[1,2,3,4,5],'num':5})
