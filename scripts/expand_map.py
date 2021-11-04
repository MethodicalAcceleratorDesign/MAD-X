from sympy import *
import numpy as np

def make_fortran_str(in_exp):
    in_str = str(in_exp)
    for i in range(0,10):
        print(in_str)
        print(str(i))
        in_str = in_str.replace(str(i),str(i)+'d0')
    return in_str


x,px,y,py,t,pt = symbols('x px y py t pt') # The canonical variables 
xf,pxf,yf,pyf,tf,ptf = symbols('xf pxf yf pyf tf ptf') # The canonical variables 
dl, beta = symbols('dl beta', positive = True, real = True)
l_pz, myzero, csq = symbols('l_pz, myzero csq', real = True)


dmapi = [x,py,y,py,t,pt]
msq = 1 + 2*pt*beta + pt**2 - px**2 - py**2
l_pz = dl / sqrt(msq)
xf = px*l_pz
pxf = px
yf = py*l_pz
pyf = py
tf = (dl*beta - (beta + pt) * l_pz)
ptf = pt
myzero = 0

dmapi = [x,px,y,py,t,pt]
dmapf = [xf,pxf,yf,pyf,tf,ptf]
te = np.empty((6,6,6),dtype=object)
re = eye(6)

re_out = open("re_out.txt", "w")
te_out = open("te_out.txt", "w")

for i in range(0,6):
    for j in range(0,6):
        re[i,j] = diff(dmapf[i],dmapi[j])
        if(re[i,j]!=0):
            m_i = i + 1
            m_j = j + 1
            re[i,j].simplify()
            tmp = re[i,j]
            tmp = tmp.subs(msq,csq)
            print(f're({m_i},{m_j}) = '+ str(tmp), file=re_out)
            for k in range(0,6):
                te[i,j,k] = diff(re[i,j],dmapi[k])
                m_k = k + 1
                if(te[i,j,k]!=myzero and te[i,j,k] is not None):
                    tmp = te[i,j,k]
                    #tmp.simplify()
                    tmp = tmp.subs(msq,csq)/2
                    print(f'te({m_i},{m_j},{m_k}) = '+ make_fortran_str(tmp) , file=te_out)

re_out.close()                    

J = Matrix([[0, 1, 0, 0, 0, 0 ],
           [-1, 0, 0, 0, 0, 0 ],
           [0 , 0, 0, 1, 0, 0 ],
           [0,  0,-1, 0, 0, 0 ],
           [0,  0, 0, 0, 0, 1 ],
           [0,  0, 0, 0,-1, 0 ]])
#t1 = te[0,:,:]
#f1 = t1*J*re + re.transpose()*J*t1
#f1.simplify()

res = re.transpose()*J*re
res.simplify()
exit()


with open('toMADX', 'w') as f:
    print('re(2,1) = '+ str(dpxdx), file=f)
    print('re(2,3) = '+ str(dpxdy), file=f)
    print('re(4,1) = '+ str(dpydx), file=f)
    print('re(4,3) = '+ str(dpydy), file=f)


fin = open("toMADX", "rt")
#output file to write the result to
fout = open("outToMADX.txt", "w")
#for each line in the input file
for line in fin:
    #read replace the string and write to output file
    a = line.replace('(L + Lint - Abs(L - Lint))','Leff')
    a = a.replace('(x**2 + y**2)','R')
    a = a.replace('(L + Lint)**2','Lsum')
    a = a.replace('(-L + Lint)**2','Ldiff')
    a = a.replace('4*x**2 + 4*y**2','4d0*R')

    fout.write(a)

#close input and output files
fin.close()
fout.close()