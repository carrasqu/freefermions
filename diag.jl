function σ(ϵ,β,μ)
  1.0./(1.0+exp(β*(ϵ-μ)))
end

# parameters in the Hamiltonian
lx=500 # lattice size
t=-1.0 # hopping
V=0.0 # strenght of parabolic potential
β=2000.0 # inverse temperature
Mmax=2 # maximum matsubara frequency
wM=(0:Mmax)*2*π/β # matsubara frequencies
μ=0.0

for lx=200:200:4000

# Generating Hamiltonian
H=zeros(lx,lx)

β=10.0*lx
for i = 1:lx
    if i<lx
     H[i,i+1]=t
     H[i+1,i]=t
    end  
    H[i,i]=V*(i*1.0-(lx/2.0+0.5))^2
end

#diagonalizing H
(ϵ,U)=eig(H)

A=t*transpose(U[2:lx,1:lx])*U[1:lx-1,1:lx]-transpose(U[1:lx-1,1:lx])*U[2:lx,1:lx]
B=-(t/lx)*diag(transpose(U[2:lx,1:lx])*U[1:lx-1,1:lx]+transpose(U[1:lx-1,1:lx])*U[2:lx,1:lx])
nf=σ(ϵ,β,μ)
kx=dot(nf,B)

Nboson=sum(nf)
#println(Nboson)
Λxx=complex(zeros(Mmax+1))

for i=2:Mmax
   for n=1:lx
       for np=1:lx
           if np!=n
              Λxx[i]=Λxx[i]+A[n,np]*A[np,n]*(nf[np]-nf[n])/(sqrt(complex(-1.))*wM[i]-(ϵ[np]-ϵ[n]))     
               
           end
       end
   end
   println(lx," ", wM[i]," ",-(-kx+Λxx[i]/lx) )
end







end # lx loop
