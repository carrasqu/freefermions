function σ(ϵ,Nboson)
  #1.0./(1.0+exp(β*(ϵ-μ)))
  lx=length(ϵ)
  nn=zeros(lx)
  nn[1:Nboson]=1.0
  return nn 
     
end



# parameters in the Hamiltonian
lx=3000 # lattice size

for lx=200:200:4000

t=-1.0 # hopping
V=0.0 # strenght of parabolic potential
β=2000.0 # inverse temperature
Mmax=2 # maximum matsubara frequency
wM=(0:Mmax)*2*π/β # matsubara frequencies
μ=0.0

Nboson=convert(Int64, lx/2)

# Generating Hamiltonian
H=zeros(lx,lx)

for i = 1:lx
    if i<lx
     H[i,i+1]=t
     H[i+1,i]=t
    end  
    H[i,i]=V*(i*1.0-(lx/2.0+0.5))^2
end
H[1,lx]=t
H[lx,1]=t


#diagonalizing H
(ϵ,U)=eig(H)

A=t* ( transpose(U[2:lx,1:lx])*U[1:lx-1,1:lx]-transpose(U[1:lx-1,1:lx])*U[2:lx,1:lx] +  transpose(reshape(U[1,:],(1,lx)) )*reshape(U[lx,:],(1,lx))                        -transpose(reshape(U[lx,:],(1,lx)) )*reshape(U[1,:],(1,lx))  ) 
B=-(t/lx)*diag(transpose(U[2:lx,1:lx])*U[1:lx-1,1:lx]+transpose(U[1:lx-1,1:lx])*U[2:lx,1:lx])
nf=σ(ϵ,Nboson)
kx=dot(nf,B)

Nboson=sum(nf)
#println(Nboson)
Λxx=complex(zeros(Mmax+1))

for i=2:Mmax
   for n=1:lx
       for np=1:lx
           if np!=n
              Λxx[i]=Λxx[i]+A[n,np]*A[np,n]*(nf[np]-nf[n])/(im*wM[i]-(ϵ[np]-ϵ[n]))     
           end
       end
   end
   println(lx," ", wM[i]," ",-(-kx+Λxx[i]/lx)," analytic ", 2.0/π )
end




end # end loop over system size



