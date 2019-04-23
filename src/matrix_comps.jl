mutable struct MQR
	Q::Any
	R::Any
	U::Any
	V::Any
end

mutable struct KALMANDECOMP
	T::Any
	t1::Any
	t2::Any
	t3::Any
	t4::Any
end

mutable struct SPLITSS 
	Gs::Any
	Gu::Any
	Gc::Any
end

mutable struct COMPLEX2REAL
	T::Any	# transform matrix: A=V*D*Vinv=V*(T*Tinv)*D*(T*Tinv)*Vinv=(V*T)*(Tinv*D*T)*(V*T)inv
	D::Any	# block diagonal matrix
	V::Any	# matrix of eigenvectors
end


"""`out=mqr(U;p=[])`

Author: Pilwon Hur, Ph.D.

Modified qr decomposition

If m>=n, then it's the same as qr()

If m<n, then generic qr does not care of the column order of Q matrix

mqr(U) will keep the column order of Q up to the level of rank.

In other words, the first r columns of Q are orthonomal vectors for the Range(U)

Within th Range(U), the order is not guaranteed to be the same as the column order of U.

However, it tends to keep the order if possible.

If you want to specify the order, then put the permutation information in p.

Ex) 
```julia 
out=mqr(U,p=[1,2])
out=mqr(U) # when you don't care about the order of the first 2 columns
out.Q
out.R
out.U
out.V
```

where `out.U=Q[1:r]`, `out.V=out.U perp`

In this case, the first 2 columns of U will be applied to the first 2 columns of Q with the same order.
"""
function mqr(U;p=[])
	m,n=size(U);
	r=rank(U);
	if m>=n
		F=qr(U)
		if r<m
			out=MQR(F.Q, F.R, F.Q[:,1:r], F.Q[:,r+1:m])
		else
			out=MQR(F.Q, F.R, F.Q[:,1:r], [])
		end
		return out;
	else 	# m<n
		F=qr(U,Val(true));	# get the independent columns and it's permuation number

		if length(p)==0	
			pnew=sort(F.p[1:m])
		else 	# when the permuation vector is provided
			pleft=setdiff(F.p,p);
			pleft=pleft[1:m-length(p)];
			pnew=vcat(p[:],pleft[:])
		end
		F=qr(U[:,pnew])
		if r<m
			out=MQR(F.Q, (F.Q)'*U, F.Q[:,1:r], F.Q[:,r+1:m])
		else
			out=MQR(F.Q, (F.Q)'*U, F.Q[:,1:r], [])
		end
		return out;
	end
end

"""`Proj=findprojector(U)`
Author: Pilwon Hur, Ph.D.

Input: a matrix U with independent basis column vectors
Output: returns a projector onto the range space of U

the following is a treatment for the case when U contains dependent vectors
"""
function findprojector(U)
	m,=size(U)

	F=mqr(U);
	# r=rank(U);
	# F=qr(U,Val(true));	# get the independent columns and it's permuation number
	# F=qr(U[:,sort(F.p[1:m])[1:r]])
	# V=F.Q[:,1:r]
	V=F.U;
	return V*inv(V'*V)*V';
end

"""`T=kalmandecomp(A,B,C,D)`
	`T=kalmandecomp(G::StateSpace)`
	`T=kalmandecomp(G::TransferFunction)`

Author: Pilwon Hur, Ph.D.

Returns the Kalman decomposition.
`T.T`: transform matrix
`T.t1`: basis for controllable and observable subspace
`T.t2`: basis for controllable and unobservable subspace
`T.t3`: basis for uncontrollable and observable subspace
`T.t4`: basis for uncontrollable and unobservable subspace
"""
function kalmandecomp(A,B,C,D)
	# Author: Pilwon Hur, Ph.D.
	#
	# n: number of states
	# m: number of outputs
	# r: number of inputs

	A=convert(Array{Float64,2},hcat(A));
	B=convert(Array{Float64,2},hcat(B));
	C=convert(Array{Float64,2},hcat(C));
	D=convert(Array{Float64,2},hcat(D));

	n,m=size(A);
	if n!=m
		error("Matrix A should be a square matrix.")
	end
	n1,r=size(B);
	if n!=n1
		error("Matrix B should have the same number of rows as the number of states.")
	end
	m,n1=size(C);
	if n!=n1
		error("Matrix C should have the same number of columns as the number of states.")
	end
	m1,r1=size(D);
	if m!=m1
		error("Matrix D should have the same number of rows as the number of outputs.")
	end
	if r!=r1
		error("Matrix D should have the same number of columns as the number of inputs")
	end

	Wc=ctrb(A,B);
	Wo=obsv(A,C);
	nc=rank(Wc);
	no=rank(Wo);

	# orthogonal controllable subspace
	# https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/
	# household based qr is not what I wanted. The order is totally different
	# F=qr(Wc,Val(true));
	F=mqr(Wc);
	cont_subspace=F.U;
	uncont_subspace=F.V;
	
	# orthogonal observable subspace
	# F=qr(Wo',Val(true));
	F=mqr(Wo');
	obsv_subspace=F.U;
	unobsv_subspace=F.V;

	# controllable and unobservable subspace
	Proj_contsubspace=findprojector(cont_subspace);
	t2=[];
	t1=cont_subspace;	
	t4=[];

	# find controllable/unobservable and uncontrollable/unobservable subspaces if unobservable subspace exists
	if no<n
		Proj_uncontsubspace=I-Proj_contsubspace;

		# if fully controllable
		if norm(Proj_uncontsubspace)<eps()*1000000;
			t2=unobsv_subspace;
			F=mqr([t2 cont_subspace],p=(1:size(t2)[2]));
			t1=F.U[:,size(t2)[2]+1:n];
		else # if not fully controllable
			coord1=nullspace(Proj_uncontsubspace*unobsv_subspace);

			# controllable/unobservable subspace
			if length(coord1)>0		# if t2 has elements
				ncontunobs,=reverse(size(coord1));
				t2=zeros(n,ncontunobs);
				[t2[:,i]=unobsv_subspace*coord1[:,i] for i=1:ncontunobs];

				F=mqr([t2 unobsv_subspace],p=(1:size(t2)[2]));	# F.U will return orthonormal basis for unobservable subspace
				t4=F.U[:,ncontunobs+1:n-no]

				if ncontunobs==nc
					t1=[];
				else
					# F=qr([t2 cont_subspace],Val(true));
					F=mqr([t2 cont_subspace],p=(1:size(t2)[2]));
					t1=F.U[:,ncontunobs+1:nc];
				end
			else 	# if t2 has no elements
				t4=unobsv_subspace;
			end
		end
	end

	ntemp=0;
	if length(t1)>0
		ntemp,=reverse(size(t1));
		temp=t1;
		if length(t2)>0
			nntemp,=reverse(size(t2));
			ntemp+=nntemp;
			temp=[temp t2];
		end
	else
		ntemp,=reverse(size(t2));
		temp=t2;
	end

	if length(t4)>0
		nntemp,=reverse(size(t4));
		ntemp+=nntemp;
		temp=[temp t4];
	end

	# temp is [t1 t2 t4]
	F=mqr(temp);
	t3=F.V;
	# if ntemp==n
	# 	t3=[];
	# else
	# 	F=qr(temp,Val(true));
	# 	t3=F.Q[:,ntemp+1:n];
	# end

	if length(t1)>0
		T=t1;
		if length(t2)>0
			T=[T t2];
		end
	else
		T=t2;
	end

	if length(t3)>0
		T=[T t3];
		if length(t4)>0
			T=[T t4];
		end
	else
		if length(t4)>0
			T=[T t4];
		end
	end

	out=KALMANDECOMP(T,t1,t2,t3,t4)

	return out
end

function kalmandecomp(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	return kalmandecomp(G.A,G.B,G.C,G.D)
end

function kalmandecomp(G::TransferFunction)
	# Author: Pilwon Hur, Ph.D.
	#
	return kalmandecomp(ss(G))
end

"""`Gs=minimumreal(G::StateSpace)`
	`Gs=minimumreal(G::TransferFunction)`
	`Gs=minimumreal(A,B,C,D)`

Author: Pilwon Hur, Ph.D.

Returns the Kalman decomposition.
`Gs`: minimum state space realization of the given system.
"""
function minimumreal(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#

	T=kalmandecomp(G);
	n,=reverse(size(T.t1));
	At=T.T\G.A*T.T;
	Bt=T.T\G.B;
	Ct=G.C*T.T;
	Dt=G.D;
	A1=At[1:n,1:n];
	B1=Bt[1:n,:];
	C1=Ct[:,1:n];
	D1=Dt;
	return ss(A1,B1,C1,D1)
end

function minimumreal(G::TransferFunction)
	# Author: Pilwon Hur, Ph.D.
	#
	return minimumreal(ss(G))
end

function minimumreal(A,B,C,D)
	# Author: Pilwon Hur, Ph.D.
	#
	return minimumreal(ss(A,B,C,D))
end


"""`Gs=mmult(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G1*G2`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function mmult(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G1*G2
	# G1: state space model of G1
	# G2: state space model of G2

	n1,=size(G1.A);
	n2,=size(G2.A);
	A=[G1.A G1.B*G2.C;zeros(n2,n1) G2.A];
	B=[G1.B*G2.D;G2.B];
	C=[G1.C G1.D*G2.C];
	D=G1.D*G2.D;
	return ss(A,B,C,D)
end

"""`Gs=madd(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G1+G2`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function madd(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G1+G2
	# G1: state space model of G1
	# G2: state space model of G2

	n1,=size(G1.A);
	n2,=size(G2.A);
	A=[G1.A zeros(n1,n2);zeros(n2,n1) G2.A];
	B=[G1.B;G2.B];
	C=[G1.C G2.C];
	D=G1.D+G2.D;
	return ss(A,B,C,D)
end

"""`Gs=msub(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G1-G2`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function msub(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G1-G2
	# G1: state space model of G1
	# G2: state space model of G2

	return madd(G1,-G2)	# note that -G2 automatically handles flipping sign of C and D.
end

"""`Gs=transp(G::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G'`
`G`: state space model of `G`
"""
function transp(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G'
	# G: state space model of G

	A=G.A';
	B=G.C';
	C=G.B';
	D=G.D';
	return ss(A,B,C,D)
end

"""`Gs=cjt(G::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G~`
`G`: state space model of `G`
"""
function cjt(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G~
	# G: state space model of G

	A=-G.A';
	B=-G.C';
	C=G.B';
	D=G.D';
	return ss(A,B,C,D)
end


"""`Gs=cj(G::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `G(-s)`
`G`: state space model of `G`
"""
function cj(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G~
	# G: state space model of G

	A=-G.A;
	B=-G.B;
	C=G.C;
	D=G.D;
	return ss(A,B,C,D)
end



"""`Gs=minv(G::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `inv(G)`
`G`: state space model of `G`
"""
function minv(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G~
	# G: state space model of G

	n,=size(G.D)
	if rank(G.D)<n
		error("D matrix is not invertiable! inv(G) does not exist!")
	end

	A=G.A-G.B*G.D\G.C;
	B=-G.B*inv(G.D);
	C=G.D\G.C;
	D=inv(G.D);
	return ss(A,B,C,D)
end

"""`Gs=mscl(G::StateSpace,alpha)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `alpha*G`
`G`: state space model of `G`
"""
function mscl(G::StateSpace,alpha)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns alpha*G
	# G: state space model of G
	# alpha: scalar

	return alpha*G
end

"""`Gs=sbs(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `[G1 G2]`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function sbs(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns [G1 G2]
	# G1: state space model of G1
	# G2: state space model of G2

	# n1,=size(G1.A);
	# n2,=size(G2.A);
	# r1,=reverse(size(G1.B));
	# r2,=reverse(size(G2.B));
	# A=[G1.A zeros(n1,n2);zeros(n2,n1) G2.A];
	# B=[G1.B zeros(n1,r2);zeros(n2,r1) G2.B];
	# C=[G1.C G2.C];
	# D=[G1.D G2.D];
	return [G1 G2]
end

"""`Gs=abv(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `[G1]`
								   `[G2]`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function abv(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns [G1]
	#         [G2]
	# G1: state space model of G1
	# G2: state space model of G2

	return [G1;G2]
end

"""`Gs=daug(G1::StateSpace,G2::StateSpace)`
Author: Pilwon Hur, Ph.D.

returns state space realization of `[G1  0]`
								   `[0  G2]`
`G1`: state space model of `G1`
`G2`: state space model of `G2`
"""
function daug(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns [G1  0]
	#         [0  G2]
	# G1: state space model of G1
	# G2: state space model of G2

	n1,=size(G1.A);
	n2,=size(G2.A);
	r1,=reverse(size(G1.B));
	r2,=reverse(size(G2.B));
	m1,=size(G1.C);
	m2,=size(G2.C);
	A=[G1.A zeros(n1,n2);zeros(n2,n1) G2.A];
	B=[G1.B zeros(n1,r2);zeros(n2,r1) G2.B];
	C=[G1.C zeros(m1,n2);zeros(m2,n1) G2.C];
	D=[G1.D zeros(m1,r2);zeros(m2,r1) G2.D];
	return ss(A,B,C,D)
end

"""`I=eye(n)`
Author: Pilwon Hur, Ph.D.

returns nxn identity matrix
"""
function eye(n)
	# Author: Pilwon Hur, Ph.D.
	#
	return Matrix{Float64}(I, n, n)
end

"""`I=eyec(n)`
Author: Pilwon Hur, Ph.D.

returns nxn complex identity matrix
"""
function eyec(n)
	# Author: Pilwon Hur, Ph.D.
	#
	return Matrix{Complex{Float64}}(I, n, n)
end

"""`Gs=hinflmi(G::StateSpace)`
	`Gs=hinflmi(G::TransferFunction)`
Author: Pilwon Hur, Ph.D.

returns hinf norm of the given system
`G`: state space model of `G`
"""
function LQRlmi(G,Q::Array{Float64,2},R::Float64)
    n2,=reverse(size(Q));
    solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=0)
    m = Model(solver=solver)
    
    if isa(G,Array{StateSpace{Float64,Array{Float64,2}},1})
        T = 2
        n1,=size(G[1].A);
        @variable(m,P[1:n1,1:n1],SDP)
        @objective(m,Min,-tr(P))
        
        for sys in G
            @SDconstraint(m,[(P*sys.A'+sys.A*P-sys.B*(R^-1)*sys.B') P;
                        P -inv(Q)]<=zeros((n1+n2),(n1+n2)))
        end
        
    elseif isa(G,StateSpace{Float64,Array{Float64,2}})
        T = 1
        n1,=size(G.A);
        @variable(m,P[1:n1,1:n1],SDP)
        @objective(m,Min,-tr(P))
        @SDconstraint(m,[(P*G.A'+G.A*P-G.B*(R^-1)*G.B') P;
                        P -inv(Q)]<=zeros((n1+n2),(n1+n2)))
    else
        print("Not in State Space or State Space list")
        return false
    end
    
    
    
    JuMP.solve(m)
    
    if T==2
        K = Array{Float64}[]
        for sys in G
            push!(K,(R^1)*sys.B'*inv(getvalue(P)))
        end
    else
        K = (R^1)*G.B'*inv(getvalue(P))
    end
    
    # returns P matrix and Gain
    return inv(getvalue(P)), K
end

function hinflmi(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	# 
	# Accepts the G, state space model
	# returns gamma and P using LMI

	A=G.A;
	B=G.B;
	C=G.C;
	D=G.D;
	n1,=size(A);
	n2,=reverse(size(B));

	# JuMP old version 0.18
	model=Model(solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=0))
	@variable(model,g2>=0)
	@variable(model,X[1:n1,1:n1],SDP) 	# symmetric positive semidefinite
	@objective(model,Min,g2)
	@SDconstraint(model,[A'*X+X*A+C'*C X*B+C'*D;(X*B+C'*D)' D'*D-g2.*eye(n2)]<=eps()*eye(n1+n2))
	JuMP.solve(model)
	return sqrt(getvalue(g2)), getvalue(X)

	# version 0.19 # For the moment, 0.19 is not stable. Downgrade to 0.18.5
    # m = Model(with_optimizer(SCS.Optimizer,eps=1e-6,max_iters=100000,verbose=0))
    # @variable(m,g2)
    # @variable(m,X[1:n1,1:n1],PSD)
    # @objective(m,Min,g2)
    # @constraint(m,g2>=eps())
    # @SDconstraint(m,[A'*X+X*A+C'*C X*B+C'*D;(X*B+C'*D)' D'*D-g2.*eye(n2)]<=zeros(n1+n2,n1+n2))
    # JuMP.optimize!(m)
    # return sqrt(JuMP.objective_value(m))


	# Convex version
	# solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=1)
	# g2=Variable(Positive())
	# X=Semidefinite(n1)
	# p=minimize(g2)
	# p.constraints+=(X in :SDP)
	# p.constraints+=(-1*[A'*X+X*A+C'*C X*B+C'*D;(X*B+C'*D)' D'*D-g2*eye(n2)] in :SDP)
	# solve!(p, solver)
	# return sqrt(p.optval), X.value
end

function hinflmi(G::TransferFunction)
	return hinflmi(ss(G))
end

"""`Gs=h2lmi(G::StateSpace)`
	`Gs=h2lmi(G::TransferFunction)`
Author: Pilwon Hur, Ph.D.

returns h2 norm of the given system
`G`: state space model of `G`
"""
function h2lmi(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	# 
	# Accepts the G, state space model
	# returns gamma and P using LMI
	# http://www.juliaopt.org/JuMP.jl/stable/refmodel.html
	
	A=G.A;
	B=G.B;
	C=G.C;
	D=G.D;
	n,=size(A);

	# JuMP version

	# old version	0.18
	solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=1)
	m=Model(solver=solver)
	@variable(m,X[1:n,1:n],SDP) 	# symmetric positive semidefinite
	@objective(m,Min,tr(B'*X*B))
	@SDconstraint(m,A*X+X*A'+C'*C<=eps()*eye(n))
	JuMP.solve(m)
	return sqrt(getobjectivevalue(m)), getvalue(X)

	# version 0.19 # For the moment, 0.19 is not stable. Downgrade to 0.18.5
    # m = Model(with_optimizer(SCS.Optimizer,eps=1e-6,max_iters=100000,verbose=0))
    # @variable(m,X[1:n,1:n],PSD)
    # @objective(m,Min,tr(B'*X*B))
    # @SDconstraint(m,A'*X+X*A+C'*C<=-eps()*eye(n) )
    # JuMP.optimize!(m)
    # return sqrt(JuMP.objective_value(m))
    # return sqrt(JuMP.objective_value(m)), JuMP.value(X)

    # Convex version
    # https://github.com/JuliaOpt/Convex.jl
    # https://convexjl.readthedocs.io/en/latest/solvers.html

 	# solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=1)
	# X=Semidefinite(n) # X=Variable(n,n)
	# p=minimize(tr(B'*X*B))
	# # p.constraints+=A'*X+X*A+C'*C<zeros(n,n)  # this is wrong. see below.
	# p.constraints+=(-(A'*X+X*A+C'*C)< in :SDP) # or p.constraints+=(A'*X+X*A+C'*C < 0) # Convex is using different syntax from JuMP
	# solve!(p, solver)
	# return sqrt(p.optval), X.value
end

function h2lmi(G::TransferFunction)
	return h2lmi(ss(G))
end

"""`Gs=h2gram(G::StateSpace)`
	`Gs=h2gram(G::TransferFunction)`
Author: Pilwon Hur, Ph.D.

returns h2 norm of the given system using Gramian
`G`: state space model of `G`
"""
function h2gram(G::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	# 
	# Accepts the G, state space model
	# returns gamma and P using LMI
	# http://www.juliaopt.org/JuMP.jl/stable/refmodel.html
	
	Wo=gram(G,:o);
	# return sqrt(getobjectivevalue(m)), getvalue(X)
	return sqrt(tr(G.B'*Wo*G.B))
end



"""`K2,Tzw=h2syn(G::StateSpace)`
	`K2,Tzw=h2gram(G::TransferFunction)`
Author: Pilwon Hur, Ph.D.

returns h2 controller of the given system G and its LFT
`G`: state space model of `G`
"""
function h2syn(G::StateSpace,ny,nu)
	# Author: Pilwon Hur, Ph.D.
	# 
	# Accepts the G, state space model
	# Accepts the G, state space model
	# returns hinf controller of the given system G and its LFT
	
	A=G.A;
	B=G.B;
	C=G.C;
	D=G.D;
	nstate=size(A,1);
	nin=size(B,2);
	nout=size(C,1);
	nw=nin-nu;
	nz=nout-ny;

	B1=reshape(B[:,1:nw],(nstate,nw));
	B2=reshape(B[:,nw+1:nin],(nstate,nu));
	C1=reshape(C[1:nz,:],(nz,nstate));
	C2=reshape(C[nz+1:nout,:],(ny,nstate));
	D11=reshape(D[1:nz,1:nw],(nz,nw));
	D12=reshape(D[1:nz,nw+1:nin],(nz,nu));
	D21=reshape(D[nz+1:nout,1:nw],(ny,nw));
	D22=reshape(D[nz+1:nout,nw+1:nin],(ny,nu));

	# controller design
	Ruu=D12'*D12;
	Ac=A-B2*inv(Ruu)*D12'*C1;
	Qc=C1'*(I-D12*inv(Ruu)*D12')*C1;
	# X2=care(Ac,B2,Qc,Ruu);
	X2=are(Ac,B2*inv(Ruu)*B2',Qc);
	F=-inv(Ruu)*(B2'*X2+D12'*C1);

	# observer design
	 Rww=D21*D21';
	Ao=A-B1*D21'*inv(Rww)*C2;
	Qo=B1*(I-D21'*inv(Rww)*D21)*B1'
	# Y2=care(Ao',C2',Qo,Rww)
	Y2=are(Ao',C2'*inv(Rww)*C2,Qo)
	L=-(Y2*C2'+B1*D21')*inv(Rww)

	# construct the H2 controller
	Ak=A+B2*F+L*C2;
	Bk=-L;
	Ck=F;
	Dk=zeros(size(Ck,1),size(Bk,2));
	K2=minimumreal(Ak,Bk,Ck,Dk);

	Tzw=lft(G,K2);

	return K2,Tzw
end


"""`Kinf,Tzw=hinfsyn(G::StateSpace,r,ny,nu)`
	`Kinf,Tzw=hinfsyn(G::TransferFunction,r,ny,nu)`
Author: Pilwon Hur, Ph.D.

returns hinf controller of the given system G and its LFT
`G`: state space model of `G`
"""
function hinfsyn(G::StateSpace,r,ny,nu)
	# Author: Pilwon Hur, Ph.D.
	# 
	# Accepts the G, state space model
	# returns hinf controller of the given system G and its LFT
	
	A=G.A;
	B=G.B;
	C=G.C;
	D=G.D;
	nstate=size(A,1);
	nin=size(B,2);
	nout=size(C,1);
	nw=nin-nu;
	nz=nout-ny;

	B1=reshape(B[:,1:nw],(nstate,nw));
	B2=reshape(B[:,nw+1:nin],(nstate,nu));
	C1=reshape(C[1:nz,:],(nz,nstate));
	C2=reshape(C[nz+1:nout,:],(ny,nstate));
	D11=reshape(D[1:nz,1:nw],(nz,nw));
	D12=reshape(D[1:nz,nw+1:nin],(nz,nu));
	D21=reshape(D[nz+1:nout,1:nw],(ny,nw));
	D22=reshape(D[nz+1:nout,nw+1:nin],(ny,nu));

	# controller design
	R1=D12'*D12;
	Ac=A-B2*inv(R1)*D12'*C1;
	Qc=C1'*(I-D12*inv(R1)*D12')'*(I-D12*inv(R1)*D12')*C1;
	R=B1*B1'/r^2-B2*inv(R1)*B2';
	
	# Xinf=care(Ac,I,Qc,inv(R));
	Xinf=are(Ac,R,Qc);
	F2=-inv(R1)*(B2'*Xinf+D12'*C1);
	
	# observer design
	R2=D21*D21';
	Ao=A-B1*D21'*inv(R2)*C2;
	Qo=B1*(I-D21'*inv(R2)*D21)*(I-D21'*inv(R2)*D21)'*B1';
	R=C1'*C1/r^2-C2'*inv(R2)*C2;
	
	# Yinf=care(Ao',I,Qo,inv(R));
	Yinf=care(Ao',R,Qo);
	Linf=-inv(I-1/r^2*Yinf*Xinf)*(Yinf*C2'+B1*D21')*inv(R2);

	F1=1/r^2 *(B1*B1'+Linf*D21*B1')*Xinf;

	

	# construct the Hinf controller
	Ak=A+B2*F2+Linf*C2+F1;
	Bk=-Linf;
	Ck=F2;
	Dk=zeros(size(Ck,1),size(Bk,2));
	Kinf=minimumreal(Ak,Bk,Ck,Dk);

	Tzw=lft(G,Kinf);

	return Kinf,Tzw
end


function lft(G1::StateSpace,G2::StateSpace)
	# if lower lft
	G2s=minimumreal(G2);
	nu,ny=size(G2s.D);
	G1s=minimumreal(G1);
	A=G1s.A;
	B=G1s.B;
	C=G1s.C;
	D=G1s.D;

	nstate=size(A,1);
	nin=size(B,2);
	nout=size(C,1);
	nw=nin-nu;
	nz=nout-ny;

	B1=reshape(B[:,1:nw],(nstate,nw));
	B2=reshape(B[:,nw+1:nin],(nstate,nu));
	C1=reshape(C[1:nz,:],(nz,nstate));
	C2=reshape(C[nz+1:nout,:],(ny,nstate));
	D11=reshape(D[1:nz,1:nw],(nz,nw));
	D12=reshape(D[1:nz,nw+1:nin],(nz,nu));
	D21=reshape(D[nz+1:nout,1:nw],(ny,nw));
	D22=reshape(D[nz+1:nout,nw+1:nin],(ny,nu));

	Ak=G2s.A;
	Bk=G2s.B;
	Ck=G2s.C;
	Dk=G2s.D;

	E=inv(I-D22*Dk);
	F=inv(I-Dk*D22);

	TzwA=[A+B2*F*Dk*C2 B2*F*Ck;
			Bk*E*C2 Ak+Bk*E*D22*Ck];
	TzwB=[B1+B2*F*Dk*D21;
			Bk*E*D21];
	TzwC=[C1+D12*F*Dk*C2 D12*F*Ck];
	TzwD=D11+D12*F*Dk*D21;

	Tzw=minimumreal(TzwA,TzwB,TzwC,TzwD);
	return Tzw
end


"""'X = ARE(A, R, Q) returns the stablizing solution (if it
    exists) to the continuous-time Riccati equation:

           A'*X + X*A - X*R*X + Q = 0

    assuming R is symmetric and nonnegative definite and Q is
    symmetric.
"""
function are(A,R,Q)
	# check the size of A matrix
	nstate,nc=size(A);
	if nstate!=nc
		error("Matrix A should be a square matrix.")
	end

	# check the size of Q matrix
	n,m=size(Q);
	if nstate!=n || nstate!=m
		error("Dimension of Matrix Q is not correct.")
	end

	# check the size of R matrix
	n,m=size(R);
	if nstate!=n || nstate!=m
		error("Dimension of Matrix R is not correct.")
	end

	# Hamiltonian matrix
	H=[A -R;-Q -A'];

	F=eigen(H);
	tol=1e-10;
	perm=Int64[];
	for i=1:2*nstate
        if abs(real(F.values[i]))<tol
           	error("No solution: Hamiltonian matrix has eigenvalues on the imaginary axis.");
           	break;
       	end
       	if real(F.values[i])<-tol
       		push!(perm,i);
       	end
	end

	X1=F.vectors[1:n,perm];
	X2=F.vectors[n+1:2*n,perm];

	# check if X1 is invertible
	L,U=lu(X1);
	if abs(prod(diag(U))) < tol
		error("No solution: Hamiltonian matrix is not complementary.");
	end

	X=X2/X1;
	return X;
end


function h2gram(G::TransferFunction)
	return h2gram(ss(G))
end

"""`Gs=hinfbis(G::StateSpace)`
	`Gs=hinfbis(G::TransferFunction)`
Author: Pilwon Hur, Ph.D.

returns hinf norm of the given system using Hamiltonian matrix and bisection
`G`: state space model of `G`
"""
function hinfbis(G::StateSpace,rl,ru)
	r=(ru+rl)/2
	rprev=0;
	A=G.A;
	B=G.B;
	C=G.C;
	D=G.D;

	while abs(r-rprev)>0.0001
    	rprev=r;
    	R=r^2*I-D'*D;
    	H=[A+B*inv(R)*D'*C B*inv(R)*B';
            -C'*(I+D*inv(R)*D')*C -(A+B*inv(R)*D'*C)'];
    	lambda=eigvals(H);
    	rutemp=ru;
    	for i=1:length(lambda)
        	if abs(real(lambda[i]))<0.0000001
            	rl=r;
            	ru=rutemp;
            	break;
        	end
        	ru=r;
    	end
    	r=(ru+rl)/2
	end
	return r
end

function hinfbis(G::TransferFunction,rl,ru)
	return hinfbis(ss(G),rl,ru)
end



function complex2real(A)
	# Convert complex eigendecomposition of real matrix
	# into eigen-like decomposition with real block matrices
	d,v=eigen(A);
	n=length(d);
	if !isreal(d)
		T=eyec(n);
		for i=1:n
			if !isreal(d[i]) && i<n
				if norm(d[i]-conj(d[i+1]))<eps()*1000
					T[i,i]=T[i+1,i]=0.5;
					T[i,i+1]=-0.5im;
					T[i+1,i+1]=0.5im;
				end
			end
		end

		# A=V*D*Vinv=V*(T*Tinv)*D*(T*Tinv)*Vinv=(V*T)*(Tinv*D*T)*(V*T)inv
		D=real(inv(T)*diagm(0=>d)*T);
		V=real(v*T);
	else
		T=eye(n);
	end

	out=COMPLEX2REAL(T,D,V);

	return out;
end

function splitSS(G::StateSpace)
	Gmin=minimumreal(G);
	d,v=eigen(Gmin.A);
	dreal=real(d);
	# p=sortperm(dreal);
	ps=(1:length(dreal))[dreal.<0];
	pu=(1:length(dreal))[dreal.>0];
	pc=(1:length(dreal))[dreal.==0];

	if isreal(d)	
		Anew=diagm(0=>d);
		Bnew=inv(v)*Gmin.B;
		Cnew=Gmin.C*v;
		Dnew=Gmin.D;		
	else 	# if complex
		T=complex2real(Gmin.A)

		Anew=T.D;
		Bnew=inv(T.V)*Gmin.B;
		Cnew=Gmin.C*T.V;
		Dnew=Gmin.D;		
	end 

	Gs=ss(Anew[ps,ps],Bnew[ps,:],Cnew[:,ps],zeros(size(Dnew)));
	Gu=ss(Anew[pu,pu],Bnew[pu,:],Cnew[:,pu],zeros(size(Dnew)));
	Gc=ss(Anew[pc,pc],Bnew[pc,:],Cnew[:,pc],Dnew);

	out=SPLITSS(Gs,Gu,Gc);
	return out;
end

function splitSS(G::TransferFunction)
	return splitSS(ss(G))
end

# function ss(G::TransferFunction,opt)
# 	if opt == "minimal"
# 		return minimumreal(G);
# 	else
# 		error("Your option should be \"minimal\".")
# 	end
# end

