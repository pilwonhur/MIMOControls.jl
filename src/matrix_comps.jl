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


function mqr(U;p=[])
	# Author: Pilwon Hur, Ph.D.
	#
	# modified qr decomposition
	# if m>=n, then it's the same as qr
	# if m<n, then generic qr does not care of the column order of Q matrix
	# mqr(U) will keep the column order of Q up to the level of rank.
	# In other words, the first r columns of Q are orthonomal vectors for the Range(U)
	# Within th Range(U), the order is not guaranteed to be the same as the column order of U.
	# However, it tends to keep the order if possible.
	# If you want to specify the order, then put the permutation information in p.
	#
	# Ex) mqr(U,p=[1,2])
	# In this case, the first 2 columns of U will be applied to the first 2 columns of Q with the same order.
	#
	# out=mqr(U)
	# out.Q, out.R, out.U, out.V
	# where out.U=Q[1:r]
	# out.V=out.U perp


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

function findprojector(U)
	# Author: Pilwon Hur, Ph.D.
	#
	# Input: a matrix U with independent basis column vectors
	# Output: returns a projector onto the range space of U

	# the following is a treatment for the case when U contains dependent vectors
	m,=size(U)

	F=mqr(U);
	# r=rank(U);
	# F=qr(U,Val(true));	# get the independent columns and it's permuation number
	# F=qr(U[:,sort(F.p[1:m])[1:r]])
	# V=F.Q[:,1:r]
	V=F.U;
	return V*inv(V'*V)*V';
end

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
		coord1=nullspace((I-Proj_contsubspace)*unobsv_subspace);

		# controllable/unobservable subspace
		if length(coord1)>0		# if t2 has elements
			ncontunobs,=reverse(size(coord1));
			t2=zeros(n,ncontunobs);
			[t2[:,i]=unobsv_subspace*coord1[:,i] for i=1:ncontunobs];

			F=mqr([t2 unobsv_subspace],p=(1:length(t2)));	# F.U will return orthonormal basis for unobservable subspace
			t4=F.U[:,ncontunobs+1:n-no]

			if ncontunobs==nc
				t1=[];
			else
				# F=qr([t2 cont_subspace],Val(true));
				F=mqr([t2 cont_subspace],p=(1:length(t2)));
				t1=F.U[:,ncontunobs+1:nc];
			end
		else 	# if t2 has no elements
			t4=unobsv_subspace;
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

function msub(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns G1-G2
	# G1: state space model of G1
	# G2: state space model of G2

	return madd(G1,-G2)	# note that -G2 automatically handles flipping sign of C and D.
end

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

function mscl(G::StateSpace,alpha)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns alpha*G
	# G: state space model of G
	# alpha: scalar

	return alpha*G
end

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

function abv(G1::StateSpace,G2::StateSpace)
	# Author: Pilwon Hur, Ph.D.
	#
	# returns [G1]
	#         [G2]
	# G1: state space model of G1
	# G2: state space model of G2

	return [G1;G2]
end

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

function eye(n)
	# Author: Pilwon Hur, Ph.D.
	#
	return Matrix{Float64}(I, n, n)
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

	solver=SCSSolver(eps=1e-6,max_iters=100000)
	m=Model(solver=solver)
	@variable(m,g2)
	@variable(m,X[1:n1,1:n1],SDP) 	# symmetric positive semidefinite
	@objective(m,Min,g2)
	@SDconstraint(m,[A'*X+X*A+C'*C X*B+C'*D;(X*B+C'*D)' D'*D-g2*eye(n2)]<=eps()*eye(n1+n2))
	JuMP.solve(m)

	return sqrt(getvalue(g2)), getvalue(X)
end

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
	
	solver=SCSSolver(eps=1e-6,max_iters=100000,verbose=0)
	m=Model(solver=solver)
	@variable(m,X[1:n,1:n],SDP) 	# symmetric positive semidefinite
	@objective(m,Min,tr(B'*X*P))
	@SDconstraint(m,A*X+X*A'+C'*C<=eps()*eye(n))
	JuMP.solve(m)

	return getobjectivevalue(m), getvalue(X)
end