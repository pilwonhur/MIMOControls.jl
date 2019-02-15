function kalmandecomp(A,B,C,D)
	# n: number of states
	# m: number of outputs
	# r: number of inputs

	n,m=size(A);
	if n!=m
		error("Matrix A should be a square matrix.")
	end
	n1,r=size(hcat(B));
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

	Wc=ctrb(A,hcat(B));
	Wo=obsv(A,C);
	nc=rank(Wc);
	no=rank(Wo);

	# orthogonal controllable subspace
	F=qr(Wc);
	cont_subspace=F.Q[:,1:nc];
	if nc<n
		uncont_subspace=F.Q[:,nc+1:n];
	else
		uncont_subspace=[];
	end

	# orthogonal observable subspace
	F=qr(Wo');
	obsv_subspace=F.Q[:,1:no];
	if no<n
		unobsv_subspace=F.Q[:,no+1:n];
	else
		unobsv_subspace=[];
	end

	# controllable and unobservable subspace
	Proj_contsubspace=cont_subspace*inv(cont_subspace'*cont_subspace)*cont_subspace';
	t2=[];
	t1=cont_subspace;	
	t4=[];
	if no<n 	# if unobservable subspace exists
		coord1=nullspace((I-Proj_contsubspace)*unobsv_subspace);
		if length(coord1)>0
			ncontunobs,=reverse(size(coord1));
			t2=zeros(n,ncontunobs);
			[t2[:,i]=unobsv_subspace*coord1[:,i] for i=1:ncontunobs];

			if ncontunobs<n-no
				F=qr([t2 unobsv_subspace]);
				t4=F.Q[:ncontunobs+1:n-no];
			end

			if ncontunobs==nc
				t1=[];
			else
				F=qr([t2 cont_subspace]);
				t1=F.Q[:,no+1:nc];
			end
		end
	end

	ntemp=0;
	if length(t1)>0
		ntemp,=size(t1);
		temp=t1;
		if length(t2)>0
			nntemp,=size(t2);
			ntemp+=nntemp;
			temp=[temp t2];
		end
	else
		ntemp,=size(t2);
		temp=t2;
	end

	if length(t4)>0
		nntemp,=size(t4);
		ntemp+=nntemp;
		temp=[temp t4];
	end

	if ntemp==n
		t3=[];
	else
		F=qr(temp);
		t3=F.Q[:,ntemp+1:n];
	end


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

	return T
end