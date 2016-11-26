function [P,Q,bs,bc,y,trueG,predG,stdInd,crsInd] = SVDPluePlus(X,Y,P,Q,bs,bc,y,l2,l1,writefile)

maxIter = 40;
[N, M] = size(X);
err = []; 
	
for iter = 1:opt.maxIter
	fprintf(1, '\n\n\n\n iter = %d\n', iter);
	fprintf(writefile, '\n\n\n\n iter = %d\n', iter);
	sumtime=0;
	avetime=0;
	tic;
	for i=1:N
		for j = 1:M
			if X(i,j)>0
				I = find(X(i, :) ~= 0);
				lenI = length(I);
				sum0 = sum(y(I,:),1);
				meani = sum(X(I))/lenI;
				lenI = sqrt(length(I));
				impsum = (P(i,:)+sum0/lenI)*Q(j,:)';
				sum1 = X(i,j) - meani - bs(1,i) - bc(1,j) - impsum;
				P(i, :) = P(i, :)-lr*((-1)*Q(j,:)*sum1+l2*P(i ,:)+l1);
				Q(j, :) = Q(j, :)-lr*((-1)*(P(i,:)+sum0/lenI)*sum1+l2*Q(j ,:)+l1);
				P(i, :) = max(0,P(i, :));
				Q(j, :) = max(0,Q(j, :));
				bs(1,i) = bs(1,i)+lr*(sum1-l2*bs(1,i)-l1);
				bc(1,j) = bc(1,j)+lr*(sum1-l2*bc(1,j)-l1);
				y(I,:) = y(I,:)+lr*(repmat(sum1*Q(j,:)/lenI,length(I), 1)-l2*y(I,:)-l1);
				U=P;
				V=Q;
			end
		end
	end
end % iter


