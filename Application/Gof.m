% Goodness of fit
SnB=zeros( nloop-nwarmup, 1);
sigma=Omega(1, 2);
E1=U(:,1);
for i=1:nloop-nwarmup
    fprintf("loop: %d\n", i)
    Weights=weights_chain(:, i);
    g_now=get_g(dim, Weights, Bdensities, a);
    g2_now=g_now;
    g1_now = get_g1(g_now, dim, TR);
    f1_now = get_f1(g1_now, TR);
    cdf1_now = get_cdf1(f1_now, TR);
    Q1_now = get_Q1(cdf1_now, TR); 
    X_now=Q1_now(U);
    
    E2=zeros(N,1);
    for ii=1:N
        x1=X_now(ii, 1);
        pdf21=@(x2) (1-sigma^2)^(-1/2)*g2_now((x2-sigma*x1).^2/(1-sigma^2)+x1^2);
        c=integral(pdf21, -Inf, Inf);
        pdf21n=@(x2) pdf21(x2)/c;
        cdf21=@(x2) min(1, integral(pdf21n, -Inf, x2));
       
        E2(ii)=cdf21(X_now(ii,2));
        
    end
 
    E=[E1 E2];
    temp=1-E.^2;

    
    SnB_now=N/3^dim-sum(exp(sum(reallog(1-E.^2+10e-6),2)))/2^(dim-1);
    %Second term
    second_term=0;
    for ii=1:N
        for jj=1:N
            summ=0;
            for k=1:dim
                summ=summ+reallog(1-max(E(ii,k),E(jj,k))+10e-6);
            end
            second_term=second_term+exp(summ);
        end
    end
    SnB_now=SnB_now+second_term/N;
    SnB(i)=SnB_now;
end


figure
histogram(SnB, 'Normalization','pdf')
xlabel("$S_n^{(B)}$")
ylabel("Density")

