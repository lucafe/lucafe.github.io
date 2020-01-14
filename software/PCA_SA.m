%%
% <latex>
% \begin{center}\vspace{-1em}
% 	\Large\textbf{Experiments on polynomial (chaos) approximation of 
%        maximum eigenvalue functions: Tutorial}\\
% 	\large\textit{Luca Fenzi\footnote{\texttt{luca.fenzi@kuleuven.be},
% \textit{Department of Computer Science, KU Leuven, Belgium}. } 
% \& Wim Michiels \normalsize \footnote{\texttt{wim.michiels@cs.kuleuven.be}, 
% \textit{Department of Computer Science, KU Leuven, Belgium}.}  }
% \end{center}
% </latex>

%%
% <latex> \vspace{-2em} </latex>

%%
% *Description.* _This tutorial describes the numerical experiments reported 
% in the article,_ Fenzi & Michiels (2018) ``Polynomial (chaos) approximation
% of maximum eigenvalue functions: efficiency and limitations",  _providing
% a template that can be modified for explorations of your own._

%%
% <latex> \vspace{-1.5em} </latex>

%%
% _The tutorial explores the polynomial approximation of smooth,
% non-differentiable and not even Lipschitz continuous benchmark functions 
% in the univariate and bivariate cases. The analyzed functions arise from 
% parameter eigenvalue problems; more in details, they are the real part of
% the rightmost eigenvalue (the so-called spectral abscissa)._

%%
% <latex> \vspace{-1.5em} </latex>

%%
% _The polynomial approximations are obtained by Galerkin and collocation 
% approaches. In the Galerkin approach, the numerical approximation of the 
% coefficients in the univariate case is achieved by extended 
% (or composite) Trapezoidal and Simpson's rules or by Gauss and Clenshaw-Curtis
% quadrature rules. For the bivariate case, the coefficients are approximated 
% by tensorial and non-tensorial Clenshaw-Curtis cubature rules, 
% based on tensor-product Chebyshev grid and Padua points, respectively. 
% The collocation approach interpolates the function on Chebyshev points 
% in the univariate case, while for the bivariate case the interpolant 
% nodes are given by tensor-product Chebyshev grid and Padua points._

%%
% <latex> \vspace{-2em} </latex>

%%
% <latex> 
% {\footnotesize\setcounter{tocdepth}{2}\tableofcontents} 
% </latex>


%%
% <latex> 
% \section{Introduction}\label{sec:intro} \vspace{-1em}
% </latex>

%%
% <latex>
% This tutorial documents the analysis performed in \cite{Fenzi2018}, 
% in a similar fashion w.r.t. \cite{Trefethen2012}. Most of the analysis can
% be carried out through the Chebfun software package. Additional MATLAB 
% functions, used for the analysis of the bivariate case 
% (Section \ref{sec:d2}), can be found at the end of this tutorial, in
% the Appendix \ref{appendix}.  
% </latex>


%%
% <latex>
% Everything here reported is fast to compute, in order that you can use it
% as a template to be modified for explorations of your own.  
% The results of the simulations used in \cite{Fenzi2018} are 
% obtained with similar codes. The present tutorial follows the same
% structure of \cite{Fenzi2018}, and we refer to the definitions and
% ambients of \cite{Fenzi2018} using the \textsc{small caps} format
% style. 
% </latex>

%%
% <latex>
% In what follows of this introduction, we explain how to 
% produce this text with the desired layout. Moreover, we define acronyms and a
% function handle to estimate the rates of convergence, 
% which are used in the upcoming sections.
% {\small\begin{description}
% \item[Chebfun download] Download Chebfun from the
% web site {\tt www.chebfun.org} and install it in
% your MATLAB path as instructed there.
% \item[PCA\_SA download] Request the MATLAB script  {\tt PCA\_SA}
% through the authors' emails. Publish this text with {\tt
% publish('PCA\_SA.m','latex')}. ({\tt PCA\_SA.tex} will appear
% in a subdirectory on your computer labeled {\tt html}.)
% \end{description}}
% </latex>

%%
% <latex>
% If you want to use the  same layout of \cite{Trefethen2012},
% then run the Chebfun script {\tt ATAPformats} before publishing
% this tutorial.
% For the layout of the figures, we set MATLAB to use \LaTeX\ to render all
% text of the images.  
% </latex>
%%
% <latex> \vspace{-2em} </latex>
set(0, 'DefaultTextInterpreter', 'LaTeX','DefaultAxesFontName', ...
    'LaTeX',  'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
    'DefaultLegendInterpreter', 'LaTeX');

%%
% <latex>
% In addition, we define the following acronyms for the MATLAB plotting options
% (cf.~\textit{e.g.}~{\tt help plot} in MATLAB) and for the more common
% $x-$ and $y-$ labels of the present text.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
LW='linewidth'; C='Color'; MS='MarkerSize';
P1='$P+1$'; M1='$M+1$'; C0='$|c_0-\tilde{c}_0^M|/|c_0|$';
AP='$\|\alpha-\alpha_P\|_\infty$'; P2='$P_d+1$';
TD='Total Degree'; MD='Maximal Degree';

%%
% The taxonomy of the spectral abscissa function is defined by the following
% terminology and highlighted by the usage of the corresponding colors:
%%
% <latex> \vspace{-2em} </latex>
cases={'SAE    ','MSSAEs ', 'MNSSAEs'};
col={[0,109,219]/255, [0,73,73]/255,[146,0,0]/255};

%%
% <latex> 
% We define a function handle to estimate the convergence rate 
% of the approximations $\{f_n\}_{n\in\mathbf{N}}$ w.r.t. $f$. If 
% $\{f_n\}_{n\in\mathbf{N}}$ has an algebraic index of convergence $r$, 
% \textit{i.e.}~$ f_n\sim\mathcal{O}(n^{-r})$ 
% (cf.~\textit{e.g.}~Section 2.3 in \cite{Boyd2001}),
% then $r$ can be estimated by 
% $$
% r=\frac{\log(\|f-f_{n_1}\|_\infty)-\log(\|f-f_{n_2}\|_\infty)}{\log(n_2)-\log(n_1)}, \quad n_1>n_2.
% $$
% Hence, the function handle is defined by the errors {\tt
% err}, which corresponds to $\|f-f_{n}\|_\infty$, and the vectors {\tt
% nn}, with more than $2$ ordered numbers $n\in\mathbf{N}$, where we want to estimate the
% convergence rates. 
% </latex>

%%
% <latex> \vspace{-2em} </latex>
rates=@(err,nn) (log(err(nn(2:end)))-log(err(nn(1:end-1))))./...
    (log(nn(1:end-1))-log(nn(2:end)));

%%
% <latex> 
% This function handle permits to empirically calculate the algebraic index of
% convergence $r$ for both numerical and approximation errors. 
% </latex>

%%
% <latex>
% {\it The present version of this tutorial is produced by Chebfun v.5.7.0 and MATLAB R2016b. 
% Ask the authors the script {\tt PCA\_SA.m}, which
% produces this text through the MATLAB command {\tt publish}. }
% </latex>

%%
% <latex> 
% \section{Parameter eigenvalue problem: Example 1}\label{sec:ex1} \vspace{-1em}
% </latex>

%%
% <latex> 
% Consider \textsc{Example 1} and the associated \textsc{Figure 1}, where the 
% spectra of the following matrices are analyzed for $x\in[-1,1]$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
A{1}=@(x) [exp(x),0;0,-1];
A{2}=@(x) [x,0;0,0];
A{3}=@(x) [0,x;1,0];

%%
% The corresponding spectral abscissa functions are
%%
% <latex> \vspace{-2em} </latex>
x=chebfun('x');
alpha{1}=exp(x);          % SAE
alpha{2}=x*(x>0);         % MSSAE
alpha{3}=sqrt(x)*(x>0);   % MNSSAE

%%
% To visually see that the spectral abscissa functions, previously defined,
% are the real part of rightmost eigenvalues of the matrices, we plot
% the real part of the spectra and the associated spectral abscissa.
%%
% <latex> \vspace{-2em} </latex>
xx=-1:0.01:1; lambda=zeros(2,length(xx));
for k=1:3
    for i=1:length(xx)
        lambda(:,i)=real(eig(A{k}(xx(i))));
    end, subplot(1,3,k)
    plot(alpha{k},C,col{k},LW, 2);   hold on,
    plot(xx, lambda, '.k',MS, 0.75); hold off, 
    title(cases{k}); ylim([-1.1, max(alpha{k})])
    xlabel('$\omega$'); ylabel('$\Re(\lambda)$'); 
end


%%
% <latex> 
% \section{Univariate polynomial approximation ($D=1$)}\label{sec:d1} \vspace{-1em}
% </latex>

%%
% <latex>
% We analyze the polynomial approximation up to order $P$ of the spectral abscissa
% functions {\tt alpha} with Galerkin and collocation approaches:
% $$ 
% \alpha\approx \alpha_P(\omega)=\sum_{i=0}^P \tilde{c}_ip_i(\omega),\quad P+1=100.$$
% </latex>
%%
% <latex> \vspace{-2em} </latex>
P=99;


%%
% <latex> 
% \subsection{Galerkin approach}\label{sec:gal1}  \vspace{-1em}
% </latex>

%%
% <latex>
% Legendre polynomials are set as polynomial basis $\{p_i\}_{i=0}^P$.
% {\tt L\{i+1\}} defines the $i$th Legendre polynomial, which
% is implemented in Chebfun by the command {\tt legpoly(i)}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
L=cell(1,P); for i=0:P, L{i+1}=legpoly(i); end

%%
% <latex>
% The $\rho$-norms of the orthogonal basis, with $\rho(\omega)=1/2$,
% are analytically known for the Legendre
% polynomial (cf.~\textit{e.g.}~formula 22.2.10 in \cite{Abramowitz1965}).
% We consider the square of these norms, {\it i.e.}~{\tt Pi(i)}$=\|
% p_i\|_\rho^2$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
Pi=@(i) 1./(2*i+1);

%%
% <latex>
% The coefficients, {\tt c(k,i+1)}$={c}_i$ of $\alpha_P(\omega)$, and
% the $\rho$-inner products, 
% {\tt a(k,i+1)}$=\langle {\tt alpha\{k\}}, {\tt L\{i+1\}} \rangle_\rho$,
% are analytically evaluated by the formula given in \textsc{Appendix A} of 
% \cite{Fenzi2018}, for {\tt alpha\{k\}} with {\tt k=1,2,3}. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
a=zeros(3,P+1); c=zeros(3,P+1);

%%
% <latex>
% First of all, we compute the coefficients of the $\rho$-inner product for
% the SAE:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
a(1,1)=(exp(1)-exp(-1))/2;c(1,1)=a(1,1)/Pi(0);
for j=2:P+1
   a(1,j)=(exp(1)+(-1)^j*exp(-1)-2*sum(c(1,j-1:-2:1)))/2;
   c(1,j)=a(1,j)/Pi(j-1);
end

%%
% <latex>
% Hence, we consider the MSSAEs and MNSSAEs cases:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for i=1:P+1
 if mod(i-1,2)==0,  j=(i-1)/2;
  a(2,i)=((-1)^j)*gamma(j-0.5)/(4*gamma(-0.5)*gamma(j+2));
  a(3,i)=((-1)^j)*gamma(j-1/4)*gamma(3/4)/...
         (4*gamma(-1/4)*gamma(j+7/4));
 else,              j=(i-2)/2; 
  a(2,i)=1/6*(i==2);
  a(3,i)=((-1)^j)*gamma(j+1/4)*gamma(5/4)/...
         (4*gamma(1/4)*gamma(j+9/4));
 end
end
c(2,:)=a(2,:)./Pi(0:P); c(3,:)=a(3,:)./Pi(0:P); close;

%%
% <latex>
% MATLAB is not a symbolic software and the operations to compute
% {\tt c} and {\tt a} can be effected by computational errors. The SAE case, 
% {\tt alpha\{1\}}, is 
% particularly effected by these errors, since the coefficients are
% evaluated by recursion (\textit{i.e.}~the computational error at iteration $i$ is
% amplified in the following iterations $j>i$). To avoid the amplification
% of computational errors due to the recursion, we observe
% $$
% \| \alpha\|_\rho^2=\sum_{i=0}^\infty c_i^2\|p_i\|_\rho^2.
% $$
% Hence, an upper bound on the coefficient $c_j$ can be derived:
% $$
% \| \alpha\|_\rho^2-\sum_{i=0}^{j-1} c_i^2\|p_i\|_\rho^2=
% \sum_{i=j}^\infty c_i^2\|p_i\|_\rho^2\geq c_j^2\|p_j\|_\rho^2,
% $$
% and applied to the SAE case, where $\| {\tt
% alpha\{1\}}\|_\rho^2=(e^{2}-e^{-2})/4$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
semilogy(1:20, abs(c(1,1:20)),'.k'); A12=(exp(2)-exp(-2))/4;
for j=2:P+1
    if Pi(j-1)*c(1,j)^2>A12-sum((c(1,1:j-1).^2).*Pi(0:j-2))
        c(1,j:end)=0;a(1,j:end)=0;break
    end
end, hold on
semilogy(1:j-1,abs(c(1,1:j-1)),'or',MS,8); hold off,xlabel(P1);
ylabel('$|c_i|$'); title('Galerkin coefficients of SAE');

%%
% <latex> 
% \subsubsection{Approximation error}\label{sec:gal1ae} \vspace{-1em}
% </latex>

%%
% <latex> 
% The error of truncating the polynomial series up to order $P$
% (\textit{i.e.}~the approximation error)
% is considered, assuming that the coefficients  (computed in 
% Section \ref{sec:gal1}) are not affected by any error.
% First the polynomial approximation {\tt alphaP\_G} is constructed, and then
% the error in $\infty$-norm, {\tt error}, is considered.
% Through Chebfun, the $L^\infty$ error for the SAE and MSSAEs cases is 
% evaluated up to machine precision. For the MNSSAEs, 
% the error can be correctly computed only in the interval $[-1,0)$; for 
% the interval $[0,1]$ the $\infty$-norm of the error is approximated by 
% the maximum error on $10^3$ equidistant point in $[0,1]$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
alphaP_G=cell(3,P+1);error=zeros(3,P+1); xxINF=linspace(0,1,1e3);
for k=1:3
 for i=1:P+1
  if i==1, alphaP_G{k,1}=c(k,1)*L{i};
  else, alphaP_G{k,i}=c(k,i)*L{i}+alphaP_G{k,i-1}; end
  if k<3, error(k,i)=norm(alpha{k}-alphaP_G{k,i},inf); else
   error(k,i)=max([norm(alphaP_G{k,i}*(x<=0),inf),...
     abs(alpha{3}(xxINF)-alphaP_G{k,i}(xxINF))]);
  end
 end,loglog(1:P+1,error(k,:),C,col{k},LW, 2); hold on
end, hold off, grid on, ylim([1e-8,10]);
xlabel(P1); ylabel(AP); title('Error Galerkin approach'); 


%%
% <latex> 
% The convergence error plot is analogous to \textsc{Figure 2}. The
% relative errors $\|\alpha-\alpha_P\|_\infty/\|\alpha\|_\infty$ can be
% obtained by the previous analysis dividing the error of the SAE case by
% $\|\alpha\|_\infty=e$, \textit{i.e.}~{\tt error(1,:)/exp(1)}.
% The convergence rates $\mathcal{O}(P^{-r})$ for MSSAEs and MNSSAEs cases
% are estimated by the median of the rates obtained by the previously 
% defined function handle {\tt rates}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=2:3
    o=rates(error(k,:),2:2:P);
    fprintf([cases{k},' %4.8f\n'], median(o));
end




%%
% <latex> 
% The polynomial approximation $\alpha_P$  with
% $P=${\tt Pmax} obtained by Galerkin approach, can be compared w.r.t.
% the original spectral abscissa functions $\alpha$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
Pmax=8;
for k=1:3, subplot(1,3,k), 
    plot(alpha{k},C,col{k},LW, 1); hold on
    plot(alphaP_G{k,Pmax},':k',LW, 2); hold off, 
    xlabel('$\omega$'); title(cases{k});
    ylabel('$\alpha(\omega)$ and $\alpha_P(\omega)$'); 
end

%%
% <latex> 
% \subsubsection{Numerical error}\label{sec:gal1ne} \vspace{-1em}
% </latex>

%%
% <latex> 
% In this section, classical integration methods and interpolatory quadrature
% rules approximate the coefficients 
% $$
% c_i=\frac{1}{2\|p_i\|_\rho^2}\int_{-1}^1\alpha(\omega)
% p_i(\omega)\mathrm{d} \omega.
% $$
% The computation of $\|p_i\|_\rho^2$ can be omitted in this analysis, 
% since $\|p_i\|_\rho^2$ cancels out in the
% evaluation of the relative errors. W.l.g. only the first coefficient
% $c_0$ is analyzed, \textit{i.e}
% </latex>
%%
% <latex> \vspace{-2em} </latex>
i=1; 



%%
% <latex>
% The classical integration methods based on equally spaced points, 
% here considered, are:
% \begin{description}
% \item[Extended Trapezoidal rule] formula 25.4.2 in \cite{Abramowitz1965}
% $$
% \int_{x_0}^{x_M}f(x)\mathrm{d} x=h\left(\frac{f_0}{2}+\sum_{i=1}^{M-1}
% f_i+\frac{f_M}{2}\right) -\frac{Mh^3}{12}f^{\prime\prime}(\xi), \ 
% \xi\in(x_0,x_M).
% $$
% \item[Extended Simpson's rule] formula 25.4.6 in \cite{Abramowitz1965}{\small 
% $$
% \int_{x_0}^{x_{2n}}f(x)\mathrm{d} x=\frac{h}{3}\left(f_0+\sum_{i=1}^{n-1}
% 2f_{2i}+4f_{2i-1}+f_{2n}\right) -\frac{nh^5}{90}f^{(4)}(\xi), \
% \xi\in(x_0,x_{2n}).
% $$}
% \end{description}
% Compute the numerical error introduced by extended Trapezoidal rule, 
% {\tt et}, using at most {\tt M}$+1$ equally spaced points.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
M=500; et=NaN(3,M);
for k=1:3
  for m=1:M, h=2/m; xx=-1:h:1; % m+1 points
    at=trapz(xx,alpha{k}(xx).*L{i}(xx))/2; et(k,m)=abs(at-a(k,i));
  end
end, close;

%%
% <latex>
% The slowest convergence rate of the error is, hence, achieved by
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 loglog((1:2:M)+1,et(k,1:2:M)/abs(a(k,i)),':',C,col{k},LW, 2);
 hold on
end, hold off, grid on, ylim([1e-6,1]);xlabel(M1); ylabel(C0)
title('Error for $\tilde{c}_0^M$ with extended Trapezoidal rule')
 set(gca, 'XLim', [2 M], 'XTick', [2,10, 100 M]);


%%
% <latex>
% Hence we consider the numerical error due to the approximation of extended
% Simpson's rule, {\tt es}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
es=NaN(3,floor((M+1)/2));  % Simpson's error
for k=1:3
 for m=2:2:M,  h=2/m; xx=-1:h:1; f=alpha{k}*L{i}; 
  as=h/6*(f(xx(1))+2*sum(f(xx(3:2:m-1)))+4*sum(f(xx(2:2:m)))+...
      f(xx(m+1)));  es(k,m/2)=abs(as-a(k,i));
 end
end

%%
% <latex>
% The slowest convergence rates are, hence, achieved by
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 loglog((2:4:M)+1,es(k,1:2:M/2)/abs(a(k,i)),'--',C,col{k},LW, 2);
 hold on
end, hold off,grid on, ylim([1e-12,1]); xlabel(M1); ylabel(C0)
title('Error for $\tilde{c}_0^M$ with extended Simpson''s rule')
 set(gca, 'XLim', [3 M], 'XTick', [3,10, 100 M]);

%%
% <latex>
% The last two convergence error plots furnish \textsc{Figure 3}.
% The corresponding orders of convergence $r$ for these two classical
% integration methods, $\mathcal{O}(M^{-r})$,  are given by:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
disp('       Trapezoidal rule       Simpson''s rule');
for k=1:3
 ot=rates(et(k,:),1:2:M); os=rates(es(k,:),(2:4:M)/2);
 fprintf([cases{k},'   %10.8f %20.8f\n'], median(ot), median(os));
end

%%
% <latex>
% The interpolatory quadrature rules, here considered, are the
% Clenshaw-Curtis quadrature rule and the Gauss quadrature rule, based on the
% Chebyshev and Legendre points, respectively. 
% The weights and the nodes of these quadrature rules can be computed by
% the Chebfun functions {\tt chebpts} and {\tt legpts}. 
% In addition, Clenshaw-Curtis quadrature
% rule can be handled by the command {\tt sum} in Chebfun, which employ
% Fast Fourier Transform; this latter method is simple and faster when few
% integrands are involved \cite{Trefethen2012}, as in our case. Compute,
% hence, the numerical errors of Clenshaw-Curtis {\tt ec} and Gauss Legendre
% {\tt el} quadrature rules. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
ec=NaN(3,M);   el=NaN(3,M);  
for k=1:3
   for m=1:M
     ac=sum(chebfun(alpha{k}*L{i},m+1))/2;
     [n,w]=legpts(m+1);  al=w*(alpha{k}(n).*L{i}(n))/2;
     ec(k,m)=abs(ac-a(k,i));  el(k,m)=abs(al-a(k,i));
   end     
end

%%
% <latex>
% \textsc{Figure 4} considers only the slowest convergence rates, as the
% following code.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
m_slow=[1,2:2:M];
for k=1:3
 loglog(m_slow+1,el(k,m_slow)/abs(a(k,i)),'--',C,col{k},LW, 2);
 hold on
 loglog(m_slow+1,ec(k,m_slow)/abs(a(k,i)),'-',C,col{k},LW, 2);
end, hold off,grid on, ylim([1e-16,1]); xlabel(M1); ylabel(C0);
title('Error for $\tilde{c}_0^M$ with interpolatory quadrature rules')
set(gca, 'XLim', [2 M], 'XTick', [2,10, 100 M]);

%%
% <latex>
% The convergence rates are analogous, and the SAE case
% convergences faster than $\mathcal{O}(M^{-r})$
% for given natural numbers $r$. However, the interpolatory quadrature rules do
% not improve the orders of convergence for the non-smooth behaviors of the
% spectral abscissa functions. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
disp('       Clenshaw-Curtis      Gauss-Legendre');
for k=2:3
 oc=rates(ec(k,:),m_slow);ol=rates(el(k,:),m_slow);
 fprintf([cases{k},'   %1.8f %20.8f\n'], median(oc), median(ol));
end

%%
% <latex>
% We compute the Galerkin polynomial approximation {\tt alphaPM\_G}, whose
% coefficients are approximated by Clenshaw-Curtis quadrature rules with 
% {\tt M}$+1=110$ points. Other than indicating the approximation errors of the
% polynomial approximation {\tt err\_PMa}, we consider also the numerical
% error {\tt err\_PMc}, \textit{i.e.}~$|c_i-\tilde{c}_i^M|$, due to the 
% approximation of the coefficient $\tilde{c}_i^M$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
err_PMc=NaN(3,P+1); err_PMa=NaN(3,P+1); M=109;
for k=1:3, alphaPM_G=[];
 for p=1:P+1
  ac=sum(chebfun(alpha{k}*L{p},M+1))/(2*Pi(p-1));
  err_PMc(k,p)=abs(ac-c(k,p));
   if p==1, alphaPM_G=ac*L{p};
   else, alphaPM_G=ac*L{p}+alphaPM_G;
   end
   if k<3, err_PMa(k,p)=norm(alpha{k}-alphaPM_G,inf); else
     err_PMa(3,p)=max([norm(alphaPM_G*(x<=0),inf),...
            abs(alpha{3}(xxINF)-alphaPM_G(xxINF))]);
   end
 end
end



%%
% <latex>
% The numerical and approximation errors for Galerkin approach with 
% Clenshaw-Curtis quadrature rule are presented.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
   loglog(1:P+1,err_PMc(k,:),'-.',C,col{k},LW, 1);  hold on
   loglog(1:P+1,err_PMa(k,:),'-',C,col{k},LW, 2);
end, hold off, grid on, ylim([1e-16,1]);
xlabel(P1);ylabel('Numerical \& approximation errors');
title('Galerkin approach with Clenshaw-Curtis quadrature rule')

%%
% <latex>
% Since the numerical errors do not dominate the approximation errors, the
% order of convergence are analogous to the ones previously obtained. In
% particular, for the non-smooth cases we have
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=2:3
    o=rates(err_PMa(k,:),1:P-20);
    fprintf([cases{k},' %4.8f\n'], median(o));
end

%%
% <latex>
% Modifying this section, you can test the advice given in \textsc{Section
% 4.1.3} of \cite{Fenzi2018} for $M<99$. 
% Furthermore, you can test the convergence rates 
% of the Galerkin polynomial approximations whose coefficients are given,
% for example, by extended Trapezoidal or Simpson's rules or by Gauss Legendre
% quadrature rule. 
% </latex>

%%
% <latex> 
% \subsection{Collocation approach}\label{sec:col1} \vspace{-1em}
% </latex>


%%
% <latex>
% The near-best polynomial approximation (\textit{i.e.}~the 
% interpolant on Chebyshev points) and the best polynomial
% approximation in $L^\infty$ sense can be easily evaluated by Chebfun 
% with the functions {\tt chebfun} and {\tt minimax}, respectively. 
% The first input of these MATLAB commands is the function that we want to
% approximate, while the second input is the number of interpolant point
% for {\tt chebfun} and the degree of the best polynomial approximation
% for {\tt minimax}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
alphaP_b=cell(3,P+1);  err_b=zeros(3,P+1);   % Best
alphaP_nb=cell(3,P+1); err_nb=zeros(3,P+1);  % Near-Best
for k=1:3
 for i=1:P+1
  alphaP_b{k,i}=minimax(alpha{k},i-1); 
  alphaP_nb{k,i}=chebfun(alpha{k},i); 
    if k<3, err_b(k,i)=norm(alpha{k}-alphaP_b{k,i},inf);
      err_nb(k,i)=norm(alpha{k}-alphaP_nb{k,i},inf); 
    else
      err_b(k,i)=max([norm(alphaP_b{k,i}*(x<=0),inf),...
         abs(alpha{3}(xxINF)-alphaP_b{k,i}(xxINF))]);
      err_nb(k,i)=max([norm(alphaP_nb{k,i}*(x<=0),inf),...
         abs(alpha{3}(xxINF)-alphaP_nb{k,i}(xxINF))]);
    end
 end
 loglog(1:P+1,err_b(k,:),':',C,col{k},LW, 1.5); hold on
 loglog(1:P+1,err_nb(k,:),C,col{k},LW, 2);
end, hold off, xlabel(P1);ylabel(AP);ylim([1e-16,1]),grid on
title('Error Best and Near-Best polynomial approximation'); 

%%
% <latex>
% The previous image corresponds to \textsc{Figure 5}.
% The convergence rates are similar, and comparable with the results 
% obtained by the Galerkin approach. Indeed, for the non-smooth cases, the
% empirical convergence rates are estimated by:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
disp('          Near-Best        Best');
for k=2:3
 ob=rates(err_b(k,:),2:2:P); onb=rates(err_nb(k,:),2:2:P);
 fprintf([cases{k},'  %10.8f %14.8f \n'], median(ob), median(onb));
end

%%
% <latex> 
% \subsection{Extra: analysis on the error function and on the mean}
% \label{sec:extra1}\vspace{-1em}
% </latex>

%%
% <latex>
% In this section, we exploit side aspects of the analysis conducted in the
% numerical experiments of \cite{Fenzi2018}. 
% </latex>

%%
% <latex>
% First of all, we compare the error curves of Galerkin approach w.r.t.
% best and near-best approximations, such that the degree of the polynomial 
% approximation is $P=${\tt Pmax}. The MATLAB warnings are disabled, since 
% the error curves, computed by Chebfun, are not accurate 
% up to machine precision. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3, warning off; subplot(1,3,k);
    plot(alpha{k}-alphaP_b{k,Pmax},'-k',LW,1); hold on
    plot(alpha{k}-alphaP_nb{k,Pmax},'-r',LW,1.5); 
    plot(alpha{k}-alphaP_G{k,Pmax},'-b',LW,2);
    xlabel('$\omega$'); ylabel('error curves'); title(cases{k});
end, warning on

%%
% <latex>
% Hence, we consider the corresponding errors in the $2$-norm and in the
% $\infty$-norm. Only the SAE and MSSAEs cases are considered since they
% can easily computed by Chebfun. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:2, b=alpha{k}-alphaP_b{k,Pmax}; 
    nb=alpha{k}-alphaP_nb{k,Pmax}; g=alpha{k}-alphaP_G{k,Pmax};
 fprintf([' *',cases{k},'  Inf-Norm         2-Norm\n']);
 fprintf('Best      %10.8e   %14.8e\n', norm(b,inf), norm(b,2));
 fprintf('Near-Best %10.8e   %14.8e\n', norm(nb,inf), norm(nb,2));
 fprintf('Galerkin  %10.8e   %14.8e\n', norm(g,inf), norm(g,2));
end 

%%
% <latex>
% Then, we compare the mean of the PC expansion associated to the best and 
% near-best polynomial approximation. Indeed, the convergence rates of the
% numerical errors in Section \ref{sec:gal1ne} can be interpreted as the convergence
% rates of the mean, evaluated by the corresponding integration method. 
% To this  end, we transform  the coefficients of the best 
% and near-best polynomial approximations from the Chebyshev to the Legendre
% bases. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
close; eb=zeros(3,P+1);  enb=zeros(3,P+1);  
for k=1:3, for p=1:P+1
  cl=cheb2leg(chebcoeffs(alphaP_b{k,p}));
  eb(k,p)=abs(cl(1)-c(k,1))/abs(c(k,1));
  cl=cheb2leg(chebcoeffs(alphaP_nb{k,p}));
  enb(k,p)=abs(cl(1)-c(k,1))/abs(c(k,1));    
           end,    
    loglog(1:P+1,eb(k,:),':',C,col{k},LW, 1.5); hold on
    loglog(1:P+1,enb(k,:),'-.',C,col{k},LW, 2);
end, hold off, grid on, xlabel(P1); ylabel(C0);ylim([1e-16,1]);
title('Error for $\tilde{c}_0$ with Best \& Near-Best approximations')
 
%%
% <latex> 
% \section{Bivariate polynomial approximation ($D=2$)}\label{sec:d2} \vspace{-1em}
% </latex>

%%
% <latex> 
% In this section we consider the polynomial approximation of the spectral
% abscissa of parameter varying eigenvalue problems with $D=2$.
% </latex>


%%
% <latex> 
% The polynomial approximation is truncated up to degree 
% $P_d=${\tt Pdmax}.
% The number of coefficients $P$ grows exponentially w.r.t.
% the multivariate (total and maximal) degree $P_d$ of
% the polynomial approximation. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
Pdmax=19; pp=1:Pdmax+1;
plot(pp,(pp+1).*(pp+2)/2,'b',pp,(pp+1).^2,'r'); 
xlabel(P2); ylabel(P1); legend(TD,MD,'Location','northwest');

%%
% <latex> 
% Set $P_d\in\mathbf{N}$, the maximal degree polynomial basis 
% $\{p_i^{[m]}\}_{i=0}^{P}$  contains the total degree basis 
% $\{p_i^{[t]}\}_{i=0}^{P}$, 
% \textit{i.e.}~$\{p_i^{[t]}\}_{i=0}^{P}\subseteq\{p_i^{[m]}\}_{i=0}^{P}$.
% These bases are constructed by the inverse of an associated pairing 
% function. The inverse of the Rosenberg-Strong pairing function, associated 
% to the maximal degree polynomial basis, is first considered.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
i1=zeros(1,(Pdmax+1)^2); i2=i1;
for i=0:length(i1)-1, t=floor(sqrt(i));
    if t>i-t^2, i1(i+1)=i-t^2; i2(i+1)=t;
    else,       i1(i+1)=t;     i2(i+1)=t^2+2*t-i;
    end
end

%%
% <latex> 
% The inverse of the Cantor pairing function is associated to the 
% total degree polynomial basis. However, we mainly consider the maximal degree 
% polynomial basis, since $\{p_i^{[t]}\}_{i=0}^{P}\subseteq\{p_i^{[m]}\}_{i=0}^{P}$,
% and we map the ordering associated to the Cantor pairing to the one obtained by
% Rosenberg-strong pairing function, {\it i.e.}
% $$
% {\tt c2r}(i)=j, \quad where\  \pi_1(i)=\pi_\infty(j).
% $$
% </latex>
%%
% <latex> \vspace{-2em} </latex>
i1c=zeros(1,nchoosek(Pdmax+2,2)); i2c=i1c; c2r=NaN(size(i1)); 
for i=0:length(i1c)-1
    w=floor((sqrt(8*i+1)-1)/2); t=(w^2+w)/2;
    i1c(i+1)=i-t;  i2c(i+1)=w-i1c(i+1);
    for j=0:length(i1)-1
        if i2c(i+1)==i2(j+1) && i1c(i+1)==i1(j+1)
            c2r(i+1)=j+1; break
        end
    end
end

%%
% <latex> 
% To test the previous pairing functions, we construct a figure similar to \textsc{Figure 6},
% which
% associates the numbers $\mathbf{N}\times\mathbf{N}$ to the corresponding
% value of the pairing functions. The construction follows the so-called {\it
% shell} of the pairing function \cite{Szudzik2017}, which is associated to the polynomial
% basis with degree equal to $P_d$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for pd=0:5
    subplot(1,2,1);  axis square; grid on   % Total
    for j=(pd+1)*pd/2+1:(pd+1)*(pd+2)/2
        text(i1(c2r(j)),i2(c2r(j)),num2str(j-1))
    end,  title('$\pi_1(i_1,i_2)=i$');    
    subplot(1,2,2); axis square; grid on    % Maximal
    for j=((pd)^2+1:(pd+1)^2)
        text(i1(j),i2(j),num2str(j-1))
    end,  title('$\pi_\infty(i_1,i_2)=i$')
end
s1=subplot(1,2,1); xlabel('$i_1$'); ylabel('$i_2$');
s2=subplot(1,2,2); xlabel('$i_1$'); ylabel('$i_2$');
set([s1,s2],'XLim',[0 5+.5],'XTick', 0:5,'YLim',[0 5+.5]);


%%
% <latex> 
% \textbf{Example 2.}
% The delay parameter varying eigenvalue problems of \textsc{Example 2} can be
% linearized by the Infinitesimal generator approach given
% $(\omega_1,\omega_2)\in\mathbf{S}$. The MATLAB function {\tt DelayedOscillator}  
% (given in the Appendix \ref{appendix}) discretizes the infinite dimensional linear eigenvalue
% problem associated to the oscillator with feedback delay system into a
% finite standard eigenvalue problem. 
% </latex>

%%
% <latex> 
% For our porpoise, we consider
% the discretization {\tt IG}$=20$ obtained by the method proposed by 
% \cite{Breda2015} (\textit{i.e.}~with {\tt stru=0}), a similar result 
% can be achieved
% by the approach of \cite{Jarlebring2010}, setting {\tt stru=1}.
% Moreover, we set the number of points ${\tt NT}^2$ for the
% approximation of the $\infty$-norm, and the linear transformation ${\tt
% Lt}:[-1,1]^2\to\mathbf{S}$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
IG=20; stru=0; NT=100;
S1=[0.9,1.1]; Lt1=@(o1) (S1(2)-S1(1))*((o1+1)/2)+S1(1);
S2=[0.1,0.2]; Lt2=@(o2) (S2(2)-S2(1))*((o2+1)/2)+S2(1);


%%
% <latex> 
% The approximation of the $\infty$-norm is achieved by ${\tt NT}^2$ 
% evaluations of the spectral abscissa in the domain $\mathbf{S}$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
XT=linspace(-1,1,NT); YT=XT; O1=Lt1(XT); O2=Lt2(YT);
[OO1,OO2]=meshgrid(O1,O2); AlphaT=zeros(NT,NT,3); 
for k=1:3
 for i=1:NT
  for j=1:NT
   [A,B]=DelayedOscillator(OO1(i,j),OO2(i,j),k,IG,stru);
   AlphaT(i,j,k)=max(real(eig(A,B)));    
  end
 end
end, close; 

%%
% <latex> 
% The behaviors of the spectral abscissa functions vary w.r.t. 
% the controllers of \textsc{Table 1}. The SAE case is a smooth bivariate function
% on $\mathbf{S}$, MSSAEs and MNSSAEs cases present non-differentiable and
% non-Lipschitz curves in the domain $\mathbf{S}$, respectively. The
% non-smooth behavior of MSSAEs is due to the crossing of two active eigenvalues,
% which do not overlap, while in the MNSSAEs a triple eigenvalue splits, 
% as shown in Figure 3 of \cite{Fenzi2017}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3,  subplot(1,3,k); meshc(OO1,OO2,AlphaT(:,:,k)); 
  xlabel('$\omega_1$'); ylabel('$\omega_2$');
  title(cases{k});zlabel('$\alpha(\omega_1,\omega_2)$'); 
end

%%
% <latex> 
% \subsection{Galerkin approach}\label{sec:gal2}\vspace{-1em}
% </latex>

%%
% <latex> 
% In this section, we compute the reference values of the
% coefficients and the corresponding polynomial approximations. Hence we
% consider the numerical errors and 
% we analyze the advice given in \textsc{Section 4.1.3} of 
% \cite{Fenzi2018}.
% </latex>

%%
% <latex> 
% \subsubsection{Approximation error}\label{sec:gal2ae} \vspace{-1em}
% </latex>


%%
% <latex> 
% The reference values are computed on {\tt mpT} Padua
% points and {\tt mcT} tensor product Chebyshev grid. The 
% corresponding total number of points are {\tt MpT}$+1$ and {\tt McT}$+1$,
% respectively. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
mpT=99; MpT=nchoosek(mpT+2,2)-1; mcT=99; McT=((mcT+1)^2)-1;

%%
% <latex> 
% The tensorial Clenshaw-Curtis cubature rule is constructed by
% tensor product of the Clenshaw-Curtis quadrature rule on {\tt mcT+1}
% Chebyshev nodes. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
tensor=@(D1) [repmat(D1,length(D1),1),repelem(D1,length(D1))];
[X,W]=chebpts(mcT+1); WcT=prod(tensor(W'),2); WcT=WcT/sum(WcT);
XYT=tensor(X); XcT=XYT(:,1); YcT=XYT(:,2);

%%
% <latex> 
% The non-tensorial Clenshaw-Curtis cubature rule is based on Padua
% points, implemented in Chebfun by the function {\tt
% paduapts}. The weights associated to the Padua points are evaluated by 
% {\tt pdwtsMM} function (cf.~\cite{Caliari2011}), included in 
% Appendix \ref{appendix}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
[XYT]=paduapts(mpT); XpT=XYT(:,1); YpT=XYT(:,2);
WpT=pdwtsMM(mpT); WpT=WpT/sum(WpT); close


%%
% <latex> 
% We observe that {\tt length(WpT)=MpT+1} and {\tt length(WcT)=McT+1}.
% </latex>

%%
% <latex>
% \textbf{Remark.} \textit{Modifying the tensor product formula given by
% the function handle {\tt tensor}, it is possible to construct integration
% rules for arbitrary dimension $D$.}
% </latex>

%%
% <latex>
% We approximate the coefficients ${c}_i^{M^\star}$ 
% associated to the bases $\{p_j^{[t]}\}_{j=0}^{P}$ and $\{p_j^{[m]}\}_{j=0}^{P}$
% by the pairing functions $\pi_1$ and $\pi_\infty$. The bases are constructed 
% by univariate Legendre polynomials {\tt L\{i+1\}}, defined in 
% Section \ref{sec:gal1}. We consider only  $p_i\in\{p_i^{[m]}\}_{i=0}^P$, 
% since  $\{p_i^{[t]}\}_{i=0}^P\subseteq\{p_i^{[m]}\}_{i=0}^P$ can be
% evaluated through the function {\tt c2r}.
% </latex>


%%
% <latex> 
% The inner products $\langle \alpha, p_i\rangle_\rho$  are first approximated 
% by non-tensorial Clenshaw-Curtis cubature
% rule, obtaining {\tt apT(i+1)}. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
apT=zeros(3,(Pdmax+1)^2);
for k=1:3, sa=zeros(MpT+1,1);
 for i=1: MpT+1
  [A,B]=DelayedOscillator(Lt1(XpT(i)),Lt2(YpT(i)),k,IG,stru);
  sa(i)=max(real(eig(A,B)));
 end
 for i=1:(Pdmax+1)^2
   apT(k,i)=sum(WpT.*sa.*L{i1(i)+1}(XpT).*L{i2(i)+1}(YpT));
 end
end


%%
% <latex> 
% Hence, we approximate the coefficients by tensorial Clenshaw-Curtis cubature
% rule, obtaining {\tt acT(i+1)}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
acT=zeros(3,(Pdmax+1)^2);   
for k=1:3, sa=zeros(McT+1,1);
 for i=1: McT+1
  [A,B]=DelayedOscillator(Lt1(XcT(i)),Lt2(YcT(i)),k,IG,stru);
  sa(i)=max(real(eig(A,B)));
 end
 for i=1:(Pdmax+1)^2
  acT(k,i)=sum(WcT.*sa.*L{i1(i)+1}(XcT).*L{i2(i)+1}(YcT));
 end
end

%%
% <latex> 
% We normalize the polynomial basis 
% $\{p_i^{[m]}\}_{i=0}^P$ by $\|p_i\|_\rho^2$.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
L2=zeros(NT,NT,(Pdmax+1)^2); [XXT,YYT]=meshgrid(XT,YT);
for i=1:(Pdmax+1)^2
 L2(:,:,i)=L{i1(i)+1}(XXT).*L{i2(i)+1}(YYT)/(Pi(i1(i))*Pi(i2(i)));
end
%%
% <latex> 
% \textbf{Remark.} \textit{ {\tt L2} is a tensor of dimension {\tt
% Nt$\times$Nt$\times$(Pdmax+1)$^2$}, it is possible to store it only if
% {\tt Nt} and {\tt Pdmax+1} are small.}
% </latex>

%%
% <latex> 
% At this point, we evaluate the polynomial approximation, obtained by 
% reference values coefficients, and the corresponding error in
% $\infty$-norm. We start with maximal degree polynomial approximation, 
% whose coefficients are evaluated by tensorial Clenshaw-Curtis cubature
% rules. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
err_cT=zeros(3,Pdmax+1);
for k=1:3
    SA=acT(k,1)*ones(size(AlphaT(:,:,k)));
    err_cT(k,1)=max(max(abs(SA-AlphaT(:,:,k))));
    for pd=1:Pdmax
        for j=((pd)^2+1:(pd+1)^2)
            SA=SA+acT(k,j)*L2(:,:,j);
        end
        err_cT(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
    end
end

%%
% <latex> 
% We compute the error of the polynomial approximation with total
% degree $P_d$ whose coefficients are approximated by non-tensorial
% Clenshaw-Curtis cubature rule. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
err_pT=zeros(3,Pdmax+1); 
for k=1:3
    SA=apT(k,1)*ones(size(AlphaT(:,:,k)));
    err_pT(k,1)=max(max(abs(SA-AlphaT(:,:,k)))); 
    for pd=1:Pdmax
        for j=(pd+1)*pd/2+1:(pd+1)*(pd+2)/2
           SA=SA+apT(k,c2r(j))*L2(:,:,c2r(j));
        end
        err_pT(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
    end
end, close

%%
% <latex> 
% The convergences of the errors can be shown by the following code.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 subplot(1,2,1),loglog(pp,err_pT(k,:),C,col{k},LW,2);hold on
 subplot(1,2,2),loglog(pp,err_cT(k,:),C,col{k},LW,2);hold on
end, AX=[1,20,1e-16,1];
subplot(1,2,1),grid on,xlabel(P2);ylabel(AP);title(TD);axis(AX)
subplot(1,2,2),grid on,xlabel(P2);ylabel(AP);title(MD);axis(AX)


%%
% <latex>
% The non-smooth cases converge as $\mathcal{O}(P_d^{-r})$, where $r$ is
% estimated by:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
disp('        Total Degree     Maximal Degree');
for k=2:3
 ot=rates(err_pT(k,:),7:2:Pdmax); om=rates(err_cT(k,:),1:2:Pdmax);
 fprintf([cases{k},'  %10.8f %14.8f \n'], median(ot), median(om));
end


%%
% <latex> 
% \subsubsection{Numerical error}\label{sec:gal2ne} \vspace{-1em}
% </latex>

%%
% <latex>
% In this section, we analyze the numerical error introduced by
% non-tensorial and tensorial Clenshaw-Curtis cubature rules in the
% computation of the coefficient ${c}_j$, in the ordering
% associated to the pairing function $\pi_1$ and $\pi_\infty$. For
% simplicity we consider {\tt j=0}, since the first coefficient is
% independent from the multivariate degree, {\it i.e.}~${c}_0={c}_0^{[t]}={c}_0^{[m]}$. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
j=0; mm=1:50; Mp=(mm+1).*(mm+2)/2-1; Mc=(mm+1).^2-1;

%%
% <latex>
% Analogously to the
% previous Section \ref{sec:gal1ne}, we consider the relative errors on the inner
% product. We start considering the numerical errors of the non-tensorial
% Clenshaw-Curtis cubature rule. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
ep=zeros(3,length(mm));
for m=mm 
    [XY]=paduapts(m); X=XY(:,1); Y=XY(:,2);
    W=pdwtsMM(m); W=W/sum(W); sa=zeros(Mp(m)+1,1);
    for k=1:3
      for i=1:Mp(m)+1
        [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
        sa(i)=max(real(eig(A,B)));
      end
      as=sum(W.*sa.*L{j+1}(X).*L{j+1}(Y));
      ep(k,m)=abs(as-apT(k,j+1));
    end
end

%%
% <latex>
% Hence, we consider the numerical errors of the tensorial
% Clenshaw-Curtis cubature rule. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
ec=zeros(3,length(mm));
for m=mm
    [X,W]=chebpts(m+1); XY=tensor(X); X=XY(:,1); Y=XY(:,2); 
    W=prod(tensor(W'),2); W=W/sum(W); sa=zeros(Mc(m)+1,1);
    for k=1:3
     for i=1:Mc(m)+1
      [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
      sa(i)=max(real(eig(A,B)));
     end
     as=sum(W.*sa.*L{j+1}(X).*L{j+1}(Y));
     ec(k,m)=abs(as-acT(k,j+1));
    end
end; close

%%
% <latex>
% The numerical errors for non-tensorial and tensorial Clenshaw-Curtis
% cubature rules are shown in the following figure, which
% corresponds to \textsc{Figure 8} in \cite{Fenzi2018}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 loglog(Mp+1,ep(k,:)/abs(apT(k,j+1)),':',C,col{k},LW,2); hold on
 loglog(Mc+1,ec(k,:)/abs(acT(k,j+1)),'--',C,col{k},LW,2); 
end, hold off, axis([Mp(1)+1,Mc(end)+1,1e-15,10]); xlabel(M1);
grid on,  ylabel('$|c_j^{M^\star}-c_j^{M}|/|c_j^{M^\star}|$')
title(['Numerical error for $c_j^M$ with $j=$',num2str(j)]);

%%
% <latex>
% The convergence rates in this case fluctuates a lot and it is not easy to
% empirally compute the convergence rates of these cubature rules. The
% estimation can be achieved by looking at the slope of the loglog plot. 
% </latex>

%%
% <latex> 
% \subsubsection{Decoupling numerical and approximation errors}\label{sec:gal2dec} \vspace{-1em}
% </latex>

%%
% <latex>
% In this section, we first evaluate the polynomial approximations which follow the 
% natural  choices given in \textsc{Section 4.1.3} in \cite{Fenzi2018}. 
% Then, as a counter check, we construct polynomial approximations
% which do not respect the advice, even though the number of points of the
% cubature rules is bigger than the number of approximated coefficients.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
mp=[25, 30]; MP=(mp+1).*(mp+2)/2-1;
mc=[25, 15]; MC=(mc+1).^2-1;


%%
% <latex>
% We can observe that {\tt mp(1)} and {\tt mc(1)} represent natural
% choices for polynomial approximations with total and maximal degree up to
% {\tt Pdmax}, respectively. {\tt mc(2)} and {\tt mp(2)} do not satisfy
% the advice for polynomial approximations with total and maximal degree up to
% {\tt Pdmax}, respectively. Indeed {\tt mc(2)<Pdmax} and 
% {\tt mp(2)<2*Pdmax}, even though {\tt MC(2)>nchoosek(Pdmax+2,2)-1} and  
% {\tt MP(2)>(Pdmax+1)$^2$}. 
% </latex>


%%
% <latex>
% We fist consider the total degree polynomial approximations, constructed 
% following the natural choice.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
[XY]=paduapts(mp(1)); W=pdwtsMM(mp(1)); W=W/sum(W);  
net=zeros(3,Pdmax+1); aet=net; X=XY(:,1); Y=XY(:,2);
for k=1:3, sa=zeros(MP(1)+1,1);
  for i=1:MP(1)+1
    [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
    sa(i)=max(real(eig(A,B)));
  end
   at=sum(W.*sa.*L{1}(X).*L{1}(Y)); net(k,1)=abs(at-apT(k,1));
   SA=at*ones(size(AlphaT(:,:,k)));
   aet(k,1)=max(max(abs(SA-AlphaT(:,:,k)))); 
 for pd=1:Pdmax
    for j=(pd+1)*pd/2+1:(pd+1)*(pd+2)/2
      at=sum(W.*sa.*L{i1(c2r(j))+1}(X).*L{i2(c2r(j))+1}(Y)); 
      net(k,pd+1)=net(k,pd+1)+abs(at-apT(k,c2r(j)));
      SA=SA+at*L2(:,:,c2r(j));
    end
    aet(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
 end
end


%%
% <latex>
% Hence, we evaluate the maximal degree polynomial approximation following
% the advice.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
[X,W]=chebpts(mc(1)+1); XY=tensor(X); W=prod(tensor(W'),2); 
nem=zeros(3,Pdmax+1); aem=nem; W=W/sum(W); X=XY(:,1); Y=XY(:,2);
for k=1:3, sa=zeros(MC(1)+1,1);
 for i=1:MC(1)+1
  [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
  sa(i)=max(real(eig(A,B)));
 end
  am=sum(W.*sa.*L{1}(X).*L{1}(Y)); nem(k,1)=abs(am-apT(k,1));
  SA=am*ones(size(AlphaT(:,:,k)));
  aem(k,1)=max(max(abs(SA-AlphaT(:,:,k)))); 
  for pd=1:Pdmax
   for j=((pd)^2+1:(pd+1)^2)
    am=sum(W.*sa.*L{i1(j)+1}(X).*L{i2(j)+1}(Y)); 
    nem(k,pd+1)=nem(k,pd+1)+abs(am-apT(k,j)); SA=SA+am*L2(:,:,j);
   end
   aem(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
  end
end

%%
% <latex>
% The  convergence of the approximation errors, obtained following the 
% natural choices, is not affected by the numerical errors, as illustrated
% in the following panes, which correspond to \textsc{Figure 7}.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 subplot(1,2,1), loglog(pp,aet(k,:),C,col{k},LW,2); 
 hold on, loglog(pp,net(k,:),'-.',C,col{k},LW,1.5); 
 subplot(1,2,2), loglog(pp,aem(k,:),C,col{k},LW,2); 
 hold on, loglog(pp,nem(k,:),'-.',C,col{k},LW,1.5);
end
subplot(1,2,1),grid on, xlabel(P2); title(TD); axis(AX)
subplot(1,2,2),grid on, xlabel(P2); title(MD); axis(AX)

%%
% <latex>
% In a similar fashion we compute the total degree polynomial
% approximation, whose coefficients are computed via {\tt mc(2)} tensorial
% Clenshaw-Curtis cubature rule. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
[X,W]=chebpts(mc(2)+1); XY=tensor(X); W=prod(tensor(W'),2); 
net=zeros(3,Pdmax+1); aet=net; W=W/sum(W); X=XY(:,1); Y=XY(:,2); 
for k=1:3, sa=zeros(MC(2)+1,1);
 for i=1:MC(2)+1
  [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
  sa(i)=max(real(eig(A,B)));
 end
 at=sum(W.*sa.*L{1}(X).*L{1}(Y)); net(k,1)=abs(at-apT(k,1));
 SA=at*ones(size(AlphaT(:,:,k)));
 aet(k,1)=max(max(abs(SA-AlphaT(:,:,k)))); 
 for pd=1:Pdmax
  for j=(pd+1)*pd/2+1:(pd+1)*(pd+2)/2
  at=sum(W.*sa.*L{i1(c2r(j))+1}(X).*L{i2(c2r(j))+1}(Y)); 
  net(k,pd+1)=net(k,pd+1)+abs(at-apT(k,c2r(j)));
  SA=SA+at*L2(:,:,c2r(j));
  end
  aet(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
 end
end

%%
% <latex>
% Hence, we consider the maximal degree polynomial
% approximation, whose coefficients are computed via {\tt mp(2)} non-tensorial
% Clenshaw-Curtis cubature rule. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
[XY]=paduapts(mp(2)); W=pdwtsMM(mp(2)); W=W/sum(W); 
nem=zeros(3,Pdmax+1); aem=nem; X=XY(:,1); Y=XY(:,2);  
for k=1:3, sa=zeros(MP(2)+1,1);
 for i=1:MP(2)+1
  [A,B]=DelayedOscillator(Lt1(X(i)),Lt2(Y(i)),k,IG,stru);
  sa(i)=max(real(eig(A,B)));
 end
 am=sum(W.*sa.*L{1}(X).*L{1}(Y)); nem(k,1)=abs(am-apT(k,1));
 SA=am*ones(size(AlphaT(:,:,k)));
 aem(k,1)=max(max(abs(SA-AlphaT(:,:,k)))); 
 for pd=1:Pdmax
  for j=((pd)^2+1:(pd+1)^2)
  am=sum(W.*sa.*L{i1(j)+1}(X).*L{i2(j)+1}(Y)); 
  nem(k,pd+1)=nem(k,pd+1)+abs(am-apT(k,j));
  SA=SA+am*L2(:,:,j);
 end
 aem(k,pd+1)=max(max(abs(SA-AlphaT(:,:,k))));
 end
end; close

%%
% <latex>
% Analogously to \textsc{Figure 9}, the approximation error is dominated by
% the numerical errors of the cubature rules, as soon as 
% the advice of \textsc{Section 4.1.3} of \cite{Fenzi2018}
% is not respected. 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 subplot(1,2,1), loglog(pp,aet(k,:),C,col{k},LW,2); 
 hold on; loglog(pp,net(k,:),'-.',C,col{k},LW,1.5);
 subplot(1,2,2), loglog(pp,aem(k,:),C,col{k},LW,2); 
 hold on; loglog(pp,nem(k,:),'-.',C,col{k},LW,1.5); 
end
subplot(1,2,1),grid on, xlabel(P2); title(TD); axis(AX)
subplot(1,2,2),grid on, xlabel(P2); title(MD); axis(AX)


%%
% <latex> 
% \subsection{Collocation approach}\label{sec:col2} \vspace{-1em}
% </latex>


%%
% <latex>
% The collocation approach in the bivariate case can be easily achieved
% by Chebfun. Let us consider the error of the interpolant on Padua
% points for total multivariate degree.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
err_t=zeros(3,Pdmax);
for pd=1:Pdmax,  XY=paduapts(pd);   
 for k=1:3,      sa=zeros(size(XY(:,1)));
  for i=1:length(sa)
   [A,B]=DelayedOscillator(Lt1(XY(i,1)),Lt2(XY(i,2)),k,IG,stru);
   sa(i)=max(real(eig(A,B)));
  end
  SA=chebfun2(sa, [-1 1 -1 1], 'padua');
  err_t(k,pd)=max(max(abs(SA(XXT,YYT)-AlphaT(:,:,k))));
 end
end



%%
% <latex>
% We analyze the error of the maximal degree polynomial interpolant
% on tensor product Chebyshev grid. (The construction of the grid differs 
% from the one previously seen.) 
% </latex>
%%
% <latex> \vspace{-2em} </latex>
err_m=zeros(3,Pdmax+1);
for pd=0:Pdmax, [XX, YY] = chebpts2(pd+1);
 for k=1:3,     sa=zeros((pd+1),(pd+1));
  for i=1:pd+1
   for j=1:pd+1
   [A,B]=DelayedOscillator(Lt1(XX(i,j)),Lt2(YY(i,j)),k,IG,stru);
   sa(i,j)=max(real(eig(A,B)));
   end
  end
  SA=chebfun2(sa);
  err_m(k,pd+1)=max(max(abs(SA(XXT,YYT)-AlphaT(:,:,k))));
 end
end, close





%%
% <latex> 
% The convergence of the errors can be shown by the following code.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
for k=1:3
 subplot(1,2,1),loglog(2:Pdmax+1,err_t(k,:),C,col{k},LW,2);hold on
 subplot(1,2,2),loglog(1:Pdmax+1,err_m(k,:),C,col{k},LW,2);hold on
end
subplot(1,2,1),grid on, xlabel(P2); ylabel(AP); title(TD); axis(AX)
subplot(1,2,2),grid on, xlabel(P2); ylabel(AP); title(MD); axis(AX)



%%
% <latex>
% The convergence rates are similar; in particular, 
% for the non-smooth cases, the
% empirically convergence rates are estimated by:
% </latex>
%%
% <latex> \vspace{-2em} </latex>
disp('        Total Degree     Maximal Degree');
for k=2:3
 ot=rates(err_t(k,:),8:2:Pdmax); om=rates(err_m(k,:),1:2:Pdmax);
 fprintf([cases{k},'  %10.8f %14.8f \n'], median(ot), median(om));
end

%%
% <latex> 
% \section*{Conclusions and Acknowledgments}\vspace{-1em}
% </latex>

%%
% <latex>
% Refer to \cite{Fenzi2018} for further analyses and conclusion.
% </latex>

%%
% <latex>
% \textbf{Acknowledgments.} This work was supported by the project 
% C14/17/072 of the KU Leuven Research Council, by the project G0A5317N of 
% the Research Foundation-Flanders (FWO - Vlaanderen), and by the project 
% UCoCoS, funded by the European Unions Horizon 2020 research and 
% innovation program under the Marie Sklodowska-Curie Grant Agreement 
% No 675080. 
% </latex>

%%
% <latex> 
% \addcontentsline{toc}{section}{References}
% \begin{thebibliography}{9}
% \bibitem{Abramowitz1965} 
% {\sc M.~Abramowitz \& I.~A.~Stegun,}
% {\em Handbook of Mathematical Functions,} 
% Dover, 1965.
% \bibitem{Breda2015} 
% {\sc D.~Breda, S.~Maset  \& R.~Vermiglio,}
% {\em Stability of Linear Delay
% Differential Equations - A numerical approach with MATLAB,} Springer, 2015.
% \bibitem{Boyd2001} 
% {\sc J.~P.~Boyd,} {\em Chebyshev and Fourier Spectral Methods,}
% Dover, 2001.
% \bibitem{Caliari2011}	
% {\sc M.~Caliari, S.~de Marchi, A.~Sommariva,  \& M.~Vianello,}
% Padua2DM: fast interpolation and cubature at the Padua points in
% MATLAB/Octave, {\em Numerical Algorithms}  56 (2011), 45--60.
% \bibitem{Fenzi2017} 
% {\sc L.~Fenzi  \& W.~Michiels,}  Robust stability optimization for linear
% delay systems in a probabilistic framework, {\em Linear Algebra and
% its Applications} 526 (2017), 1--26.
% \bibitem{Fenzi2018} 
% {\sc L.~Fenzi \& W.~Michiels,} Polynomial (chaos) approximation of
% maximum eigenvalue function: efficiency and limitations, {\em arXiv,}
% 1804.03881 (2018).
% \bibitem{Jarlebring2010} 
% {\sc E.~Jarlebring, K.~Meerbergen  \& W.~Michiels,} A Krilov method for
% the delay eigenvalue problem, {\em SIAM Journal of Computer Science,} 32 (2010),
% 3278--3300.
% \bibitem{Szudzik2017} 
% {\sc M.~P.~Szudzik,} The Rosenberg-Strong pairing function, {\em arXiv,}
% 1706.04129 (2017).
% \bibitem{Trefethen2000} 
% {\sc L.~N.~Trefethen,} {\em Spectral Methods in MATLAB,} SIAM, 2000.
% \bibitem{Trefethen2012} 
% {\sc L.~N.~Trefethen,} {\em Approximation Theory and Approximation Practice,}
% Cambridge University Press, 2012.
% \end{thebibliography}	
% </latex>

%%
% <latex> 
% \appendix
% \section{MATLAB functions for bivariate analysis}\label{appendix}\vspace{-1em}
% </latex>

%%
% <latex>
% This appendix reports the additional MATLAB functions requested for the
% analysis of the bivariate case, Section \ref{sec:d2}. Before each function, 
% the help is given.  
% </latex>


%%
% <latex>
% {\tt [A,B]=DelayedOscillator(nu,xi,RE,N,stru)} constructs the
% approximation of  the delay eigenvalue problem associated to 
% the oscillator with feedback delay, \textsc{Example 2} in \cite{Fenzi2018}.
% The non-linear finite dimensional eigenvalue problem is turned into a 
% discretized version of the infinite dimensional linear problem,
% by the Infinitesimal Generator approach \cite{Breda2015,Jarlebring2010}.
% The oscillator with feedback delay system is
% $$ \ddot{x}(t)=-{\tt nu}^2 x(t)-2{\tt nu}\cdot{\tt xi}\ \dot{x}(t)+ 
% {\tt K(1)}\ x(t-1)+{\tt K(2)}\ \dot{x}(t-1). $$
% \textbf{Input:}
% \begin{description} 
% \item[{\tt nu}] damping ratio, which is studied in $[0.9,1.1]$. 
% \item[{\tt xi}] angular frequency, which is studied in $[0.1,0.2]$. 
% \item[{\tt RE}] specifies the controller parameters and consequently the 
% behavior of the {\tt R}ightmost {\tt E}igenvalues, \textsc{Table 1} in \cite{Fenzi2018}. 
% The SAE, MSSAEs and MNSSAEs cases can be obtained by setting {\tt RE}
% equals to the numbers {\tt 1}, {\tt 2}, and {\tt 3}, respectively.
% \item[{\tt N}] discretization of the Infinitesimal Generator. It is an 
% integer which defines the dimension of the final characteristic matrix,
% \textit{i.e.}~$2({\tt N}+1)\times2({\tt N}+1)$. 
% \item[{\tt stru}] structure of the discretized infinitesimal generator. 
% Setting {\tt stru=0}, the eigenvalue problem does not present a particular
% structure as in \cite{Breda2015}, while setting {\tt stru=1}, we get a 
% structured eigenvalue problem by the approach \cite{Jarlebring2010}.
% \end{description}
% \textbf{Output:} {\tt (A,B)} generalized eigenvalue problem
% (cf.~\textit{e.g.}~{\tt help eig} for further information).\par
% </latex>
%%
% <latex> \vspace{-2em} </latex>
function [A,B]=DelayedOscillator(nu,xi,RE,N,stru)
% Controller parameters - (Table 1 in [6])
if RE==1,      K=[0.2, 0.2];          % SAE
elseif RE==2,  K=[0.5105, -0.0918];   % MSSAEs
elseif RE==3,  K=[0.6179, -0.0072];   % MNSSAEs
else, error('MyComponent:incorrectType', ['Error.\n'...
 'RE must be a number\n 1 - SAE,\n 2 - MSSAEs,\n 3 - MNSSAEs.'])
end
% Definition of the oscillator with feedback delay system
n=2;            % Dimension of the system
h=2; TAU=[0,1]; % Delays 
E=eye(n);       % Leading matrix 
A=cell(1,h); A{1}=[0,1; -nu^2, -2*nu*xi];  
             A{h}=[0,0;  K(1), K(2)]; 
% INFINITESIMAL GENERATOR APPROACH
if stru==0      % No stucture, approach in [2]
% Differentiation matrix D (pag 54 in [9])
  x=(TAU(h)/2)*(cos(pi*(0:N)/N)'-1); % Chebyshev nodes in [tau,0]
  c=[2; ones(N-1,1); 2].*(-1).^(0:N)';
  X=repmat(x,1,N+1);   dX=X-X';
  D=(c*(1./c)')./(dX+(eye(N+1))); D=D-diag(sum(D,2));
  DN=kron(D,eye(n));
% Definition of the eigenvalue problem (Example 5.1 in [2])   
  AN=zeros(n*(N+1));  AN(1:n,1:n)=A{1}; AN(1:n,n*N+1:end)=A{2};
  AN(n+1:end,:)=DN(n+1:end,:); A=AN; B=eye(size(A));
elseif stru==1    % Structured eigenvalue problem [7]
% Band Matrix L_N (Section 2.3 in [7])  
  L_N=zeros(N+1); L_N(2,1:3)=[2 0 -1];
  for i=3:1:N, L_N(i,i-1:i+1)=[1/(i-1) 0 -1/(i-1)];
  end
  L_N(N+1,N:N+1)=[1/N 0]; L_N=(TAU(h)/4)*L_N;
% Pi_N & Sigma_N matrices (Theorem 2.1 in [7].)
  Pi_N=kron(L_N,eye(n)); Pi_N(1:n,:)=kron(ones(1,N+1),E);
  Sigma_N=eye(n*(N+1));  Sigma_N(1:n,1:n)=zeros(n); 
  for i=1:N+1, xx=(-2*TAU/TAU(h)+1); 
 % Chebishev polynomial of degree i-1 evaluated in xx
   if i==1, CS=ones(1,h);
   elseif i==2, CS=xx;  C1=xx; C2=ones(1,h);       
   else,        CS=2*xx.*C1-C2; C2=C1; C1=CS;
   end
   for j=1:h % First row block of Sigma_N
      Sigma_N(1:n,((i-1)*n+1):(i*n))=...
         Sigma_N(1:n,((i-1)*n+1):(i*n))+A{j}*CS(j);
   end
  end, A=Sigma_N; B=Pi_N; % Definition of the eigenvalue problem
else
 error('MyComponent:incorrectType',['Error.\nstru must be a',...
  'number\n  0 - No Structure,\n  1 - Structured eigenproblem.'])   
end
end

%%
% <latex>
% {\tt W= pdwtsMM(n)} computes the cubature weights {\tt W} by Matrix
% Multiplication (MM) so that, if {\tt Pad} is the matrix of Padua points 
% computed in Chebfun by {\tt Pad = paduadpts(n)},
% then the cubature of a function {\tt funct} is given by 
% {\tt W'*funct(Pad(:,1),Pad(:,2))}.
% The interested reader is referred to \cite{Caliari2011} for further 
% information on the topic and on the algorithm. \par 
% \textbf{Input:} {\tt n} interpolation degree. \par 
% \textbf{Output:} {\tt W} cubature weights associated to the matrix {\tt
% Pad} of Padua points.
% </latex>
%%
% <latex> \vspace{-2em} </latex>
function W= pdwtsMM(n)
if n == 0,  W = 4; % degree 0
else
  argn1=linspace(0,pi,n+1); argn2=linspace(0,pi,n+2);
  k=(0:2:n)'; l=(n-mod(n,2))/2+1; lp=(n+mod(n,2))/2+1;
% even-degree Chebyshev polynomials on the subgrids
  TE1=cos(k*argn1(1:2:n+1)); TE1(2:l,:)=TE1(2:l,:)*sqrt(2);
  TO1=cos(k*argn1(2:2:n+1)); TO1(2:l,:)=TO1(2:l,:)*sqrt(2);
  TE2=cos(k*argn2(1:2:n+2)); TE2(2:l,:)=TE2(2:l,:)*sqrt(2);
  TO2=cos(k*argn2(2:2:n+2)); TO2(2:l,:)=TO2(2:l,:)*sqrt(2);
% even,even moments matrix
  mom=2*sqrt(2)./(1-k.^2); mom(1)=2; [M1,M2]=meshgrid(mom);
  M0 = fliplr(triu(fliplr(M1.*M2)));
% interpolation weights matrices
  W1=2*ones(l)/(n*(n+1));       W1(:,1)=W1(:,1)/2; 
  W2=2*ones(lp,lp-1)/(n*(n+1)); W2(1,:)=W2(1,:)/2;
  if mod(n,2)==0, i=n/2+1;
    M0(i,1) = M0(i,1)/2; W1(:,i)=W1(:,i)/2; W1(i,:)=W1(i,:)/2;
  else, i=(n+1)/2;
    W2(i+1,:)=W2(i+1,:)/2; W2(:,i)=W2(:,i)/2;
  end
% cubature weights as matrices on the subgrids.
  L1=W1.*(TE1'*M0*TO2)';  L2=W2.*(TO1'*M0*TE2)';
  if mod(n,2) == 0 % W=zeros(n/2+1,n+1);
    W(:,1:2:n+1)=L1; W(:,2:2:n+1)=L2; W=W(:);
  else
    W=[L1',L2']';                     W=W(:);
  end
end
end


