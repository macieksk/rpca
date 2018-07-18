\encoding{UTF-8}
\name{trpca}
\alias{trpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO Decompose a matrix into a low-rank component and a sparse component by solving Principal Components Pursuit
}
\description{
TODO This function decomposes a rectangular matrix \var{M} into a low-rank component, and a sparse component, by solving a convex program called Principal Component Pursuit.
}
\usage{
trpca(M,     #k,
             k.start=1,
             lambda = 1/sqrt(max(dim(M))), #This is ok only for dense matrices
             lambda2 = 100*lambda, #TODO needs proper L1 sparse vs L2 noise weight setting
             L2noise = TRUE, #Do decompose into M=L+S+E, or just M=L+S if FALSE
             mu = prod( dim(M)) / (4*sum(abs(M)) ), #This is ok only for dense matrices
             mu.max = mu*100,  #Stops mu from getting to large too fast 
                               #(i.e. from caring too much for constraint than objective.function)
             mu.min = mu/200,   #If smallest computed SV is larger than 1/mu.min we increase k.current
                               #and compute one more SV in next iteration
             mu.growth.ratio=1.1,
             term.delta=10^(-7),
             max.iter=5000,
             trace=FALSE,
             message.iter=100,
             n.iter.without.L2noise=5, #Number of start iterations without decomposing L2noise
             #thresh.nuclear.fun=trpca.thresh.nuclear.sparse2,
             #thresh.l1.fun=thresh.l1.sparse,
             #zero.matrix.fun=zero.matrix.sparse,             
             thresh.nuclear.fun=trpca.thresh.nuclear,
             thresh.l1.fun=thresh.l1,
             zero.matrix.fun=zero.matrix,                
             F2norm.fun=F2norm)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{M}{ a rectangular matrix that is to be decomposed into a low-rank component and a sparse component
           \eqn{M = L + S} .
}
  \item{lambda}{
   parameter of the convex problem \eqn{\|L\|_{*} + \lambda \|S\|_{1} }{||L||_{*} + lambda ||S||_{1} }
   which is minimized in the Principal Components Pursuit algorithm. 
   The default value is the one suggested in Candès, E. J., section 1.4, 
   and together with reasonable assumptions about \var{L} and \var{S} 
   guarantees that a correct decomposition is obtained.
}
  \item{mu}{
    parameter from the augumented Lagrange multiplier formulation of the PCP, Candès, E. J., section 5. 
    Default value is the one suggested in references.
}
  \item{term.delta}{
   The algorithm terminates when \eqn{\|M-L-S\|_{F} \leq \delta \|M\|_{F}}{||M-L-S||_F<=delta||M||_F} 
   where \eqn{\|\ \|_{F}}{|| ||_F} is Frobenius norm of a matrix. 
   %The default value is lower for GPU version, as it operates on single precision floats.
}
  \item{max.iter}{
    Maximal number of iterations of the augumented Lagrange multiplier algorithm. 
    A warning is issued if the algorithm does not converge by then.
}
  \item{trace}{
    Print out information with every iteration. 
}
%  \item{gpu.to.choose}{
%   If not \code{NULL}, the argument to \code{chooseGpu} function from package \code{gputools}.
%}
  \item{thresh.nuclear.fun, thresh.l1.fun, F2norm.fun}{
   Arguments for internal use only.
}
}
\details{
TODO, documentation for original rpca below.

These functions decompose a rectangular matrix \var{M} into a low-rank component, and a sparse component, by solving a convex program called Principal Component Pursuit:
    \deqn{ \textrm{minimize}\quad   \|L\|_{*} + \lambda \|S\|_{1} }{%
            minimize   ||L||_{*} + lambda ||S||_1 }                                                                    
    \deqn{ \textrm{subject to}\quad   L+S = M }{%
            subject to    L + S = M }
where \eqn{\|L\|_{*}}{||L||_{*}} is the nuclear norm of \var{L} (sum of singular values).
}
\value{
The function returns two matrices \code{S} and \code{L}, which have the property that 
\eqn{L+S \simeq M}{L + S ~= M}, where the quality of the approximation depends on the argument \code{term.delta},
and the convergence of the algorithm.
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{S }{The sparse component of the matrix decomposition.}
 \item{L }{The low-rank component of the matrix decomposition.}
 \item{L.svd}{The singular value decomposition of \code{L}, as returned by the function \code{La.svd} .}
 \item{convergence$converged}{\code{TRUE} if the algorithm converged with respect to \code{term.delta}.}
 \item{convergence$iterations}{Number of performed iterations.}
 \item{convergence$final.delta}{The final iteration \code{delta} which is compared with \code{term.delta}.}
 \item{convergence$all.delta}{All \code{delta} from all iterations.}
}
\references{
Candès, E. J., Li, X., Ma, Y., & Wright, J. (2011). Robust principal component analysis?. Journal of the ACM (JACM), 58(3), 11.

Yuan, X., & Yang, J. (2009). Sparse and low-rank matrix decomposition via alternating direction methods. preprint, 12.
}
\author{
\packageAuthor{rpca}
}
%\note{
%%%  ~~further notes~~
%\code{rpca.gpu} requires package \code{gputools} in version which provides function \code{gpuSvd}, and version 0.26 %does provide it when compiled with CULA library. This can be achieved by downloading the required version and %installing it with the command: \cr
%\code{install.packages("gputools_0.26.tar.gz",}\cr
%\code{  configure.args="--with-cuda-home=/opt/cuda --with-cula-home=/opt/cula")}
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ rpca }
\keyword{ trpca }
\keyword{ robust pca }
\keyword{ low-rank and sparse components }
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## TODO original rpca examples below

data(iris)
M <- as.matrix(iris[,1:4])
Mcent <- sweep(M,2,colMeans(M))

res <- rpca(Mcent)

## Check convergence and number of iterations
with(res$convergence,list(converged,iterations))
## Final delta F2 norm divided by F2norm(Mcent)
with(res$convergence,final.delta)

## Check properites of the decomposition
with(res,c(
all(abs( L+S - Mcent ) < 10^-5),
all( L == L.svd$u\%*\%(L.svd$d*L.svd$vt) )
))
# [1] TRUE TRUE

## The low rank component has rank 2
length(res$L.svd$d)
## However, the sparse component is not sparse 
## - thus this data set is not the best example here.
mean(res$S==0)

## Plot the first (the only) two principal components
## of the low-rank component L
rpc<-res$L.svd$u\%*\%diag(res$L.svd$d)
plot(jitter(rpc[,1:2],amount=.001),col=iris[,5])

## Compare with classical principal components
pc <- prcomp(M,center=TRUE)
plot(pc$x[,1:2],col=iris[,5])
points(rpc[,1:2],col=iris[,5],pch="+")

## "Sparse" elements distribution
plot(density(abs(res$S),from=0))
curve(dexp(x,rate=1/mean(abs(res$S))),add=TRUE,lty=2)

## Plot measurements against measurements corrected by sparse components
par(mfcol=c(2,2))
for(i in 1:4) {
plot(M[,i],M[,i]-res$S[,i],col=iris[,5],xlab=colnames(M)[i])
}
}
