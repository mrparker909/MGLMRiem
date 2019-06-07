#' @export
randspd <- function(n, c=3, udist=3) {
# RANDSPD generates n by b random symmatrix positive matrix P.
#
#   P = randspd(n)
#   P = randspd(n,c)
#   P = randspd(n,c,udit)
#
#   c is parameter for variance. Bigger c has bigger variance.
#   udist is the upper bound of distance from I to P w.r.t. GL-invariant
#   measure.
#
#   See also SYNTH_DTI_DATA

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $

  P = c*matrix(runif(n*n), ncol=n) # c*(rand(n)-0.5);
  P = P%*%t(P)

  while(dist_M_spd(P,diag(rep(1,n))) > udist) {
    P = c*matrix(runif(n*n), ncol=n) # c*(rand(n)-0.5);
    P = P%*%t(P)
  }
  return(P)
}
