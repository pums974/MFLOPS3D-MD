      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
     * iwrk,kwrk,ier)
c  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   x     : real array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : real array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1)+my*(ky+1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
c           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  if(my-1) 100,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
  50  continue
  60  ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),
     * iwrk(1),iwrk(mx+1))
 100  return
      end
      subroutine clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the ordered set of m points x(i) in the idim-dimensional space
c  with x(1)=x(m), and given also a corresponding set of strictly in-
c  creasing values u(i) and the set of positive numbers w(i),i=1,2,...,m
c  subroutine clocur determines a smooth approximating closed spline
c  curve s(u), i.e.
c      x1 = s1(u)
c      x2 = s2(u)       u(1) <= u <= u(m)
c      .........
c      xidim = sidim(u)
c  with sj(u),j=1,2,...,idim periodic spline functions of degree k with
c  common knots t(j),j=1,2,...,n.
c  if ipar=1 the values u(i),i=1,2,...,m must be supplied by the user.
c  if ipar=0 these values are chosen automatically by clocur as
c      v(1) = 0
c      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
c      u(i) = v(i)/v(m) ,i=1,2,...,m
c  if iopt=-1 clocur calculates the weighted least-squares closed spline
c  curve according to a given set of knots.
c  if iopt>=0 the number of knots of the splines sj(u) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(u) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(u) is given in the b-spline representation and can be
c  evaluated by means of subroutine curev.
c
c  calling sequence:
c     call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,
c    * fp,wrk,lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares closed spline curve (iopt=-1) or a smoothing
c           closed spline curve (iopt=0 or 1) must be determined. if
c           iopt=0 the routine will start with an initial set of knots
c           t(i)=u(1)+(u(m)-u(1))*(i-k-1),i=1,2,...,2*k+2. if iopt=1 the
c           routine will continue with the knots found at the last call.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   ipar  : integer flag. on entry ipar must specify whether (ipar=1)
c           the user will supply the parameter values u(i),or whether
c           (ipar=0) these values are to be calculated by clocur.
c           unchanged on exit.
c   idim  : integer. on entry idim must specify the dimension of the
c           curve. 0 < idim < 11.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > 1. unchanged on exit.
c   u     : real array of dimension at least (m). in case ipar=1,before
c           entry, u(i) must be set to the i-th value of the parameter
c           variable u for i=1,2,...,m. these values must then be
c           supplied in strictly ascending order and will be unchanged
c           on exit. in case ipar=0, on exit,the array will contain the
c           values u(i) as determined by clocur.
c   mx    : integer. on entry mx must specify the actual dimension of
c           the array x as declared in the calling (sub)program. mx must
c           not be too small (see x). unchanged on exit.
c   x     : real array of dimension at least idim*m.
c           before entry, x(idim*(i-1)+j) must contain the j-th coord-
c           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
c           idim. since first and last data point must coincide it
c           means that x(j)=x(idim*(m-1)+j),j=1,2,...,idim.
c           unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. w(m) is not used.
c           unchanged on exit. see also further comments.
c   k     : integer. on entry k must specify the degree of the splines.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the splines returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+2*k, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the smoothing spline curve returned
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline curve,i.e. the position of the interior knots t(k+2),
c           t(k+3),..,t(n-k-1) as well as the position of the additional
c           t(1),t(2),..,t(k+1)=u(1) and u(m)=t(n-k),...,t(n) needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   nc    : integer. on entry nc must specify the actual dimension of
c           the array c as declared in the calling (sub)program. nc
c           must not be too small (see c). unchanged on exit.
c   c     : real array of dimension at least (nest*idim).
c           on succesful exit, this array will contain the coefficients
c           in the b-spline representation of the spline curve s(u),i.e.
c           the b-spline coefficients of the spline sj(u) will be given
c           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
c   fp    : real. unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline curve returned.
c   wrk   : real array of dimension at least m*(k+1)+nest*(7+idim+5*k).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the close curve returned has a residual
c             sum of squares fp such that abs(fp-s)/s <= tol with tol a
c             relative tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the curve returned is an interpolating
c             spline curve (fp=0).
c    ier=-2 : normal return. the curve returned is the weighted least-
c             squares point,i.e. each spline sj(u) is a constant. in
c             this extreme case fp gives the upper bound fp0 for the
c             smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares closed
c             curve according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing curve with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing curve
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,2,...,m
c             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(7+idim+5*k),
c             nc>=nest*idim, x(j)=x(idim*(m-1)+j), j=1,2,...,idim
c             if ipar=0: sum j=1,idim (x(i*idim+j)-x((i-1)*idim+j))**2>0
c                        i=1,2,...,m-1.
c             if ipar=1: u(1)<u(2)<...<u(m)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
c                         u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
c                            (u(1)=0 and u(m)=1 in case ipar=0)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points uu(j) with
c                       uu(j) = u(i) or u(i)+(u(m)-u(1)) such that
c                         t(j) < uu(j) < t(j+k+1), j=k+1,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+2*k
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the curve will be too smooth and signal will be
c   lost ; if s is too small the curve will pick up too much noise. in
c   the extreme cases the program will return an interpolating curve if
c   s=0 and the weighted least-squares point if s is very large.
c   between these extremes, a properly chosen s will result in a good
c   compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in x(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the weighted
c   least-squares point and the upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximating curve shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if clocur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   curve underlying the data. but, if the computation mode iopt=1 is
c   used, the knots returned may also depend on the s-values at previous
c   calls (if these were smaller). therefore, if after a number of
c   trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   clocur once more with the selected value for s but now with iopt=0.
c   indeed, clocur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c   the form of the approximating curve can strongly be affected  by
c   the choice of the parameter values u(i). if there is no physical
c   reason for choosing a particular parameter u, often good results
c   will be obtained with the choice of clocur(in case ipar=0), i.e.
c        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
c   where
c        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
c   other possibilities for q(i) are
c        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
c        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= 1
c
c
c  other subroutines required:
c    fpbacp,fpbspl,fpchep,fpclos,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real s,fp
      integer iopt,ipar,idim,m,mx,k,nest,n,nc,lwrk,ier
c  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real per,tol,dist
      integer i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,lwest,
     * maxit,m1,nmin,ncc,j
c  ..function references..
      real sqrt
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(ipar.lt.0 .or. ipar.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      k2 = k1+1
      nmin = 2*k1
      if(m.lt.2 .or. nest.lt.nmin) go to 90
      ncc = nest*idim
      if(mx.lt.m*idim .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(7+idim+5*k)
      if(lwrk.lt.lwest) go to 90
      i1 = idim
      i2 = m*idim
      do 5 j=1,idim
         if(x(i1).ne.x(i2)) go to 90
         i1 = i1-1
         i2 = i2-1
   5  continue
      if(ipar.ne.0 .or. iopt.gt.0) go to 40
      i1 = 0
      i2 = idim
      u(1) = 0.
      do 20 i=2,m
         dist = 0.
         do 10 j1=1,idim
            i1 = i1+1
            i2 = i2+1
            dist = dist+(x(i2)-x(i1))**2
  10     continue
         u(i) = u(i-1)+sqrt(dist)
  20  continue
      if(u(m).le.0.) go to 90
      do 30 i=2,m
         u(i) = u(i)/u(m)
  30  continue
      u(m) = 0.1e+01
  40  if(w(1).le.0.) go to 90
      m1 = m-1
      do 50 i=1,m1
         if(u(i).ge.u(i+1) .or. w(i).le.0.) go to 90
  50  continue
      if(iopt.ge.0) go to 70
      if(n.le.nmin .or. n.gt.nest) go to 90
      per = u(m)-u(1)
      j1 = k1
      t(j1) = u(1)
      i1 = n-k
      t(i1) = u(m)
      j2 = j1
      i2 = i1
      do 60 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  60  continue
      call fpchep(u,m,t,n,k,ier)
      if(ier) 90,80,90
  70  if(s.lt.0.) go to 90
      if(s.eq.0. .and. nest.lt.(m+2*k)) go to 90
      ier = 0
c we partition the working space and determine the spline approximation.
  80  ifp = 1
      iz = ifp+nest
      ia1 = iz+ncc
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpclos(iopt,idim,m,u,mx,x,w,k,s,nest,tol,maxit,k1,k2,n,t,
     * ncc,c,fp,wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1),
     * wrk(ig2),wrk(iq),iwrk,ier)
  90  return
      end
      subroutine cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,
     * lwrk,iwrk,kwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m, subroutine cocosp determines the weighted
c  least-squares cubic spline s(x) with given knots t(j),j=1,2,...,n
c  which satisfies the following concavity/convexity conditions
c      s''(t(j+3))*e(j) <= 0, j=1,2,...n-6
c  the fit is given in the b-spline representation( b-spline coef-
c  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,
c    * lwrk,iwrk,kwrk,ier)
c
c  parameters:
c    m   : integer. on entry m must specify the number of data points.
c          m > 3. unchanged on exit.
c    x   : real array of dimension at least (m). before entry, x(i)
c          must be set to the i-th value of the independent variable x,
c          for i=1,2,...,m. these values must be supplied in strictly
c          ascending order. unchanged on exit.
c    y   : real array of dimension at least (m). before entry, y(i)
c          must be set to the i-th value of the dependent variable y,
c          for i=1,2,...,m. unchanged on exit.
c    w   : real array of dimension at least (m). before entry, w(i)
c          must be set to the i-th value in the set of weights. the
c          w(i) must be strictly positive. unchanged on exit.
c    n   : integer. on entry n must contain the total number of knots
c          of the cubic spline. m+4>=n>=8. unchanged on exit.
c    t   : real array of dimension at least (n). before entry, this
c          array must contain the knots of the spline, i.e. the position
c          of the interior knots t(5),t(6),...,t(n-4) as well as the
c          position of the boundary knots t(1),t(2),t(3),t(4) and t(n-3)
c          t(n-2),t(n-1),t(n) needed for the b-spline representation.
c          unchanged on exit. see also the restrictions (ier=10).
c    e   : real array of dimension at least (n). before entry, e(j)
c          must be set to 1 if s(x) must be locally concave at t(j+3),
c          to (-1) if s(x) must be locally convex at t(j+3) and to 0
c          if no convexity constraint is imposed at t(j+3),j=1,2,..,n-6.
c          e(n-5),...,e(n) are not used. unchanged on exit.
c  maxtr : integer. on entry maxtr must contain an over-estimate of the
c          total number of records in the used tree structure, to indic-
c          ate the storage space available to the routine. maxtr >=1
c          in most practical situation maxtr=100 will be sufficient.
c          always large enough is
c                         n-5       n-6
c              maxtr =  (     ) + (     )  with l the greatest
c                          l        l+1
c          integer <= (n-6)/2 . unchanged on exit.
c  maxbin: integer. on entry maxbin must contain an over-estimate of the
c          number of knots where s(x) will have a zero second derivative
c          maxbin >=1. in most practical situation maxbin = 10 will be
c          sufficient. always large enough is maxbin=n-6.
c          unchanged on exit.
c    c   : real array of dimension at least (n).
c          on succesful exit, this array will contain the coefficients
c          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
c    sq  : real. on succesful exit, sq contains the weighted sum of
c          squared residuals of the spline approximation returned.
c    sx  : real array of dimension at least m. on succesful exit
c          this array will contain the spline values s(x(i)),i=1,...,m
c   bind : logical array of dimension at least (n). on succesful exit
c          this array will indicate the knots where s''(x)=0, i.e.
c                s''(t(j+3)) .eq. 0 if  bind(j) = .true.
c                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6
c   wrk  : real array of dimension at least  m*4+n*7+maxbin*(maxbin+n+1)
c          used as working space.
c   lwrk : integer. on entry,lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.lwrk
c          must not be too small (see wrk). unchanged on exit.
c   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
c          used as working space.
c   kwrk : integer. on entry,kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program. kwrk
c          must not be too small (see iwrk). unchanged on exit.
c   ier   : integer. error flag
c      ier=0 : succesful exit.
c      ier>0 : abnormal termination: no approximation is returned
c        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
c                 probably causes : maxbin too small.
c        ier=2  : the number of records in the tree structure exceeds
c                 maxtr.
c                 probably causes : maxtr too small.
c        ier=3  : the algoritm finds no solution to the posed quadratic
c                 programming problem.
c                 probably causes : rounding errors.
c        ier=10 : on entry, the input data are controlled on validity.
c                 the following restrictions must be satisfied
c                   m>3, maxtr>=1, maxbin>=1, 8<=n<=m+4,w(i) > 0,
c                   x(1)<x(2)<...<x(m), t(1)<=t(2)<=t(3)<=t(4)<=x(1),
c                   x(1)<t(5)<t(6)<...<t(n-4)<x(m)<=t(n-3)<=...<=t(n),
c                   kwrk>=maxtr*4+2*(maxbin+1),
c                   lwrk>=m*4+n*7+maxbin*(maxbin+n+1),
c                   the schoenberg-whitney conditions, i.e. there must
c                   be a subset of data points xx(j) such that
c                     t(j) < xx(j) < t(j+4), j=1,2,...,n-4
c                 if one of these restrictions is found to be violated
c                 control is immediately repassed to the calling program
c
c
c  other subroutines required:
c    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno,fpchec
c
c  references:
c   dierckx p. : an algorithm for cubic spline fitting with convexity
c                constraints, computing 24 (1980) 349-371.
c   dierckx p. : an algorithm for least-squares cubic spline fitting
c                with convexity and concavity constraints, report tw39,
c                dept. computer science, k.u.leuven, 1978.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c   p. dierckx
c   dept. computer science, k.u.leuven
c   celestijnenlaan 200a, b-3001 heverlee, belgium.
c   e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : march 1978
c  latest update : march 1987.
c
c  ..
c  ..scalar arguments..
      real sq
      integer m,n,maxtr,maxbin,lwrk,kwrk,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(n),e(n),c(n),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(n)
c  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,ji,jib,jjb,jl,jr,ju,kwest,
     * lwest,mb,nm,n6
      real one
c  ..
c  set constant
      one = 0.1e+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(m.lt.4 .or. n.lt.8) go to 40
      if(maxtr.lt.1 .or. maxbin.lt.1) go to 40
      lwest = 7*n+m*4+maxbin*(1+n+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 40
      if(w(1).le.0.) go to 40
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 40
  10  continue
      call fpchec(x,m,t,n,3,ier)
      if(ier) 40,20,40
c  set numbers e(i)
  20  n6 = n-6
      do 30 i=1,n6
        if(e(i).gt.0.) e(i) = one
        if(e(i).lt.0.) e(i) = -one
  30  continue
c  we partition the working space and determine the spline approximation
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia),
     * wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji),
     * iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
  40  return
      end
      subroutine concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,
     * sx,bind,wrk,lwrk,iwrk,kwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i), i=1,2,...,m,subroutine concon determines a cubic spline
c  approximation s(x) which satisfies the following local convexity
c  constraints  s''(x(i))*v(i) <= 0, i=1,2,...,m.
c  the number of knots n and the position t(j),j=1,2,...n is chosen
c  automatically by the routine in a way that
c       sq = sum((w(i)*(y(i)-s(x(i))))**2) be <= s.
c  the fit is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...n-4) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c
c     call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,
c    * sx,bind,wrk,lwrk,iwrk,kwrk,ier)
c
c  parameters:
c    iopt: integer flag.
c          if iopt=0, the routine will start with the minimal number of
c          knots to guarantee that the convexity conditions will be
c          satisfied. if iopt=1, the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately
c          preceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c    m   : integer. on entry m must specify the number of data points.
c          m > 3. unchanged on exit.
c    x   : real array of dimension at least (m). before entry, x(i)
c          must be set to the i-th value of the independent variable x,
c          for i=1,2,...,m. these values must be supplied in strictly
c          ascending order. unchanged on exit.
c    y   : real array of dimension at least (m). before entry, y(i)
c          must be set to the i-th value of the dependent variable y,
c          for i=1,2,...,m. unchanged on exit.
c    w   : real array of dimension at least (m). before entry, w(i)
c          must be set to the i-th value in the set of weights. the
c          w(i) must be strictly positive. unchanged on exit.
c    v   : real array of dimension at least (m). before entry, v(i)
c          must be set to 1 if s(x) must be locally concave at x(i),
c          to (-1) if s(x) must be locally convex at x(i) and to 0
c          if no convexity constraint is imposed at x(i).
c    s   : real. on entry s must specify an over-estimate for the
c          the weighted sum of squared residuals sq of the requested
c          spline. s >=0. unchanged on exit.
c   nest : integer. on entry nest must contain an over-estimate of the
c          total number of knots of the spline returned, to indicate
c          the storage space available to the routine. nest >=8.
c          in most practical situation nest=m/2 will be sufficient.
c          always large enough is  nest=m+4. unchanged on exit.
c  maxtr : integer. on entry maxtr must contain an over-estimate of the
c          total number of records in the used tree structure, to indic-
c          ate the storage space available to the routine. maxtr >=1
c          in most practical situation maxtr=100 will be sufficient.
c          always large enough is
c                         nest-5      nest-6
c              maxtr =  (       ) + (        )  with l the greatest
c                           l          l+1
c          integer <= (nest-6)/2 . unchanged on exit.
c  maxbin: integer. on entry maxbin must contain an over-estimate of the
c          number of knots where s(x) will have a zero second derivative
c          maxbin >=1. in most practical situation maxbin = 10 will be
c          sufficient. always large enough is maxbin=nest-6.
c          unchanged on exit.
c    n   : integer.
c          on exit with ier <=0, n will contain the total number of
c          knots of the spline approximation returned. if the comput-
c          ation mode iopt=1 is used this value of n should be left
c          unchanged between subsequent calls.
c    t   : real array of dimension at least (nest).
c          on exit with ier<=0, this array will contain the knots of the
c          spline,i.e. the position of the interior knots t(5),t(6),...,
c          t(n-4) as well as the position of the additional knots
c          t(1)=t(2)=t(3)=t(4)=x(1) and t(n-3)=t(n-2)=t(n-1)=t(n)=x(m)
c          needed for the the b-spline representation.
c          if the computation mode iopt=1 is used, the values of t(1),
c          t(2),...,t(n) should be left unchanged between subsequent
c          calls.
c    c   : real array of dimension at least (nest).
c          on succesful exit, this array will contain the coefficients
c          c(1),c(2),..,c(n-4) in the b-spline representation of s(x)
c    sq  : real. unless ier>0 , sq contains the weighted sum of
c          squared residuals of the spline approximation returned.
c    sx  : real array of dimension at least m. on exit with ier<=0
c          this array will contain the spline values s(x(i)),i=1,...,m
c          if the computation mode iopt=1 is used, the values of sx(1),
c          sx(2),...,sx(m) should be left unchanged between subsequent
c          calls.
c    bind: logical array of dimension at least nest. on exit with ier<=0
c          this array will indicate the knots where s''(x)=0, i.e.
c                s''(t(j+3)) .eq. 0 if  bind(j) = .true.
c                s''(t(j+3)) .ne. 0 if  bind(j) = .false., j=1,2,...,n-6
c          if the computation mode iopt=1 is used, the values of bind(1)
c          ,...,bind(n-6) should be left unchanged between subsequent
c          calls.
c   wrk  : real array of dimension at least (m*4+nest*8+maxbin*(maxbin+
c          nest+1)). used as working space.
c   lwrk : integer. on entry,lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.lwrk
c          must not be too small (see wrk). unchanged on exit.
c   iwrk : integer array of dimension at least (maxtr*4+2*(maxbin+1))
c          used as working space.
c   kwrk : integer. on entry,kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program. kwrk
c          must not be too small (see iwrk). unchanged on exit.
c   ier   : integer. error flag
c      ier=0 : normal return, s(x) satisfies the concavity/convexity
c              constraints and sq <= s.
c      ier<0 : abnormal termination: s(x) satisfies the concavity/
c              convexity constraints but sq > s.
c        ier=-3 : the requested storage space exceeds the available
c                 storage space as specified by the parameter nest.
c                 probably causes: nest too small. if nest is already
c                 large (say nest > m/2), it may also indicate that s
c                 is too small.
c                 the approximation returned is the least-squares cubic
c                 spline according to the knots t(1),...,t(n) (n=nest)
c                 which satisfies the convexity constraints.
c        ier=-2 : the maximal number of knots n=m+4 has been reached.
c                 probably causes: s too small.
c        ier=-1 : the number of knots n is less than the maximal number
c                 m+4 but concon finds that adding one or more knots
c                 will not further reduce the value of sq.
c                 probably causes : s too small.
c      ier>0 : abnormal termination: no approximation is returned
c        ier=1  : the number of knots where s''(x)=0 exceeds maxbin.
c                 probably causes : maxbin too small.
c        ier=2  : the number of records in the tree structure exceeds
c                 maxtr.
c                 probably causes : maxtr too small.
c        ier=3  : the algoritm finds no solution to the posed quadratic
c                 programming problem.
c                 probably causes : rounding errors.
c        ier=4  : the minimum number of knots (given by n) to guarantee
c                 that the concavity/convexity conditions will be
c                 satisfied is greater than nest.
c                 probably causes: nest too small.
c        ier=5  : the minimum number of knots (given by n) to guarantee
c                 that the concavity/convexity conditions will be
c                 satisfied is greater than m+4.
c                 probably causes: strongly alternating convexity and
c                 concavity conditions. normally the situation can be
c                 coped with by adding n-m-4 extra data points (found
c                 by linear interpolation e.g.) with a small weight w(i)
c                 and a v(i) number equal to zero.
c        ier=10 : on entry, the input data are controlled on validity.
c                 the following restrictions must be satisfied
c                   0<=iopt<=1, m>3, nest>=8, s>=0, maxtr>=1, maxbin>=1,
c                   kwrk>=maxtr*4+2*(maxbin+1), w(i)>0, x(i) < x(i+1),
c                   lwrk>=m*4+nest*8+maxbin*(maxbin+nest+1)
c                 if one of these restrictions is found to be violated
c                 control is immediately repassed to the calling program
c
c  further comments:
c    as an example of the use of the computation mode iopt=1, the
c    following program segment will cause concon to return control
c    each time a spline with a new set of knots has been computed.
c     .............
c     iopt = 0
c     s = 0.1e+60  (s very large)
c     do 10 i=1,m
c       call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
c    *  bind,wrk,lwrk,iwrk,kwrk,ier)
c       ......
c       s = sq
c       iopt=1
c 10  continue
c     .............
c
c  other subroutines required:
c    fpcoco,fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
c
c  references:
c   dierckx p. : an algorithm for cubic spline fitting with convexity
c                constraints, computing 24 (1980) 349-371.
c   dierckx p. : an algorithm for least-squares cubic spline fitting
c                with convexity and concavity constraints, report tw39,
c                dept. computer science, k.u.leuven, 1978.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c   p. dierckx
c   dept. computer science, k.u.leuven
c   celestijnenlaan 200a, b-3001 heverlee, belgium.
c   e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : march 1978
c  latest update : march 1987.
c
c  ..
c  ..scalar arguments..
      real s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
c  ..array arguments..
      real x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),wrk(lwrk)
      integer iwrk(kwrk)
      logical bind(nest)
c  ..local scalars..
      integer i,lwest,kwest,ie,iw,lww
      real one
c  ..
c  set constant
      one = 0.1e+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.0 .or. iopt.gt.1) go to 30
      if(m.lt.4 .or. nest.lt.8) go to 30
      if(s.lt.0.) go to 30
      if(maxtr.lt.1 .or. maxbin.lt.1) go to 30
      lwest = 8*nest+m*4+maxbin*(1+nest+maxbin)
      kwest = 4*maxtr+2*(maxbin+1)
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 30
      if(iopt.gt.0) go to 20
      if(w(1).le.0.) go to 30
      if(v(1).gt.0.) v(1) = one
      if(v(1).lt.0.) v(1) = -one
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 30
         if(v(i).gt.0.) v(i) = one
         if(v(i).lt.0.) v(i) = -one
  10  continue
  20  ier = 0
c  we partition the working space and determine the spline approximation
      ie = 1
      iw = ie+nest
      lww = lwrk-nest
      call fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
     * bind,wrk(ie),wrk(iw),lww,iwrk,kwrk,ier)
  30  return
      end
      subroutine concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,
     * nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
c  given the ordered set of m points x(i) in the idim-dimensional space
c  and given also a corresponding set of strictly increasing values u(i)
c  and the set of positive numbers w(i),i=1,2,...,m, subroutine concur
c  determines a smooth approximating spline curve s(u), i.e.
c      x1 = s1(u)
c      x2 = s2(u)      ub = u(1) <= u <= u(m) = ue
c      .........
c      xidim = sidim(u)
c  with sj(u),j=1,2,...,idim spline functions of odd degree k with
c  common knots t(j),j=1,2,...,n.
c  in addition these splines will satisfy the following boundary
c  constraints        (l)
c      if ib > 0 :  sj   (u(1)) = db(idim*l+j) ,l=0,1,...,ib-1
c  and                (l)
c      if ie > 0 :  sj   (u(m)) = de(idim*l+j) ,l=0,1,...,ie-1.
c  if iopt=-1 concur calculates the weighted least-squares spline curve
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the splines sj(u) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(u) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(u) is given in the b-spline representation and can be
c  evaluated by means of subroutine curev.
c
c  calling sequence:
c     call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,
c    * t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline curve (iopt=-1) or a smoothing spline
c           curve (iopt=0 or 1) must be determined.if iopt=0 the routine
c           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
c           i=1,2,...,k+1. if iopt=1 the routine will continue with the
c           knots found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   idim  : integer. on entry idim must specify the dimension of the
c           curve. 0 < idim < 11.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k-max(ib-1,0)-max(ie-1,0). unchanged on exit.
c   u     : real array of dimension at least (m). before entry,
c           u(i) must be set to the i-th value of the parameter variable
c           u for i=1,2,...,m. these values must be supplied in
c           strictly ascending order and will be unchanged on exit.
c   mx    : integer. on entry mx must specify the actual dimension of
c           the arrays x and xx as declared in the calling (sub)program
c           mx must not be too small (see x). unchanged on exit.
c   x     : real array of dimension at least idim*m.
c           before entry, x(idim*(i-1)+j) must contain the j-th coord-
c           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
c           idim. unchanged on exit.
c   xx    : real array of dimension at least idim*m.
c           used as working space. on exit xx contains the coordinates
c           of the data points to which a spline curve with zero deriv-
c           ative constraints has been determined.
c           if the computation mode iopt =1 is used xx should be left
c           unchanged between calls.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   ib    : integer. on entry ib must specify the number of derivative
c           constraints for the curve at the begin point. 0<=ib<=(k+1)/2
c           unchanged on exit.
c   db    : real array of dimension nb. before entry db(idim*l+j) must
c           contain the l-th order derivative of sj(u) at u=u(1) for
c           j=1,2,...,idim and l=0,1,...,ib-1 (if ib>0).
c           unchanged on exit.
c   nb    : integer, specifying the dimension of db. nb>=max(1,idim*ib)
c           unchanged on exit.
c   ie    : integer. on entry ie must specify the number of derivative
c           constraints for the curve at the end point. 0<=ie<=(k+1)/2
c           unchanged on exit.
c   de    : real array of dimension ne. before entry de(idim*l+j) must
c           contain the l-th order derivative of sj(u) at u=u(m) for
c           j=1,2,...,idim and l=0,1,...,ie-1 (if ie>0).
c           unchanged on exit.
c   ne    : integer, specifying the dimension of de. ne>=max(1,idim*ie)
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the splines.
c           k=1,3 or 5.
c           unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the splines returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+k+1+max(0,ib-1)+max(0,ie-1),
c           the number of knots needed for interpolation (s=0).
c           unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the smoothing spline curve returned
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline curve,i.e. the position of the interior knots t(k+2),
c           t(k+3),..,t(n-k-1) as well as the position of the additional
c           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   nc    : integer. on entry nc must specify the actual dimension of
c           the array c as declared in the calling (sub)program. nc
c           must not be too small (see c). unchanged on exit.
c   c     : real array of dimension at least (nest*idim).
c           on succesful exit, this array will contain the coefficients
c           in the b-spline representation of the spline curve s(u),i.e.
c           the b-spline coefficients of the spline sj(u) will be given
c           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
c   cp    : real array of dimension at least 2*(k+1)*idim.
c           on exit cp will contain the b-spline coefficients of a
c           polynomial curve which satisfies the boundary constraints.
c           if the computation mode iopt =1 is used cp should be left
c           unchanged between calls.
c   np    : integer. on entry np must specify the actual dimension of
c           the array cp as declared in the calling (sub)program. np
c           must not be too small (see cp). unchanged on exit.
c   fp    : real. unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline curve returned.
c   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the curve returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the curve returned is an interpolating
c             spline curve, satisfying the constraints (fp=0).
c    ier=-2 : normal return. the curve returned is the weighted least-
c             squares polynomial curve of degree k, satisfying the
c             constraints. in this extreme case fp gives the upper
c             bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares spline
c             curve according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline curve
c             with fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing curve
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, k = 1,3 or 5, m>k-max(0,ib-1)-max(0,ie-1),
c             nest>=2k+2, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
c             nc >=nest*idim ,u(1)<u(2)<...<u(m),w(i)>0 i=1,2,...,m,
c             mx>=idim*m,0<=ib<=(k+1)/2,0<=ie<=(k+1)/2,nb>=1,ne>=1,
c             nb>=ib*idim,ne>=ib*idim,np>=2*(k+1)*idim,
c             if iopt=-1:2*k+2<=n<=min(nest,mmax) with mmax = m+k+1+
c                        max(0,ib-1)+max(0,ie-1)
c                        u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points uu(j) such that
c                         t(j) < uu(j) < t(j+k+1), j=1+max(0,ib-1),...
c                                                   ,n+k-1-max(0,ie-1)
c             if iopt>=0: s>=0
c                         if s=0 : nest >=mmax (see above)
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the curve will be too smooth and signal will be
c   lost ; if s is too small the curve will pick up too much noise. in
c   the extreme cases the program will return an interpolating curve if
c   s=0 and the least-squares polynomial curve of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in x(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial curve and the upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximating curve shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if concur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   curve underlying the data. but, if the computation mode iopt=1 is
c   used, the knots returned may also depend on the s-values at previous
c   calls (if these were smaller). therefore, if after a number of
c   trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   concur once more with the selected value for s but now with iopt=0.
c   indeed, concur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c   the form of the approximating curve can strongly be affected by
c   the choice of the parameter values u(i). if there is no physical
c   reason for choosing a particular parameter u, often good results
c   will be obtained with the choice
c        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
c   where
c        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
c   other possibilities for q(i) are
c        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
c        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= 1
c
c  other subroutines required:
c    fpback,fpbspl,fpched,fpcons,fpdisc,fpgivs,fpknot,fprati,fprota
c    curev,fppocu,fpadpo,fpinst
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real s,fp
      integer iopt,idim,m,mx,ib,nb,ie,ne,k,nest,n,nc,np,lwrk,ier
c  ..array arguments..
      real u(m),x(mx),xx(mx),db(nb),de(ne),w(m),t(nest),c(nc),wrk(lwrk)
      real cp(np)
      integer iwrk(nest)
c  ..local scalars..
      real tol,dist
      integer i,ib1,ie1,ja,jb,jfp,jg,jq,jz,j,k1,k2,lwest,maxit,nmin,
     * ncc,kk,mmin,nmax,mxx
c ..function references
      integer max0
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      kk = k1/2
      if(kk*2.ne.k1) go to 90
      k2 = k1+1
      if(ib.lt.0 .or. ib.gt.kk) go to 90
      if(ie.lt.0 .or. ie.gt.kk) go to 90
      nmin = 2*k1
      ib1 = max0(0,ib-1)
      ie1 = max0(0,ie-1)
      mmin = k1-ib1-ie1
      if(m.lt.mmin .or. nest.lt.nmin) go to 90
      if(nb.lt.(idim*ib) .or. ne.lt.(idim*ie)) go to 90
      if(np.lt.(2*k1*idim)) go to 90
      mxx = m*idim
      ncc = nest*idim
      if(mx.lt.mxx .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if(lwrk.lt.lwest) go to 90
      if(w(1).le.0.) go to 90
      do 10 i=2,m
         if(u(i-1).ge.u(i) .or. w(i).le.0.) go to 90
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 90
      j = n
      do 20 i=1,k1
         t(i) = u(1)
         t(j) = u(m)
         j = j-1
  20  continue
      call fpched(u,m,t,n,k,ib,ie,ier)
      if(ier) 90,40,90
  30  if(s.lt.0.) go to 90
      nmax = m+k1+ib1+ie1
      if(s.eq.0. .and. nest.lt.nmax) go to 90
      ier = 0
      if(iopt.gt.0) go to 70
c  we determine a polynomial curve satisfying the boundary constraints.
  40  call fppocu(idim,k,u(1),u(m),ib,db,nb,ie,de,ne,cp,np)
c  we generate new data points which will be approximated by a spline
c  with zero derivative constraints.
      j = nmin
      do 50 i=1,k1
        wrk(i) = u(1)
        wrk(j) = u(m)
        j = j-1
  50  continue
c  evaluate the polynomial curve
      call curev(idim,wrk,nmin,cp,np,k,u,m,xx,mxx,ier)
c  substract from the old data, the values of the polynomial curve
      do 60 i=1,mxx
        xx(i) = x(i)-xx(i)
  60  continue
c we partition the working space and determine the spline curve.
  70  jfp = 1
      jz = jfp+nest
      ja = jz+ncc
      jb = ja+nest*k1
      jg = jb+nest*k2
      jq = jg+nest*k2
      call fpcons(iopt,idim,m,u,mxx,xx,w,ib,ie,k,s,nest,tol,maxit,k1,
     * k2,n,t,ncc,c,fp,wrk(jfp),wrk(jz),wrk(ja),wrk(jb),wrk(jg),wrk(jq),
     * iwrk,ier)
c  add the polynomial curve to the calculated spline.
      call fpadpo(idim,t,n,c,ncc,k,cp,np,wrk(jz),wrk(ja),wrk(jb))
  90  return
      end
      subroutine cualde(idim,t,n,c,nc,k1,u,d,nd,ier)
c  subroutine cualde evaluates at the point u all the derivatives
c                     (l)
c     d(idim*l+j) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim
c  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim
c  given in its b-spline representation.
c
c  calling sequence:
c     call cualde(idim,t,n,c,nc,k1,u,d,nd,ier)
c
c  input parameters:
c    idim : integer, giving the dimension of the spline curve.
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(u).
c    c    : array,length nc, which contains the b-spline coefficients.
c    nc   : integer, giving the total number of coefficients of s(u).
c    k1   : integer, giving the order of s(u) (order=degree+1).
c    u    : real, which contains the point where the derivatives must
c           be evaluated.
c    nd   : integer, giving the dimension of the array d. nd >= k1*idim
c
c  output parameters:
c    d    : array,length nd,giving the different curve derivatives.
c           d(idim*l+j) will contain the j-th coordinate of the l-th
c           derivative of the curve at the point u.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    nd >= k1*idim
c    t(k1) <= u <= t(n-k1+1)
c
c  further comments:
c    if u coincides with a knot, right derivatives are computed
c    ( left derivatives if u = t(n-k1+1) ).
c
c  other subroutines required: fpader.
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer idim,n,nc,k1,nd,ier
      real u
c  ..array arguments..
      real t(n),c(nc),d(nd)
c  ..local scalars..
      integer i,j,kk,l,m,nk1
c  ..local array..
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nd.lt.(k1*idim)) go to 500
      nk1 = n-k1
      if(u.lt.t(k1) .or. u.gt.t(nk1+1)) go to 500
c  search for knot interval t(l) <= u < t(l+1)
      l = k1
 100  if(u.lt.t(l+1) .or. l.eq.nk1) go to 200
      l = l+1
      go to 100
 200  if(t(l).ge.t(l+1)) go to 500
      ier = 0
c  calculate the derivatives.
      j = 1
      do 400 i=1,idim
        call fpader(t,n,c(j),k1,u,l,h)
        m = i
        do 300 kk=1,k1
          d(m) = h(kk)
          m = m+idim
 300    continue
        j = j+n
 400  continue
 500  return
      end
      subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
c  subroutine curev evaluates in a number of points u(i),i=1,2,...,m
c  a spline curve s(u) of degree k and dimension idim, given in its
c  b-spline representation.
c
c  calling sequence:
c     call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
c
c  input parameters:
c    idim : integer, giving the dimension of the spline curve.
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(u).
c    c    : array,length nc, which contains the b-spline coefficients.
c    nc   : integer, giving the total number of coefficients of s(u).
c    k    : integer, giving the degree of s(u).
c    u    : array,length m, which contains the points where s(u) must
c           be evaluated.
c    m    : integer, giving the number of points where s(u) must be
c           evaluated.
c    mx   : integer, giving the dimension of the array x. mx >= m*idim
c
c  output parameters:
c    x    : array,length mx,giving the value of s(u) at the different
c           points. x(idim*(i-1)+j) will contain the j-th coordinate
c           of the i-th point on the curve.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    mx >= m*idim
c    t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer idim,n,nc,k,m,mx,ier
c  ..array arguments..
      real t(n),c(nc),u(m),x(mx)
c  ..local scalars..
      integer i,j,jj,j1,k1,l,ll,l1,mm,nk1
      real arg,sp,tb,te
c  ..local array..
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mx.lt.(m*idim)) go to 100
      ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      mm = 0
      do 80 i=1,m
c  fetch a new u-value arg.
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(u) at u=arg.
        ll = l-k1
        do 70 j1=1,idim
          jj = ll
          sp = 0.
          do 60 j=1,k1
            jj = jj+1
            sp = sp+c(jj)*h(j)
  60      continue
          mm = mm+1
          x(mm) = sp
          ll = ll+n
  70    continue
  80  continue
 100  return
      end
      subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
c  approximation of degree k on the interval xb <= x <= xe.
c  if iopt=-1 curfit calculates the weighted least-squares spline
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
c           k+1. if iopt=1 the routine will continue with the knots
c           found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   x     : real array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. unchanged on exit.
c   y     : real array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m. unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   xb,xe : real values. on entry xb and xe must specify the boundaries
c           of the approximation interval. xb<=x(1), xe>=x(m).
c           unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is  nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier =10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : real. unless ier=10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program.lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares polynomial of degree k. in this extreme case fp
c             gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the weighted least-squares
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) such that
c                         t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if curfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   curfit once more with the selected value for s but now with iopt=0.
c   indeed, curfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : an algorithm for smoothing, differentiation and integ-
c                ration of experimental data using spline functions,
c                j.comp.appl.maths 1 (1975) 165-184.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : an improved algorithm for curve fitting with spline
c                functions, report tw54, dept. computer science,k.u.
c                leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real xb,xe,s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real tol
      integer i,ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,maxit,nmin
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(7+3*k)
      if(lwrk.lt.lwest) go to 50
      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
      do 10 i=2,m
         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 50
      j = n
      do 20 i=1,k1
         t(i) = xb
         t(j) = xe
         j = j-1
  20  continue
      call fpchec(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia = iz+nest
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
  50  return
      end
      real function dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
c  function dblint calculates the double integral
c         / xe  / ye
c        |     |      s(x,y) dx dy
c    xb /  yb /
c  with s(x,y) a bivariate spline of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   xb,xe : real values, containing the boundaries of the integration
c   yb,ye   domain. s(x,y) is considered to be identically zero out-
c           side the rectangle (tx(kx+1),tx(nx-kx))*(ty(ky+1),ty(ny-ky))
c
c  output parameters:
c   aint  : real , containing the double integral of s(x,y).
c   wrk   : real array of dimension at least (nx+ny-kx-ky-2).
c           used as working space.
c           on exit, wrk(i) will contain the integral
c                / xe
c               | ni,kx+1(x) dx , i=1,2,...,nx-kx-1
c           xb /
c           with ni,kx+1(x) the normalized b-spline defined on
c           the knots tx(i),...,tx(i+kx+1)
c           wrk(j+nx-kx-1) will contain the integral
c                / ye
c               | nj,ky+1(y) dy , j=1,2,...,ny-ky-1
c           yb /
c           with nj,ky+1(y) the normalized b-spline defined on
c           the knots ty(j),...,ty(j+ky+1)
c
c  other subroutines required: fpintb
c
c  references :
c    gaffney p.w. : the calculation of indefinite integrals of b-splines
c                   j. inst. maths applics 17 (1976) 37-41.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nx,ny,kx,ky
      real xb,xe,yb,ye
c  ..array arguments..
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),wrk(nx+ny-kx-ky-2)
c  ..local scalars..
      integer i,j,l,m,nkx1,nky1
      real res
c  ..
      nkx1 = nx-kx-1
      nky1 = ny-ky-1
c  we calculate the integrals of the normalized b-splines ni,kx+1(x)
      call fpintb(tx,nx,wrk,nkx1,xb,xe)
c  we calculate the integrals of the normalized b-splines nj,ky+1(y)
      call fpintb(ty,ny,wrk(nkx1+1),nky1,yb,ye)
c  calculate the integral of s(x,y)
      dblint = 0.
      do 200 i=1,nkx1
        res = wrk(i)
        if(res.eq.0.) go to 200
        m = (i-1)*nky1
        l = nkx1
        do 100 j=1,nky1
          m = m+1
          l = l+1
          dblint = dblint+res*wrk(l)*c(m)
 100    continue
 200  continue
      return
      end
      real function evapol(tu,nu,tv,nv,c,rad,x,y)
c  function program evacir evaluates the function f(x,y) = s(u,v),
c  defined through the transformation
c      x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v)
c  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given
c  in its standard b-spline representation.
c
c  calling sequence:
c     f = evapol(tu,nu,tv,nv,c,rad,x,y)
c
c  input parameters:
c   tu    : real array, length nu, which contains the position of the
c           knots in the u-direction.
c   nu    : integer, giving the total number of knots in the u-direction
c   tv    : real array, length nv, which contains the position of the
c           knots in the v-direction.
c   nv    : integer, giving the total number of knots in the v-direction
c   c     : real array, length (nu-4)*(nv-4), which contains the
c           b-spline coefficients.
c   rad   : real function subprogram, defining the boundary of the
c           approximation domain. must be declared external in the
c           calling (sub)-program
c   x,y   : real values.
c           before entry x and y must be set to the co-ordinates of
c           the point where f(x,y) must be evaluated.
c
c  output parameter:
c   f     : real
c           on exit f contains the value of f(x,y)
c
c  other subroutines required:
c    bispev,fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nu,nv
      real x,y
c  ..array arguments..
      real tu(nu),tv(nv),c((nu-4)*(nv-4))
c  ..user specified function
      real rad
c  ..local scalars..
      integer ier
      real u,v,r,f,one,dist
c  ..local arrays
      real wrk(8)
      integer iwrk(2)
c  ..function references
      real atan2,sqrt
c  ..
c  calculate the (u,v)-coordinates of the given point.
      one = 1
      u = 0.
      v = 0.
      dist = x**2+y**2
      if(dist.le.0.) go to 10
      v = atan2(y,x)
      r = rad(v)
      if(r.le.0.) go to 10
      u = sqrt(dist)/r
      if(u.gt.one) u = one
c  evaluate s(u,v)
  10  call bispev(tu,nu,tv,nv,c,3,3,u,1,v,1,f,wrk,8,iwrk,2,ier)
      evapol = f
      return
      end

      subroutine fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
c  subroutine fourco calculates the integrals
c                    /t(n-3)
c    ress(i) =      !        s(x)*sin(alfa(i)*x) dx    and
c              t(4)/
c                    /t(n-3)
c    resc(i) =      !        s(x)*cos(alfa(i)*x) dx, i=1,...,m,
c              t(4)/
c  where s(x) denotes a cubic spline which is given in its
c  b-spline representation.
c
c  calling sequence:
c     call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
c
c  input parameters:
c    t    : real array,length n, containing the knots of s(x).
c    n    : integer, containing the total number of knots. n>=10.
c    c    : real array,length n, containing the b-spline coefficients.
c    alfa : real array,length m, containing the parameters alfa(i).
c    m    : integer, specifying the number of integrals to be computed.
c    wrk1 : real array,length n. used as working space
c    wrk2 : real array,length n. used as working space
c
c  output parameters:
c    ress : real array,length m, containing the integrals ress(i).
c    resc : real array,length m, containing the integrals resc(i).
c    ier  : error flag:
c      ier=0 : normal return.
c      ier=10: invalid input data (see restrictions).
c
c  restrictions:
c    n >= 10
c    t(4) < t(5) < ... < t(n-4) < t(n-3).
c    t(1) <= t(2) <= t(3) <= t(4).
c    t(n-3) <= t(n-2) <= t(n-1) <= t(n).
c
c  other subroutines required: fpbfou,fpcsin
c
c  references :
c    dierckx p. : calculation of fouriercoefficients of discrete
c                 functions using cubic splines. j. computational
c                 and applied mathematics 3 (1977) 207-209.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,m,ier
c  ..array arguments..
      real t(n),c(n),wrk1(n),wrk2(n),alfa(m),ress(m),resc(m)
c  ..local scalars..
      integer i,j,n4
      real rs,rc
c  ..
      n4 = n-4
c  before starting computations a data check is made. in the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(n.lt.10) go to 50
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 50
        if(t(j).lt.t(j-1)) go to 50
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 50
  20  continue
      ier = 0
c  main loop for the different alfa(i).
      do 40 i=1,m
c  calculate the integrals
c    wrk1(j) = integral(nj,4(x)*sin(alfa*x))    and
c    wrk2(j) = integral(nj,4(x)*cos(alfa*x)),  j=1,2,...,n-4,
c  where nj,4(x) denotes the normalised cubic b-spline defined on the
c  knots t(j),t(j+1),...,t(j+4).
         call fpbfou(t,n,alfa(i),wrk1,wrk2)
c  calculate the integrals ress(i) and resc(i).
         rs = 0.
         rc = 0.
         do 30 j=1,n4
            rs = rs+c(j)*wrk1(j)
            rc = rc+c(j)*wrk2(j)
  30     continue
         ress(i) = rs
         resc(i) = rc
  40  continue
  50  return
      end
      subroutine fpader(t,n,c,k1,x,l,d)
c  subroutine fpader calculates the derivatives
c             (j-1)
c     d(j) = s     (x) , j=1,2,...,k1
c  of a spline of order k1 at the point t(l)<=x<t(l+1), using the
c  stable recurrence scheme of de boor
c  ..
c  ..scalar arguments..
      real x
      integer n,k1,l
c  ..array arguments..
      real t(n),c(n),d(k1)
c  ..local scalars..
      integer i,ik,j,jj,j1,j2,ki,kj,li,lj,lk
      real ak,fac,one
c  ..local array..
      real h(6)
c  ..
      one = 0.1e+01
      lk = l-k1
      do 100 i=1,k1
        ik = i+lk
        h(i) = c(ik)
 100  continue
      kj = k1
      fac = one
      do 700 j=1,k1
        ki = kj
        j1 = j+1
        if(j.eq.1) go to 300
        i = k1
        do 200 jj=j,k1
          li = i+lk
          lj = li+kj
          h(i) = (h(i)-h(i-1))/(t(lj)-t(li))
          i = i-1
 200    continue
 300    do 400 i=j,k1
          d(i) = h(i)
 400    continue
        if(j.eq.k1) go to 600
        do 500 jj=j1,k1
          ki = ki-1
          i = k1
          do 500 j2=jj,k1
            li = i+lk
            lj = li+ki
            d(i) = ((x-t(li))*d(i)+(t(lj)-x)*d(i-1))/(t(lj)-t(li))
            i = i-1
 500    continue
 600    d(j) = d(k1)*fac
        ak = k1-j
        fac = fac*ak
        kj = kj-1
 700  continue
      return
      end
      subroutine fpadno(maxtr,up,left,right,info,count,merk,jbind,
     * n1,ier)
c  subroutine fpadno adds a branch of length n1 to the triply linked
c  tree,the information of which is kept in the arrays up,left,right
c  and info. the information field of the nodes of this new branch is
c  given in the array jbind. in linking the new branch fpadno takes
c  account of the property of the tree that
c    info(k) < info(right(k)) ; info(k) < info(left(k))
c  if necessary the subroutine calls subroutine fpfrno to collect the
c  free nodes of the tree. if no computer words are available at that
c  moment, the error parameter ier is set to 1.
c  ..
c  ..scalar arguments..
      integer maxtr,count,merk,n1,ier
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),jbind(n1)
c  ..local scalars..
      integer k,niveau,point
      logical bool
c  ..subroutine references..
c    fpfrno
c  ..
      point = 1
      niveau = 1
  10  k = left(point)
      bool = .true.
  20  if(k.eq.0) go to 50
      if(info(k)-jbind(niveau)) 30,40,50
  30  point = k
      k = right(point)
      bool = .false.
      go to 20
  40  point = k
      niveau = niveau+1
      go to 10
  50  if(niveau.gt.n1) go to 90
      count = count+1
      if(count.le.maxtr) go to 60
      call fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier)
      if(ier.ne.0) go to 100
  60  info(count) = jbind(niveau)
      left(count) = 0
      right(count) = k
      if(bool) go to 70
      bool = .true.
      right(point) = count
      up(count) = up(point)
      go to 80
  70  up(count) = point
      left(point) = count
  80  point = count
      niveau = niveau+1
      k = 0
      go to 50
  90  ier = 0
 100  return
      end
      subroutine fpadpo(idim,t,n,c,nc,k,cp,np,cc,t1,t2)
c  given a idim-dimensional spline curve of degree k, in its b-spline
c  representation ( knots t(j),j=1,...,n , b-spline coefficients c(j),
c  j=1,...,nc) and given also a polynomial curve in its b-spline
c  representation ( coefficients cp(j), j=1,...,np), subroutine fpadpo
c  calculates the b-spline representation (coefficients c(j),j=1,...,nc)
c  of the sum of the two curves.
c
c  other subroutine required : fpinst
c
c  ..
c  ..scalar arguments..
      integer idim,k,n,nc,np
c  ..array arguments..
      real t(n),c(nc),cp(np),cc(nc),t1(n),t2(n)
c  ..local scalars..
      integer i,ii,j,jj,k1,l,l1,n1,n2,nk1,nk2
c  ..
      k1 = k+1
      nk1 = n-k1
c  initialization
      j = 1
      l = 1
      do 20 jj=1,idim
        l1 = j
        do 10 ii=1,k1
          cc(l1) = cp(l)
          l1 = l1+1
          l = l+1
  10    continue
        j = j+n
        l = l+k1
  20  continue
      if(nk1.eq.k1) go to 70
      n1 = k1*2
      j = n
      l = n1
      do 30 i=1,k1
        t1(i) = t(i)
        t1(l) = t(j)
        l = l-1
        j = j-1
  30  continue
c  find the b-spline representation of the given polynomial curve
c  according to the given set of knots.
      nk2 = nk1-1
      do 60 l=k1,nk2
        l1 = l+1
        j = 1
        do 40 i=1,idim
          call fpinst(0,t1,n1,cc(j),k,t(l1),l,t2,n2,cc(j),n)
          j = j+n
  40    continue
        do 50 i=1,n2
          t1(i) = t2(i)
  50    continue
        n1 = n2
  60  continue
c  find the b-spline representation of the resulting curve.
  70  j = 1
      do 90 jj=1,idim
        l = j
        do 80 i=1,nk1
          c(l) = cc(l)+c(l)
          l = l+1
  80    continue
        j = j+n
  90  continue
      return
      end
      subroutine fpback(a,z,n,k,c,nest)
c  subroutine fpback calculates the solution of the system of
c  equations a*c = z with a a n x n upper triangular matrix
c  of bandwidth k.
c  ..
c  ..scalar arguments..
      integer n,k,nest
c  ..array arguments..
      real a(nest,k),z(n),c(n)
c  ..local scalars..
      real store
      integer i,i1,j,k1,l,m
c  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
      subroutine fpbacp(a,b,z,n,k,c,k1,nest)
c  subroutine fpbacp calculates the solution of the system of equations
c  g * c = z  with g  a n x n upper triangular matrix of the form
c            ! a '   !
c        g = !   ' b !
c            ! 0 '   !
c  with b a n x k matrix and a a (n-k) x (n-k) upper triangular
c  matrix of bandwidth k1.
c  ..
c  ..scalar arguments..
      integer n,k,k1,nest
c  ..array arguments..
      real a(nest,k1),b(nest,k),z(n),c(n)
c  ..local scalars..
      integer i,i1,j,l,l0,l1,n2
      real store
c  ..
      n2 = n-k
      l = n
      do 30 i=1,k
        store = z(l)
        j = k+2-i
        if(i.eq.1) go to 20
        l0 = l
        do 10 l1=j,k
          l0 = l0+1
          store = store-c(l0)*b(l,l1)
  10    continue
  20    c(l) = store/b(l,j-1)
        l = l-1
        if(l.eq.0) go to 80
  30  continue
      do 50 i=1,n2
        store = z(i)
        l = n2
        do 40 j=1,k
          l = l+1
          store = store-c(l)*b(i,j)
  40    continue
        c(i) = store
  50  continue
      i = n2
      c(i) = c(i)/a(i,1)
      if(i.eq.1) go to 80
      do 70 j=2,n2
        i = i-1
        store = c(i)
        i1 = k
        if(j.le.k) i1=j-1
        l = i
        do 60 l0=1,i1
          l = l+1
          store = store-c(l)*a(i,l0+1)
  60    continue
        c(i) = store/a(i,1)
  70  continue
  80  return
      end
      subroutine fpbfou(t,n,par,ress,resc)
c  subroutine fpbfou calculates the integrals
c                    /t(n-3)
c    ress(j) =      !        nj,4(x)*sin(par*x) dx    and
c              t(4)/
c                    /t(n-3)
c    resc(j) =      !        nj,4(x)*cos(par*x) dx ,  j=1,2,...n-4
c              t(4)/
c  where nj,4(x) denotes the cubic b-spline defined on the knots
c  t(j),t(j+1),...,t(j+4).
c
c  calling sequence:
c     call fpbfou(t,n,par,ress,resc)
c
c  input parameters:
c    t    : real array,length n, containing the knots.
c    n    : integer, containing the number of knots.
c    par  : real, containing the value of the parameter par.
c
c  output parameters:
c    ress  : real array,length n, containing the integrals ress(j).
c    resc  : real array,length n, containing the integrals resc(j).
c
c  restrictions:
c    n >= 10, t(4) < t(5) < ... < t(n-4) < t(n-3).
c  ..
c  ..scalar arguments..
      integer n
      real par
c  ..array arguments..
      real t(n),ress(n),resc(n)
c  ..local scalars..
      integer i,ic,ipj,is,j,jj,jp1,jp4,k,li,lj,ll,nmj,nm3,nm7
      real ak,beta,con1,con2,c1,c2,delta,eps,fac,f1,f2,f3,one,quart,
     * sign,six,s1,s2,term
c  ..local arrays..
      real co(5),si(5),hs(5),hc(5),rs(3),rc(3)
c  ..function references..
      real cos,sin,abs
c  ..
c  initialization.
      one = 0.1e+01
      six = 0.6e+01
      eps = 0.1e-07
      quart = 0.25e0
      con1 = 0.5e-01
      con2 = 0.12e+03
      nm3 = n-3
      nm7 = n-7
      if(par.ne.0.) term = six/par
      beta = par*t(4)
      co(1) = cos(beta)
      si(1) = sin(beta)
c  calculate the integrals ress(j) and resc(j), j=1,2,3 by setting up
c  a divided difference table.
      do 30 j=1,3
        jp1 = j+1
        jp4 = j+4
        beta = par*t(jp4)
        co(jp1) = cos(beta)
        si(jp1) = sin(beta)
        call fpcsin(t(4),t(jp4),par,si(1),co(1),si(jp1),co(jp1),
     *  rs(j),rc(j))
        i = 5-j
        hs(i) = 0.
        hc(i) = 0.
        do 10 jj=1,j
          ipj = i+jj
          hs(ipj) = rs(jj)
          hc(ipj) = rc(jj)
  10    continue
        do 20 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = jp4
          do 20 ll=i,4
            lj = li-jj
            fac = t(li)-t(lj)
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  20    continue
        ress(j) = hs(5)-hs(4)
        resc(j) = hc(5)-hc(4)
  30  continue
      if(nm7.lt.4) go to 160
c  calculate the integrals ress(j) and resc(j),j=4,5,...,n-7.
      do 150 j=4,nm7
        jp4 = j+4
        beta = par*t(jp4)
        co(5) = cos(beta)
        si(5) = sin(beta)
        delta = t(jp4)-t(j)
c  the way of computing ress(j) and resc(j) depends on the value of
c  beta = par*(t(j+4)-t(j)).
        beta = delta*par
        if(abs(beta).le.one) go to 60
c  if !beta! > 1 the integrals are calculated by setting up a divided
c  difference table.
        do 40 k=1,5
          hs(k) = si(k)
          hc(k) = co(k)
  40    continue
        do 50 jj=1,3
          k = 5
          li = jp4
          do 50 ll=jj,4
            lj = li-jj
            fac = par*(t(li)-t(lj))
            hs(k) = (hs(k)-hs(k-1))/fac
            hc(k) = (hc(k)-hc(k-1))/fac
            k = k-1
            li = li-1
  50    continue
        s2 = (hs(5)-hs(4))*term
        c2 = (hc(5)-hc(4))*term
        go to 130
c  if !beta! <= 1 the integrals are calculated by evaluating a series
c  expansion.
  60    f3 = 0.
        do 70 i=1,4
          ipj = i+j
          hs(i) = par*(t(ipj)-t(j))
          hc(i) = hs(i)
          f3 = f3+hs(i)
  70    continue
        f3 = f3*con1
        c1 = quart
        s1 = f3
        if(abs(f3).le.eps) go to 120
        sign = one
        fac = con2
        k = 5
        is = 0
        do 110 ic=1,20
          k = k+1
          ak = k
          fac = fac*ak
          f1 = 0.
          f3 = 0.
          do 80 i=1,4
            f1 = f1+hc(i)
            f2 = f1*hs(i)
            hc(i) = f2
            f3 = f3+f2
  80      continue
          f3 = f3*six/fac
          if(is.eq.0) go to 90
          is = 0
          s1 = s1+f3*sign
          go to 100
  90      sign = -sign
          is = 1
          c1 = c1+f3*sign
 100      if(abs(f3).le.eps) go to 120
 110    continue
 120    s2 = delta*(co(1)*s1+si(1)*c1)
        c2 = delta*(co(1)*c1-si(1)*s1)
 130    ress(j) = s2
        resc(j) = c2
        do 140 i=1,4
          co(i) = co(i+1)
          si(i) = si(i+1)
 140    continue
 150  continue
c  calculate the integrals ress(j) and resc(j),j=n-6,n-5,n-4 by setting
c  up a divided difference table.
 160  do 190 j=1,3
        nmj = nm3-j
        i = 5-j
        call fpcsin(t(nm3),t(nmj),par,si(4),co(4),si(i-1),co(i-1),
     *  rs(j),rc(j))
        hs(i) = 0.
        hc(i) = 0.
        do 170 jj=1,j
          ipj = i+jj
          hc(ipj) = rc(jj)
          hs(ipj) = rs(jj)
 170    continue
        do 180 jj=1,3
          if(i.lt.jj) i = jj
          k = 5
          li = nmj
          do 180 ll=i,4
            lj = li+jj
            fac = t(lj)-t(li)
            hs(k) = (hs(k-1)-hs(k))/fac
            hc(k) = (hc(k-1)-hc(k))/fac
            k = k-1
            li = li+1
 180    continue
        ress(nmj) = hs(4)-hs(5)
        resc(nmj) = hc(4)-hc(5)
 190  continue
      return
      end
      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
c  ..array arguments..
      integer lx(mx),ly(my)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wx(mx,kx+1),wy(my,ky+1)
c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real arg,sp,tb,te
c  ..local arrays..
      real h(6)
c  ..subroutine references..
c    fpbspl
c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end

      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      real x
      integer n,k,l
c  ..array arguments..
      real t(n),h(6)
c  ..local scalars..
      real f,one
      integer i,j,li,lj
c  ..local arrays..
      real hh(5)
c  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
      subroutine fpchec(x,m,t,n,k,ier)
c  subroutine fpchec verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
c  and the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      real x(m),t(n)
c  ..local scalars..
      integer i,j,k1,k2,l,nk1,nk2,nk3
      real tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
      subroutine fpched(x,m,t,n,k,ib,ie,ier)
c  subroutine fpched verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a spline of degree k,with ib derative constraints
c  at x(1) and ie constraints at x(m), in relation to the number and
c  the position of the data points x(i),i=1,2,...,m. if all of the
c  following conditions are fulfilled, the error parameter ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m + max(0,ib-1) + max(0,ie-1)
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=1+ib1,2+ib1,...,n-k-1-ie1
c               with ib1 = max(0,ib-1), ie1 = max(0,ie-1)
c  ..
c  ..scalar arguments..
      integer m,n,k,ib,ie,ier
c  ..array arguments..
      real x(m),t(n)
c  ..local scalars..
      integer i,ib1,ie1,j,jj,k1,k2,l,nk1,nk2,nk3
      real tj,tl
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ib1 = ib-1
      if(ib1.lt.0) ib1 = 0
      ie1 = ie-1
      if(ie1.lt.0) ie1 = 0
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.(m+ib1+ie1)) go to 80
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
c  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      jj = 2+ib1
      l = jj+k
      nk3 = nk1-1-ie1
      if(nk3.lt.jj) go to 70
      do 60 j=jj,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
      subroutine fpchep(x,m,t,n,k,ier)
c  subroutine fpchep verifies the number and the position of the knots
c  t(j),j=1,2,...,n of a periodic spline of degree k, in relation to
c  the number and the position of the data points x(i),i=1,2,...,m.
c  if all of the following conditions are fulfilled, ier is set
c  to zero. if one of the conditions is violated ier is set to ten.
c      1) k+1 <= n-k-1 <= m+k-1
c      2) t(1) <= t(2) <= ... <= t(k+1)
c         t(n-k) <= t(n-k+1) <= ... <= t(n)
c      3) t(k+1) < t(k+2) < ... < t(n-k)
c      4) t(k+1) <= x(i) <= t(n-k)
c      5) the conditions specified by schoenberg and whitney must hold
c         for at least one subset of data points, i.e. there must be a
c         subset of data points y(j) such that
c             t(j) < y(j) < t(j+k+1), j=k+1,...,n-k-1
c  ..
c  ..scalar arguments..
      integer m,n,k,ier
c  ..array arguments..
      real x(m),t(n)
c  ..local scalars..
      integer i,i1,i2,j,j1,k1,k2,l,l1,l2,mm,m1,nk1,nk2
      real per,tj,tl,xi
c  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      m1 = m-1
      ier = 10
c  check condition no 1
      if(nk1.lt.k1 .or. n.gt.m+2*k) go to 130
c  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 130
        if(t(j).lt.t(j-1)) go to 130
        j = j-1
  20  continue
c  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 130
  30  continue
c  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 130
c  check condition no 5
      l1 = k1
      l2 = 1
      do 50 l=1,m
         xi = x(l)
  40     if(xi.lt.t(l1+1) .or. l.eq.nk1) go to 50
         l1 = l1+1
         l2 = l2+1
         if(l2.gt.k1) go to 60
         go to 40
  50  continue
      l = m
  60  per = t(nk2)-t(k1)
      do 120 i1=2,l
         i = i1-1
         mm = i+m1
         do 110 j=k1,nk1
            tj = t(j)
            j1 = j+k1
            tl = t(j1)
  70        i = i+1
            if(i.gt.mm) go to 120
            i2 = i-m1
            if(i2) 80,80,90
  80        xi = x(i)
            go to 100
  90        xi = x(i2)+per
 100        if(xi.le.tj) go to 70
            if(xi.ge.tl) go to 120
 110     continue
         ier = 0
         go to 130
 120  continue
 130  return
      end
      subroutine fpclos(iopt,idim,m,u,mx,x,w,k,s,nest,tol,maxit,k1,k2,
     * n,t,nc,c,fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real s,tol,fp
      integer iopt,idim,m,mx,k,nest,maxit,k1,k2,n,nc,ier
c  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),z(nc),a1(nest,k1),
     * a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,cos,d1,fac,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv,
     * p1,p2,p3,sin,store,term,ui,wi,rn,one,con1,con4,con9,half
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jj,jk,jper,j1,j2,kk,
     * kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1,
     * nrint,n10,n11,n7,n8
c  ..local arrays..
      real h(6),h1(7),h2(6),xi(10)
c  ..function references..
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares closed curve      c
c  sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    c
c  if iopt=-1 sinf(u) is the requested curve                           c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares curve until finally fp<=s.         c
c  the initial choice of knots depends on the value of s and iopt.     c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+2*k.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial curve of   c
c      degree k; n = nmin = 2*k+2. since s(u) must be periodic we      c
c      find that s(u) reduces to a fixed point.                        c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial curve.         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
c  determine the length of the period of the splines.
      per = u(m)-u(1)
      if(iopt.lt.0) go to 50
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s.gt.0. .or. nmax.eq.nmin) go to 30
c  if s=0, s(u) is an interpolating curve.
      n = nmax
c  test whether the required storage space exceeds the available one.
      if(n.gt.nest) go to 620
c  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 .eq.k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = u(i)
  10  continue
      if(s.gt.0.) go to 50
      kk = k-1
      kk1 = k
      if(kk.gt.0) go to 50
      t(1) = t(m)-per
      t(2) = u(1)
      t(m+1) = u(m)
      t(m+2) = t(3)+per
      jj = 0
      do 15 i=1,m1
        j = i
        do 12 j1=1,idim
          jj = jj+1
          c(j) = x(jj)
          j = j+n
  12    continue
  15  continue
      jj = 1
      j = m
      do 17 j1=1,idim
        c(j) = c(jj)
        j = j+n
        jj = jj+n
  17  continue
      fp = 0.
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
         j = i+k
         t(j) = (u(i)+u(i-1))*half
  25  continue
      go to 50
c  if s > 0 our initial choice depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial curve. (i.e. a constant point).
c  if iopt=1 and fp0>s we start computing the least-squares closed
c  curve according the set of knots found at the last call of the
c  routine.
  30  if(iopt.eq.0) go to 35
      if(n.eq.nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 50
c  the case that s(u) is a fixed point is treated separetely.
c  fp0 denotes the corresponding sum of squared residuals.
  35  fp0 = 0.
      d1 = 0.
      do 37 j=1,idim
        z(j) = 0.
  37  continue
      jj = 0
      do 45 it=1,m1
        wi = w(it)
        call fpgivs(wi,d1,cos,sin)
        do 40 j=1,idim
          jj = jj+1
          fac = wi*x(jj)
          call fprota(cos,sin,fac,z(j))
          fp0 = fp0+fac**2
  40    continue
  45  continue
      do 47 j=1,idim
        z(j) = z(j)/d1
  47  continue
c  test whether that fixed point is a solution of our problem.
      fpms = fp0-s
      if(fpms.lt.acc .or. nmax.eq.nmin) go to 640
      fpold = fp0
c  test whether the required storage space exceeds the available one.
      if(n.ge.nest) go to 620
c  start computing the least-squares closed curve with one
c  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = u(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
c  main loop for the different sets of knots. m is a save upper
c  bound for the number of trials.
  50  do 340 iter=1,m
c  find nrint, the number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(u). if we take
c      t(k+1) = u(1), t(n-k) = u(m)
c      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
c      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
c  then s(u) will be a smooth closed curve if the b-spline
c  coefficients satisfy the following conditions
c      c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
c  with n7=n-2*k-1.
        t(k1) = u(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = u(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
c  compute the b-spline coefficients of the least-squares closed curve
c  sinf(u). the observation matrix a is built up row by row while
c  taking into account condition (**) and is reduced to triangular
c  form by givens transformations .
c  at the same time fp=f(p=inf) is computed.
c  the n7 x n7 triangularised upper matrix a has the form
c            ! a1 '    !
c        a = !    ' a2 !
c            ! 0  '    !
c  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
c  matrix of bandwith k+1 ( n10 = n7-k).
c  initialization.
        do 65 i=1,nc
          z(i) = 0.
  65    continue
        do 70 i=1,nk1
          do 70 j=1,kk1
            a1(i,j) = 0.
  70    continue
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = 0.
        l = k1
        jj = 0
        do 290 it=1,m1
c  fetch the current data point u(it),x(it)
          ui = u(it)
          wi = w(it)
          do 75 j=1,idim
            jj = jj+1
            xi(j) = x(jj)*wi
  75      continue
c  search for knot interval t(l) <= ui < t(l+1).
  80      if(ui.lt.t(l+1)) go to 85
          l = l+1
          go to 80
c  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  85      call fpbspl(t,n,k,ui,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
c  test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
          if(l5.lt.n10) go to 285
          if(jper.ne.0) go to 160
c  initialize the matrix a2.
          do 95 i=1,n7
          do 95 j=1,kk
              a2(i,j) = 0.
  95      continue
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik.le.0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
c  if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
c  we take account of condition (**) for setting up the new row
c  of the observation matrix a. this row is stored in the arrays h1
c  (the part with respect to a1) and h2 (the part with
c  respect to a2).
 160      do 170 i=1,kk
            h1(i) = 0.
            h2(i) = 0.
 170      continue
          h1(kk1) = 0.
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1.le.0) go to 200
            if(l1.le.n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
c  rotate the new row of the observation matrix into triangle
c  by givens transformations.
          if(n10.le.0) go to 250
c  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv.ne.0.) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.
            go to 240
c  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
c  transformation to the right hand side.
            j1 = j
            do 217 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 217        continue
c  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j.eq.n10) go to 250
            i2 = min0(n10-j,kk)
c  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = 0.
 240      continue
c  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij.le.0) go to 270
            piv = h2(j)
            if(piv.eq.0.) go to 270
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
c  transformations to right hand side.
            j1 = ij
            do 255 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 255        continue
            if(j.eq.kk) go to 280
            j1 = j+1
c  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 280      do 282 j2=1,idim
            fp = fp+xi(j2)**2
 282      continue
          go to 290
c  rotation of the new row of the observation matrix into
c  triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
c  at ui.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 140
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
c  transformations to right hand side.
            j1 = j
            do 125 j2=1,idim
              call fprota(cos,sin,xi(j2),z(j1))
              j1 = j1+n
 125        continue
            if(i.eq.kk1) go to 150
            i2 = 1
            i3 = i+1
c  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 150      do 155 j2=1,idim
            fp = fp+xi(j2)**2
 155      continue
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients .
        j1 = 1
        do 292 j2=1,idim
           call fpbacp(a1,a2,z(j1),n7,kk,c(j1),kk1,nest)
           j1 = j1+n
 292    continue
c  calculate from condition (**) the remaining coefficients.
        do 297 i=1,k
          j1 = i
          do 295 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 295      continue
 297    continue
        if(iopt.lt.0) go to 660
c  test whether the approximation sinf(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 350
c  if n=nmax, sinf(u) is an interpolating curve.
        if(n.eq.nmax) go to 630
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of the
c  storage capacity limitation.
        if(n.eq.nest) go to 620
c  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
c  compute the sum of squared residuals for each knot interval
c  t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k1
        jj = 0
        do 320 it=1,m1
          if(u(it).lt.t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.
          l0 = l-k2
          do 310 j2=1,idim
            fac = 0.
            j1 = l0
            do 305 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 305        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 310      continue
          fpart = fpart+term
          if(new.eq.0) go to 320
          if(l.gt.k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
c  add a new knot
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation
          if(n.eq.nmax) go to 5
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 340
 330    continue
c  restart the computations with the new set of knots.
 340  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing closed curve sp(u).          c
c  **********************************************************          c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing curve     c
c  sp(u). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(u) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that f(p),c
c  the sum of squared residuals be = s. we already know that the least-c
c  squares polynomial curve corresponds to p=0, and that the least-    c
c  squares periodic spline curve corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l.eq.0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
c  form the matrix g  as the matrix a extended by the rows of matrix b.
c  the rows of matrix b with weight 1/p are rotated into
c  the triangularised observation matrix a.
c  after triangularisation our n7 x n7 matrix g takes the form
c            ! g1 '    !
c        g = !    ' g2 !
c            ! 0  '    !
c  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
c  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
c  store matrix a into g
        do 358 i=1,nc
          c(i) = z(i)
 358    continue
        do 360 i=1,n7
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = 0.
          g2(i,1) = 0.
          do 360 j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
 360    continue
        l = n10
        do 370 j=1,k1
          if(l.le.0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
c  fetch a new row of matrix b and store it in the arrays h1 (the part
c  with respect to g1) and h2 (the part with respect to g2).
          do 380 j=1,idim
            xi(j) = 0.
 380      continue
          do 385 i=1,k1
            h1(i) = 0.
            h2(i) = 0.
 385      continue
          h1(k2) = 0.
          if(it.gt.n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0.eq.n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1.le.0) go to 450
            if(l1.le.n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11.le.0) go to 510
c  rotate this row into triangle by givens transformations
c  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
c  transformation to right hand side.
            j1 = j
            do 475 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 475        continue
c  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j.eq.n11) go to 510
            i2 = min0(n11-j,k1)
c  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.
 500      continue
c  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij.le.0) go to 530
            piv = h2(j)
c  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
c  transformation to the right hand side.
            j1 = ij
            do 515 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 515        continue
            if(j.eq.k1) go to 540
            j1 = j+1
c  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
c  backward substitution to obtain the b-spline coefficients
        j1 = 1
        do 542 j2=1,idim
          call fpbacp(g1,g2,c(j1),n7,k1,c(j1),k2,nest)
          j1 = j1+n
 542    continue
c  calculate from condition (**) the remaining b-spline coefficients.
        do 547 i=1,k
          j1 = i
          do 545 j=1,idim
            j2 = j1+n7
            c(j2) = c(j1)
            j1 = j1+n
 545      continue
 547    continue
c  computation of f(p).
        fp = 0.
        l = k1
        jj = 0
        do 570 it=1,m1
          if(u(it).lt.t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.
          do 565 j2=1,idim
            fac = 0.
            j1 = l0
            do 560 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 560        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 565      continue
          fp = fp+term*w(it)**2
 570    continue
c  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 600
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 580
        if((f2-f3) .gt. acc) go to 575
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2.lt.0.) ich3 = 1
 580    if(ich1.ne.0) go to 590
        if((f1-f2) .gt. acc) go to 585
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 595
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 590    if(f2.ge.f1 .or. f2.le.f3) go to 610
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 595  continue
c  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
c  the point (z(1),z(2),...,z(idim)) is a solution of our problem.
c  a constant function is a spline of degree k with all b-spline
c  coefficients equal to that constant.
      do 650 i=1,k1
        rn = k1-i
        t(i) = u(1)-rn*per
        j = i+k1
        rn = i-1
        t(j) = u(m)+rn*per
 650  continue
      n = nmin
      j1 = 0
      do 658 j=1,idim
        fac = z(j)
        j2 = j1
        do 654 i=1,k1
          j2 = j2+1
          c(j2) = fac
 654    continue
        j1 = j1+n
 658  continue
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
 660  return
      end
      subroutine fpcoco(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,
     * bind,e,wrk,lwrk,iwrk,kwrk,ier)
c  ..scalar arguments..
      real s,sq
      integer iopt,m,nest,maxtr,maxbin,n,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real x(m),y(m),w(m),v(m),t(nest),c(nest),sx(m),e(nest),wrk(lwrk)
      logical bind(nest)
c  ..local scalars..
      integer i,ia,ib,ic,iq,iu,iz,izz,i1,j,k,l,l1,m1,nmax,nr,n4,n6,n8,
     * ji,jib,jjb,jl,jr,ju,mb,nm
      real sql,sqmax,term,tj,xi,half
c  ..subroutine references..
c    fpcosp,fpbspl,fpadno,fpdeno,fpseno,fpfrno
c  ..
c  set constant
      half = 0.5e0
c  determine the maximal admissible number of knots.
      nmax = m+4
c  the initial choice of knots depends on the value of iopt.
c    if iopt=0 the program starts with the minimal number of knots
c    so that can be guarantied that the concavity/convexity constraints
c    will be satisfied.
c    if iopt = 1 the program will continue from the point on where she
c    left at the foregoing call.
      if(iopt.gt.0) go to 80
c  find the minimal number of knots.
c  a knot is located at the data point x(i), i=2,3,...m-1 if
c    1) v(i) ^= 0    and
c    2) v(i)*v(i-1) <= 0  or  v(i)*v(i+1) <= 0.
      m1 = m-1
      n = 4
      do 20 i=2,m1
        if(v(i).eq.0. .or. (v(i)*v(i-1).gt.0. .and.
     *  v(i)*v(i+1).gt.0.)) go to 20
        n = n+1
c  test whether the required storage space exceeds the available one.
        if(n+4.gt.nest) go to 200
        t(n) = x(i)
  20  continue
c  find the position of the knots t(1),...t(4) and t(n-3),...t(n) which
c  are needed for the b-spline representation of s(x).
      do 30 i=1,4
        t(i) = x(1)
        n = n+1
        t(n) = x(m)
  30  continue
c  test whether the minimum number of knots exceeds the maximum number.
      if(n.gt.nmax) go to 210
c  main loop for the different sets of knots.
c  find corresponding values e(j) to the knots t(j+3),j=1,2,...n-6
c    e(j) will take the value -1,1, or 0 according to the requirement
c    that s(x) must be locally convex or concave at t(j+3) or that the
c    sign of s''(x) is unrestricted at that point.
  40  i= 1
      xi = x(1)
      j = 4
      tj = t(4)
      n6 = n-6
      do 70 l=1,n6
  50    if(xi.eq.tj) go to 60
        i = i+1
        xi = x(i)
        go to 50
  60    e(l) = v(i)
        j = j+1
        tj = t(j)
  70  continue
c  we partition the working space
      nm = n+maxbin
      mb = maxbin+1
      ia = 1
      ib = ia+4*n
      ic = ib+nm*maxbin
      iz = ic+n
      izz = iz+n
      iu = izz+n
      iq = iu+maxbin
      ji = 1
      ju = ji+maxtr
      jl = ju+maxtr
      jr = jl+maxtr
      jjb = jr+maxtr
      jib = jjb+mb
c  given the set of knots t(j),j=1,2,...n, find the least-squares cubic
c  spline which satisfies the imposed concavity/convexity constraints.
      call fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,wrk(ia),
     * wrk(ib),wrk(ic),wrk(iz),wrk(izz),wrk(iu),wrk(iq),iwrk(ji),
     * iwrk(ju),iwrk(jl),iwrk(jr),iwrk(jjb),iwrk(jib),ier)
c  if sq <= s or in case of abnormal exit from fpcosp, control is
c  repassed to the driver program.
      if(sq.le.s .or. ier.gt.0) go to 300
c  calculate for each knot interval t(l-1) <= xi <= t(l) the
c  sum((wi*(yi-s(xi)))**2).
c  find the interval t(k-1) <= x <= t(k) for which this sum is maximal
c  on the condition that this interval contains at least one interior
c  data point x(nr) and that s(x) is not given there by a straight line.
  80  sqmax = 0.
      sql = 0.
      l = 5
      nr = 0
      i1 = 1
      n4 = n-4
      do 110 i=1,m
        term = (w(i)*(sx(i)-y(i)))**2
        if(x(i).lt.t(l) .or. l.gt.n4) go to 100
        term = term*half
        sql = sql+term
        if(i-i1.le.1 .or. (bind(l-4).and.bind(l-3))) go to 90
        if(sql.le.sqmax) go to 90
        k = l
        sqmax = sql
        nr = i1+(i-i1)/2
  90    l = l+1
        i1 = i
        sql = 0.
 100    sql = sql+term
 110  continue
      if(m-i1.le.1 .or. (bind(l-4).and.bind(l-3))) go to 120
      if(sql.le.sqmax) go to 120
      k = l
      nr = i1+(m-i1)/2
c  if no such interval is found, control is repassed to the driver
c  program (ier = -1).
 120  if(nr.eq.0) go to 190
c  if s(x) is given by the same straight line in two succeeding knot
c  intervals t(l-1) <= x <= t(l) and t(l) <= x <= t(l+1),delete t(l)
      n8 = n-8
      l1 = 0
      if(n8.le.0) go to 150
      do 140 i=1,n8
        if(.not. (bind(i).and.bind(i+1).and.bind(i+2))) go to 140
        l = i+4-l1
        if(k.gt.l) k = k-1
        n = n-1
        l1 = l1+1
        do 130 j=l,n
          t(j) = t(j+1)
 130    continue
 140  continue
c  test whether we cannot further increase the number of knots.
 150  if(n.eq.nmax) go to 180
      if(n.eq.nest) go to 170
c  locate an additional knot at the point x(nr).
      j = n
      do 160 i=k,n
        t(j+1) = t(j)
        j = j-1
 160  continue
      t(k) = x(nr)
      n = n+1
c  restart the computations with the new set of knots.
      go to 40
c  error codes and messages.
 170  ier = -3
      go to 300
 180  ier = -2
      go to 300
 190  ier = -1
      go to 300
 200  ier = 4
      go to 300
 210  ier = 5
 300  return
      end
      subroutine fpcons(iopt,idim,m,u,mx,x,w,ib,ie,k,s,nest,tol,maxit,
     * k1,k2,n,t,nc,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real s,tol,fp
      integer iopt,idim,m,mx,ib,ie,k,nest,maxit,k1,k2,n,nc,ier
c  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),
     * z(nc),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,con1,con4,con9,cos,fac,fpart,fpms,fpold,fp0,f1,f2,f3,
     * half,one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,ui,wi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,jb,je,jj,j1,j2,j3,kbe,
     * l,li,lj,l0,mb,me,mm,new,nk1,nmax,nmin,nn,nplus,npl1,nrint,n8
c  ..local arrays..
      real h(7),xi(10)
c  ..function references
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares curve sinf(u),    c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(u) is the requested curve.                          c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares curve until finally fp<=s.         c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1-max(0,ib-1)-max(0,ie-1)                 c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial curve of   c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the polynomial curve of degree k.           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
c  find which data points are to be concidered.
      mb = 2
      jb = ib
      if(ib.gt.0) go to 10
      mb = 1
      jb = 1
  10  me = m-1
      je = ie
      if(ie.gt.0) go to 20
      me = m
      je = 1
  20  if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      kbe = k1-jb-je
      mmin = kbe+2
      mm = m-mmin
      nmax = nmin+mm
      if(s.gt.0.) go to 40
c  if s=0, s(u) is an interpolating curve.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
      if(mm.eq.0) go to 60
  25  i = k2
      j = 3-jb+k/2
      do 30 l=1,mm
        t(i) = u(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial curve which is a spline curve without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline curve
c  according to the set of knots found at the last call of the routine.
  40  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(u).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = u(1)
          t(i) = u(m)
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline curve
c  sinf(u). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.
c  nn denotes the dimension of the splines
        nn = nk1-ib-ie
c  initialize the b-spline coefficients and the observation matrix a.
        do 75 i=1,nc
          z(i) = 0.
          c(i) = 0.
  75    continue
        if(me.lt.mb) go to 134
        if(nn.eq.0) go to 82
        do 80 i=1,nn
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
  82    l = k1
        jj = (mb-1)*idim
        do 130 it=mb,me
c  fetch the current data point u(it),x(it).
          ui = u(it)
          wi = w(it)
          do 84 j=1,idim
             jj = jj+1
             xi(j) = x(jj)*wi
  84      continue
c  search for knot interval t(l) <= ui < t(l+1).
  86      if(ui.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 86
c  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  90      call fpbspl(t,n,k,ui,l,h)
          do 92 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  92      continue
c  take into account that certain b-spline coefficients must be zero.
          lj = k1
          j = nk1-l-ie
          if(j.ge.0) go to 94
          lj = lj+j
  94      li = 1
          j = l-k1-ib
          if(j.ge.0) go to 96
          li = li-j
          j = 0
  96      if(li.gt.lj) go to 120
c  rotate the new row of the observation matrix into triangle.
          do 110 i=li,lj
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            j1 = j
            do 98 j2 =1,idim
               call fprota(cos,sin,xi(j2),z(j1))
               j1 = j1+n
  98        continue
            if(i.eq.lj) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,lj
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      do 125 j2=1,idim
             fp  = fp+xi(j2)**2
 125      continue
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        if(nn.eq.0) go to 134
        j1 = 1
        do 132 j2=1,idim
           j3 = j1+ib
           call fpback(a,z(j1),nn,k1,c(j3),nest)
           j1 = j1+n
 132    continue
c  test whether the approximation sinf(u) is an acceptable solution.
 134    if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
c  if n = nmax, sinf(u) is an interpolating spline curve.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum of squared residuals for each knot interval
c  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        jj = (mb-1)*idim
        do 180 it=mb,me
          if(u(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 175 j2=1,idim
            fac = 0.
            j1 = l0
            do 170 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 170        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 175      continue
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation
          if(n.eq.nmax) go to 25
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial curve is a
c  solution of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline curve sp(u).          c
c  **********************************************************          c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing curve     c
c  sp(u). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(u) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that f(p),c
c  the sum of squared residuals be = s. we already know that the least c
c  squares kth degree polynomial curve corresponds to p=0, and that    c
c  the least-squares spline curve corresponds to p=infinity. the       c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 252 i=1,nn
         p = p+a(i,1)
 252  continue
      rn = nn
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 255 i=1,nc
          c(i) = z(i)
 255    continue
        do 260 i=1,nn
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 264 i=1,k2
            h(i) = b(it,i)*pinv
 264      continue
          do 268 j=1,idim
            xi(j) = 0.
 268      continue
c  take into account that certain b-spline coefficients must be zero.
          if(it.gt.ib) go to 274
          j1 = ib-it+2
          j2 = 1
          do 270 i=j1,k2
            h(j2) = h(i)
            j2 = j2+1
 270      continue
          do 272 i=j2,k2
            h(i) = 0.
 272      continue
 274      jj = max0(1,it-ib)
          do 290 j=jj,nn
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            j1 = j
            do 277 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 277        continue
            if(j.eq.nn) go to 300
            i2 = min0(nn-j,k1)
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        j1 = 1
        do 308 j2=1,idim
          j3 = j1+ib
          call fpback(g,c(j1),nn,k2,c(j3),nest)
          if(ib.eq.0) go to 306
          j3 = j1
          do 304 i=1,ib
            c(j3) = 0.
            j3 = j3+1
 304      continue
 306      j1 =j1+n
 308    continue
c  computation of f(p).
        fp = 0.
        l = k2
        jj = (mb-1)*idim
        do 330 it=mb,me
          if(u(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 325 j2=1,idim
            fac = 0.
            j1 = l0
            do 320 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 320        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 325      continue
          fp = fp+term*w(it)**2
 330    continue
c  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
      subroutine fpcosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,nm,mb,a,
     * b,const,z,zz,u,q,info,up,left,right,jbind,ibind,ier)
c  ..
c  ..scalar arguments..
      real sq
      integer m,n,maxtr,maxbin,nm,mb,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(n),e(n),c(n),sx(m),a(n,4),b(nm,maxbin),
     * const(n),z(n),zz(n),u(maxbin),q(m,4)
      integer info(maxtr),up(maxtr),left(maxtr),right(maxtr),jbind(mb),
     * ibind(mb)
      logical bind(n)
c  ..local scalars..
      integer count,i,i1,j,j1,j2,j3,k,kdim,k1,k2,k3,k4,k5,k6,
     * l,lp1,l1,l2,l3,merk,nbind,number,n1,n4,n6
      real f,wi,xi
c  ..local array..
      real h(4)
c  ..subroutine references..
c    fpbspl,fpadno,fpdeno,fpfrno,fpseno
c  ..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  if we use the b-spline representation of s(x) our approximation     c
c  problem results in a quadratic programming problem:                 c
c    find the b-spline coefficients c(j),j=1,2,...n-4 such that        c
c        (1) sumi((wi*(yi-sumj(cj*nj(xi))))**2),i=1,2,...m is minimal  c
c        (2) sumj(cj*n''j(t(l+3)))*e(l) <= 0, l=1,2,...n-6.            c
c  to solve this problem we use the theil-van de panne procedure.      c
c  if the inequality constraints (2) are numbered from 1 to n-6,       c
c  this algorithm finds a subset of constraints ibind(1)..ibind(nbind) c
c  such that the solution of the minimization problem (1) with these   c
c  constraints in equality form, satisfies all constraints. such a     c
c  feasible solution is optimal if the lagrange parameters associated  c
c  with that problem with equality constraints, are all positive.      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine n6, the number of inequality constraints.
      n6 = n-6
c  fix the parameters which determine these constraints.
      do 10 i=1,n6
        const(i) = e(i)*(t(i+4)-t(i+1))/(t(i+5)-t(i+2))
  10  continue
c  initialize the triply linked tree which is used to find the subset
c  of constraints ibind(1),...ibind(nbind).
      count = 1
      info(1) = 0
      left(1) = 0
      right(1) = 0
      up(1) = 1
      merk = 1
c  set up the normal equations  n'nc=n'y  where n denotes the m x (n-4)
c  observation matrix with elements ni,j = wi*nj(xi)  and y is the
c  column vector with elements yi*wi.
c  from the properties of the b-splines nj(x),j=1,2,...n-4, it follows
c  that  n'n  is a (n-4) x (n-4)  positive definit bandmatrix of
c  bandwidth 7. the matrices n'n and n'y are built up in a and z.
      n4 = n-4
c  initialization
      do 20 i=1,n4
        z(i) = 0.
        do 20 j=1,4
          a(i,j) = 0.
  20  continue
      l = 4
      lp1 = l+1
      do 70 i=1,m
c  fetch the current row of the observation matrix.
        xi = x(i)
        wi = w(i)**2
c  search for knot interval  t(l) <= xi < t(l+1)
  30    if(xi.lt.t(lp1) .or. l.eq.n4) go to 40
        l = lp1
        lp1 = l+1
        go to 30
c  evaluate the four non-zero cubic b-splines nj(xi),j=l-3,...l.
  40    call fpbspl(t,n,3,xi,l,h)
c  store in q these values h(1),h(2),...h(4).
        do 50 j=1,4
          q(i,j) = h(j)
  50    continue
c  add the contribution of the current row of the observation matrix
c  n to the normal equations.
        l3 = l-3
        k1 = 0
        do 60 j1 = l3,l
          k1 = k1+1
          f = h(k1)
          z(j1) = z(j1)+f*wi*y(i)
          k2 = k1
          j2 = 4
          do 60 j3 = j1,l
            a(j3,j2) = a(j3,j2)+f*wi*h(k2)
            k2 = k2+1
            j2 = j2-1
  60    continue
  70  continue
c  since n'n is a symmetric matrix it can be factorized as
c        (3)  n'n = (r1)'(d1)(r1)
c  with d1 a diagonal matrix and r1 an (n-4) x (n-4)  unit upper
c  triangular matrix of bandwidth 4. the matrices r1 and d1 are built
c  up in a. at the same time we solve the systems of equations
c        (4)  (r1)'(z2) = n'y
c        (5)  (d1) (z1) = (z2)
c  the vectors z2 and z1 are kept in zz and z.
      do 140 i=1,n4
        k1 = 1
        if(i.lt.4) k1 = 5-i
        k2 = i-4+k1
        k3 = k2
        do 100 j=k1,4
          k4 = j-1
          k5 = 4-j+k1
          f = a(i,j)
          if(k1.gt.k4) go to 90
          k6 = k2
          do 80 k=k1,k4
            f = f-a(i,k)*a(k3,k5)*a(k6,4)
            k5 = k5+1
            k6 = k6+1
  80      continue
  90      if(j.eq.4) go to 110
          a(i,j) = f/a(k3,4)
          k3 = k3+1
 100    continue
 110    a(i,4) = f
        f = z(i)
        if(i.eq.1) go to 130
        k4 = i
        do 120 j=k1,3
          k = k1+3-j
          k4 = k4-1
          f = f-a(i,k)*z(k4)*a(k4,4)
 120    continue
 130    z(i) = f/a(i,4)
        zz(i) = f
 140  continue
c  start computing the least-squares cubic spline without taking account
c  of any constraint.
      nbind = 0
      n1 = 1
      ibind(1) = 0
c  main loop for the least-squares problems with different subsets of
c  the constraints (2) in equality form. the resulting b-spline coeff.
c  c and lagrange parameters u are the solution of the system
c            ! n'n  b' ! ! c !   ! n'y !
c        (6) !         ! !   ! = !     !
c            !  b   0  ! ! u !   !  0  !
c  z1 is stored into array c.
 150  do 160 i=1,n4
        c(i) = z(i)
 160  continue
c  if there are no equality constraints, compute the coeff. c directly.
      if(nbind.eq.0) go to 370
c  initialization
      kdim = n4+nbind
      do 170 i=1,nbind
        do 170 j=1,kdim
          b(j,i) = 0.
 170  continue
c  matrix b is built up,expressing that the constraints nrs ibind(1),...
c  ibind(nbind) must be satisfied in equality form.
      do 180 i=1,nbind
        l = ibind(i)
        b(l,i) = e(l)
        b(l+1,i) = -(e(l)+const(l))
        b(l+2,i) = const(l)
 180  continue
c  find the matrix (b1) as the solution of the system of equations
c        (7)  (r1)'(d1)(b1) = b'
c  (b1) is built up in the upper part of the array b(rows 1,...n-4).
      do 220 k1=1,nbind
        l = ibind(k1)
        do 210 i=l,n4
          f = b(i,k1)
          if(i.eq.1) go to 200
          k2 = 3
          if(i.lt.4) k2 = i-1
          do 190 k3=1,k2
            l1 = i-k3
            l2 = 4-k3
            f = f-b(l1,k1)*a(i,l2)*a(l1,4)
 190      continue
 200      b(i,k1) = f/a(i,4)
 210    continue
 220  continue
c  factorization of the symmetric matrix  -(b1)'(d1)(b1)
c        (8)  -(b1)'(d1)(b1) = (r2)'(d2)(r2)
c  with (d2) a diagonal matrix and (r2) an nbind x nbind unit upper
c  triangular matrix. the matrices r2 and d2 are built up in the lower
c  part of the array b (rows n-3,n-2,...n-4+nbind).
      do 270 i=1,nbind
        i1 = i-1
        do 260 j=i,nbind
          f = 0.
          do 230 k=1,n4
            f = f+b(k,i)*b(k,j)*a(k,4)
 230      continue
          k1 = n4+1
          if(i1.eq.0) go to 250
          do 240 k=1,i1
            f = f+b(k1,i)*b(k1,j)*b(k1,k)
            k1 = k1+1
 240      continue
 250      b(k1,j) = -f
          if(j.eq.i) go to 260
          b(k1,j) = b(k1,j)/b(k1,i)
 260    continue
 270  continue
c  according to (3),(7) and (8) the system of equations (6) becomes
c         ! (r1)'    0  ! ! (d1)    0  ! ! (r1)  (b1) ! ! c !   ! n'y !
c    (9)  !             ! !            ! !            ! !   ! = !     !
c         ! (b1)'  (r2)'! !   0   (d2) ! !   0   (r2) ! ! u !   !  0  !
c  backward substitution to obtain the b-spline coefficients c(j),j=1,..
c  n-4 and the lagrange parameters u(j),j=1,2,...nbind.
c  first step of the backward substitution: solve the system
c             ! (r1)'(d1)      0     ! ! (c1) !   ! n'y !
c        (10) !                      ! !      ! = !     !
c             ! (b1)'(d1)  (r2)'(d2) ! ! (u1) !   !  0  !
c  from (4) and (5) we know that this is equivalent to
c        (11)  (c1) = (z1)
c        (12)  (r2)'(d2)(u1) = -(b1)'(z2)
      do 310 i=1,nbind
        f = 0.
        do 280 j=1,n4
          f = f+b(j,i)*zz(j)
 280    continue
        i1 = i-1
        k1 = n4+1
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          f = f+u(j)*b(k1,i)*b(k1,j)
          k1 = k1+1
 290    continue
 300    u(i) = -f/b(k1,i)
 310  continue
c  second step of the backward substitution: solve the system
c             ! (r1)  (b1) ! ! c !   ! c1 !
c        (13) !            ! !   ! = !    !
c             !   0   (r2) ! ! u !   ! u1 !
      k1 = nbind
      k2 = kdim
c  find the lagrange parameters u.
      do 340 i=1,nbind
        f = u(k1)
        if(i.eq.1) go to 330
        k3 = k1+1
        do 320 j=k3,nbind
          f = f-u(j)*b(k2,j)
 320    continue
 330    u(k1) = f
        k1 = k1-1
        k2 = k2-1
 340  continue
c  find the b-spline coefficients c.
      do 360 i=1,n4
        f = c(i)
        do 350 j=1,nbind
          f = f-u(j)*b(i,j)
 350    continue
        c(i) = f
 360  continue
 370  k1 = n4
      do 390 i=2,n4
        k1 = k1-1
        f = c(k1)
        k2 = 1
        if(i.lt.5) k2 = 5-i
        k3 = k1
        l = 3
        do 380 j=k2,3
          k3 = k3+1
          f = f-a(k3,l)*c(k3)
          l = l-1
 380    continue
        c(k1) = f
 390  continue
c  test whether the solution of the least-squares problem with the
c  constraints ibind(1),...ibind(nbind) in equality form, satisfies
c  all of the constraints (2).
      k = 1
c  number counts the number of violated inequality constraints.
      number = 0
      do 440 j=1,n6
        l = ibind(k)
        k = k+1
        if(j.eq.l) go to 440
        k = k-1
c  test whether constraint j is satisfied
        f = e(j)*(c(j)-c(j+1))+const(j)*(c(j+2)-c(j+1))
        if(f.le.0.) go to 440
c  if constraint j is not satisfied, add a branch of length nbind+1
c  to the tree. the nodes of this branch contain in their information
c  field the number of the constraints ibind(1),...ibind(nbind) and j,
c  arranged in increasing order.
        number = number+1
        k1 = k-1
        if(k1.eq.0) go to 410
        do 400 i=1,k1
          jbind(i) = ibind(i)
 400    continue
 410    jbind(k) = j
        if(l.eq.0) go to 430
        do 420 i=k,nbind
          jbind(i+1) = ibind(i)
 420    continue
 430    call fpadno(maxtr,up,left,right,info,count,merk,jbind,n1,ier)
c  test whether the storage space which is required for the tree,exceeds
c  the available storage space.
        if(ier.ne.0) go to 560
 440  continue
c  test whether the solution of the least-squares problem with equality
c  constraints is a feasible solution.
      if(number.eq.0) go to 470
c  test whether there are still cases with nbind constraints in
c  equality form to be considered.
 450  if(merk.gt.1) go to 460
      nbind = n1
c  test whether the number of knots where s''(x)=0 exceeds maxbin.
      if(nbind.gt.maxbin) go to 550
      n1 = n1+1
      ibind(n1) = 0
c  search which cases with nbind constraints in equality form
c  are going to be considered.
      call fpdeno(maxtr,up,left,right,nbind,merk)
c  test whether the quadratic programming problem has a solution.
      if(merk.eq.1) go to 570
c  find a new case with nbind constraints in equality form.
 460  call fpseno(maxtr,up,left,right,info,merk,ibind,nbind)
      go to 150
c  test whether the feasible solution is optimal.
 470  ier = 0
      do 480 i=1,n6
        bind(i) = .false.
 480  continue
      if(nbind.eq.0) go to 500
      do 490 i=1,nbind
        if(u(i).le.0.) go to 450
        j = ibind(i)
        bind(j) = .true.
 490  continue
c  evaluate s(x) at the data points x(i) and calculate the weighted
c  sum of squared residual right hand sides sq.
 500  sq = 0.
      l = 4
      lp1 = 5
      do 530 i=1,m
 510    if(x(i).lt.t(lp1) .or. l.eq.n4) go to 520
        l = lp1
        lp1 = l+1
        go to 510
 520    sx(i) = c(l-3)*q(i,1)+c(l-2)*q(i,2)+c(l-1)*q(i,3)+c(l)*q(i,4)
        sq = sq+(w(i)*(y(i)-sx(i)))**2
 530  continue
      go to 600
c  error codes and messages.
 550  ier = 1
      go to 600
 560  ier = 2
      go to 600
 570  ier = 3
 600  return
      end
      subroutine fpcsin(a,b,par,sia,coa,sib,cob,ress,resc)
c  fpcsin calculates the integrals ress=integral((b-x)**3*sin(par*x))
c  and resc=integral((b-x)**3*cos(par*x)) over the interval (a,b),
c  given sia=sin(par*a),coa=cos(par*a),sib=sin(par*b) and cob=cos(par*b)
c  ..
c  ..scalar arguments..
      real a,b,par,sia,coa,sib,cob,ress,resc
c  ..local scalars..
      integer i,j
      real ab,ab4,ai,alfa,beta,b2,b4,eps,fac,f1,f2,one,quart,six,
     * three,two
c  ..function references..
      real abs
c  ..
      one = 0.1e+01
      two = 0.2e+01
      three = 0.3e+01
      six = 0.6e+01
      quart = 0.25e+0
      eps = 0.1e-09
      ab = b-a
      ab4 = ab**4
      alfa = ab*par
c the way of calculating the integrals ress and resc depends on
c the value of alfa = (b-a)*par.
      if(abs(alfa).le.one) go to 100
c integration by parts.
      beta = one/alfa
      b2 = beta**2
      b4 = six*b2**2
      f1 = three*b2*(one-two*b2)
      f2 = beta*(one-six*b2)
      ress = ab4*(coa*f2+sia*f1+sib*b4)
      resc = ab4*(coa*f1-sia*f2+cob*b4)
      go to 400
c ress and resc are found by evaluating a series expansion.
 100  fac = quart
      f1 = fac
      f2 = 0.
      i = 4
      do 200 j=1,5
        i = i+1
        ai = i
        fac = fac*alfa/ai
        f2 = f2+fac
        if(abs(fac).le.eps) go to 300
        i = i+1
        ai = i
        fac = -fac*alfa/ai
        f1 = f1+fac
        if(abs(fac).le.eps) go to 300
 200  continue
 300  ress = ab4*(coa*f2+sia*f1)
      resc = ab4*(coa*f1-sia*f2)
 400  return
      end
      subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,
     * n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real xb,xe,s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),fpint(nest),
     * z(nest),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,con1,con4,con9,cos,half,fpart,fpms,fpold,fp0,f1,f2,f3,
     * one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,wi,xi,yi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,k3,l,l0,
     * mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
c  ..local arrays..
      real h(7)
c  ..function references
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares spline sinf(x),   c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(x) is the requested approximation.                  c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares polynomial of degree k.   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.) go to 45
c  if s=0, s(x) is an interpolating spline.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = x(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree k which is a spline without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline
c  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = xb
          t(i) = xe
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline
c  sinf(x). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.
c  initialize the observation matrix a.
        do 80 i=1,nk1
          z(i) = 0.
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
        l = k1
        do 130 it=1,m
c  fetch the current data point x(it),y(it).
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  85      if(xi.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  90      call fpbspl(t,n,k,xi,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
c  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      fp = fp+yi**2
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        call fpback(a,z,nk1,k1,c,nest)
c  test whether the approximation sinf(x) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
c  if n = nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
c  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        do 180 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 170 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 170      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 10
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial is a solution
c  of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline sp(x).                c
c  ***************************************************                 c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
c  the least-squares kth degree polynomial corresponds to p=0, and     c
c  that the least-squares spline corresponds to p=infinity. the        c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 255 i=1,nk1
         p = p+a(i,1)
 255  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 260 i=1,nk1
          c(i) = z(i)
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          yi = 0.
          do 290 j=it,nk1
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,c(j))
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        call fpback(g,c,nk1,k2,c,nest)
c  computation of f(p).
        fp = 0.
        l = k2
        do 330 it=1,m
          if(x(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 320 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 320      continue
          fp = fp+(w(it)*(term-y(it)))**2
 330    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
      subroutine fpcuro(a,b,c,d,x,n)
c  subroutine fpcuro finds the real zeros of a cubic polynomial
c  p(x) = a*x**3+b*x**2+c*x+d.
c
c  calling sequence:
c     call fpcuro(a,b,c,d,x,n)
c
c  input parameters:
c    a,b,c,d: real values, containing the coefficients of p(x).
c
c  output parameters:
c    x      : real array,length 3, which contains the real zeros of p(x)
c    n      : integer, giving the number of real zeros of p(x).
c  ..
c  ..scalar arguments..
      real a,b,c,d
      integer n
c  ..array argument..
      real x(3)
c  ..local scalars..
      integer i
      real a1,b1,c1,df,disc,d1,e3,f,four,half,ovfl,pi3,p3,q,r,
     * step,tent,three,two,u,u1,u2,y
c  ..function references..
      real abs,amax1,atan,atan2,cos,sign,sqrt
c  set constants
      two = 0.2e+01
      three = 0.3e+01
      four = 0.4e+01
      ovfl =0.1e+05
      half = 0.5e+0
      tent = 0.1e+0
      e3 = tent/0.3e0
      pi3 = atan(0.1e+01)/0.75e0
      a1 = abs(a)
      b1 = abs(b)
      c1 = abs(c)
      d1 = abs(d)
c  test whether p(x) is a third degree polynomial.
      if(amax1(b1,c1,d1).lt.a1*ovfl) go to 300
c  test whether p(x) is a second degree polynomial.
      if(amax1(c1,d1).lt.b1*ovfl) go to 200
c  test whether p(x) is a first degree polynomial.
      if(d1.lt.c1*ovfl) go to 100
c  p(x) is a constant function.
      n = 0
      go to 800
c  p(x) is a first degree polynomial.
 100  n = 1
      x(1) = -d/c
      go to 500
c  p(x) is a second degree polynomial.
 200  disc = c*c-four*b*d
      n = 0
      if(disc.lt.0.) go to 800
      n = 2
      u = sqrt(disc)
      b1 = b+b
      x(1) = (-c+u)/b1
      x(2) = (-c-u)/b1
      go to 500
c  p(x) is a third degree polynomial.
 300  b1 = b/a*e3
      c1 = c/a
      d1 = d/a
      q = c1*e3-b1*b1
      r = b1*b1*b1+(d1-b1*c1)*half
      disc = q*q*q+r*r
      if(disc.gt.0.) go to 400
      u = sqrt(abs(q))
      if(r.lt.0.) u = -u
      p3 = atan2(sqrt(-disc),abs(r))*e3
      u2 = u+u
      n = 3
      x(1) = -u2*cos(p3)-b1
      x(2) = u2*cos(pi3-p3)-b1
      x(3) = u2*cos(pi3+p3)-b1
      go to 500
 400  u = sqrt(disc)
      u1 = -r+u
      u2 = -r-u
      n = 1
      x(1) = sign(abs(u1)**e3,u1)+sign(abs(u2)**e3,u2)-b1
c  apply a newton iteration to improve the accuracy of the roots.
 500  do 700 i=1,n
        y = x(i)
        f = ((a*y+b)*y+c)*y+d
        df = (three*a*y+two*b)*y+c
        step = 0.
        if(abs(f).lt.abs(df)*tent) step = f/df
        x(i) = y-step
 700  continue
 800  return
      end
      subroutine fpcyt1(a,n,nn)
c (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
c elements stored as follows
c
c    | a(1,2) a(1,3)                                    a(1,1)  |
c    | a(2,1) a(2,2) a(2,3)                                     |
c    |        a(3,1) a(3,2) a(3,3)                              |
c    |               ...............                            |
c    |                               a(n-1,1) a(n-1,2) a(n-1,3) |
c    | a(n,3)                                  a(n,1)   a(n,2)  |
c
c  ..
c  ..scalar arguments..
      integer n,nn
c  ..array arguments..
      real a(nn,6)
c  ..local scalars..
      real aa,beta,gamma,sum,teta,v,one
      integer i,n1,n2
c  ..
c  set constant
      one = 1
      n2 = n-2
      beta = one/a(1,2)
      gamma = a(n,3)
      teta = a(1,1)*beta
      a(1,4) = beta
      a(1,5) = gamma
      a(1,6) = teta
      sum = gamma*teta
      do 10 i=2,n2
         v = a(i-1,3)*beta
         aa = a(i,1)
         beta = one/(a(i,2)-aa*v)
         gamma = -gamma*v
         teta = -teta*aa*beta
         a(i,4) = beta
         a(i,5) = gamma
         a(i,6) = teta
         sum = sum+gamma*teta
  10  continue
      n1 = n-1
      v = a(n2,3)*beta
      aa = a(n1,1)
      beta = one/(a(n1,2)-aa*v)
      gamma = a(n,1)-gamma*v
      teta = (a(n1,3)-teta*aa)*beta
      a(n1,4) = beta
      a(n1,5) = gamma
      a(n1,6) = teta
      a(n,4) = one/(a(n,2)-(sum+gamma*teta))
      return
      end
      subroutine fpcyt2(a,n,b,c,nn)
c subroutine fpcyt2 solves a linear n x n system
c         a * c = b
c where matrix a is a cyclic tridiagonal matrix, decomposed
c using subroutine fpsyt1.
c  ..
c  ..scalar arguments..
      integer n,nn
c  ..array arguments..
      real a(nn,6),b(n),c(n)
c  ..local scalars..
      real cc,sum
      integer i,j,j1,n1
c  ..
      c(1) = b(1)*a(1,4)
      sum = c(1)*a(1,5)
      n1 = n-1
      do 10 i=2,n1
         c(i) = (b(i)-a(i,1)*c(i-1))*a(i,4)
         sum = sum+c(i)*a(i,5)
  10  continue
      cc = (b(n)-sum)*a(n,4)
      c(n) = cc
      c(n1) = c(n1)-cc*a(n1,6)
      j = n1
      do 20 i=3,n
         j1 = j-1
         c(j1) = c(j1)-c(j)*a(j1,3)*a(j1,4)-cc*a(j1,6)
         j = j1
  20  continue
      return
      end
      subroutine fpdeno(maxtr,up,left,right,nbind,merk)
c  subroutine fpdeno frees the nodes of all branches of a triply linked
c  tree with length < nbind by putting to zero their up field.
c  on exit the parameter merk points to the terminal node of the
c  most left branch of length nbind or takes the value 1 if there
c  is no such branch.
c  ..
c  ..scalar arguments..
      integer maxtr,nbind,merk
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr)
c  ..local scalars ..
      integer i,j,k,l,niveau,point
c  ..
      i = 1
      niveau = 0
  10  point = i
      i = left(point)
      if(i.eq.0) go to 20
      niveau = niveau+1
      go to 10
  20  if(niveau.eq.nbind) go to 70
  30  i = right(point)
      j = up(point)
      up(point) = 0
      k = left(j)
      if(point.ne.k) go to 50
      if(i.ne.0) go to 40
      niveau = niveau-1
      if(niveau.eq.0) go to 80
      point = j
      go to 30
  40  left(j) = i
      go to 10
  50  l = right(k)
      if(point.eq.l) go to 60
      k = l
      go to 50
  60  right(k) = i
      point = k
  70  i = right(point)
      if(i.ne.0) go to 10
      i = up(point)
      niveau = niveau-1
      if(niveau.eq.0) go to 80
      point = i
      go to 70
  80  k = 1
      l = left(k)
      if(up(l).eq.0) return
  90  merk = k
      k = left(k)
      if(k.ne.0) go to 90
      return
      end
      subroutine fpdisc(t,n,k2,b,nest)
c  subroutine fpdisc calculates the discontinuity jumps of the kth
c  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
c  ..scalar arguments..
      integer n,k2,nest
c  ..array arguments..
      real t(n),b(nest,k2)
c  ..local scalars..
      real an,fac,prod
      integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
c  ..local array..
      real h(12)
c  ..
      k1 = k2-1
      k = k1-1
      nk1 = n-k1
      nrint = nk1-k
      an = nrint
      fac = an/(t(nk1+1)-t(k1))
      do 40 l=k2,nk1
        lmk = l-k1
        do 10 j=1,k1
          ik = j+k1
          lj = l+j
          lk = lj-k2
          h(j) = t(l)-t(lk)
          h(ik) = t(l)-t(lj)
  10    continue
        lp = lmk
        do 30 j=1,k2
          jk = j
          prod = h(j)
          do 20 i=1,k
            jk = jk+1
            prod = prod*h(jk)*fac
  20      continue
          lk = lp+k1
          b(lmk,j) = (t(lk)-t(lp))/prod
          lp = lp+1
  30    continue
  40  continue
      return
      end
      subroutine fpfrno(maxtr,up,left,right,info,point,merk,n1,
     * count,ier)
c  subroutine fpfrno collects the free nodes (up field zero) of the
c  triply linked tree the information of which is kept in the arrays
c  up,left,right and info. the maximal length of the branches of the
c  tree is given by n1. if no free nodes are found, the error flag
c  ier is set to 1.
c  ..
c  ..scalar arguments..
      integer maxtr,point,merk,n1,count,ier
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr)
c  ..local scalars
      integer i,j,k,l,n,niveau
c  ..
      ier = 1
      if(n1.eq.2) go to 140
      niveau = 1
      count = 2
  10  j = 0
      i = 1
  20  if(j.eq.niveau) go to 30
      k = 0
      l = left(i)
      if(l.eq.0) go to 110
      i = l
      j = j+1
      go to 20
  30  if(i-count) 110,100,40
  40  if(up(count).eq.0) go to 50
      count = count+1
      go to 30
  50  up(count) = up(i)
      left(count) = left(i)
      right(count) = right(i)
      info(count) = info(i)
      if(merk.eq.i) merk = count
      if(point.eq.i) point = count
      if(k.eq.0) go to 60
      right(k) = count
      go to 70
  60  n = up(i)
      left(n) = count
  70  l = left(i)
  80  if(l.eq.0) go to 90
      up(l) = count
      l = right(l)
      go to 80
  90  up(i) = 0
      i = count
 100  count = count+1
 110  l = right(i)
      k = i
      if(l.eq.0) go to 120
      i = l
      go to 20
 120  l = up(i)
      j = j-1
      if(j.eq.0) go to 130
      i = l
      go to 110
 130  niveau = niveau+1
      if(niveau.le.n1) go to 10
      if(count.gt.maxtr) go to 140
      ier = 0
 140  return
      end
      subroutine fpgivs(piv,ww,cos,sin)
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real piv,ww,cos,sin
c  ..local scalars..
      real dd,one,store
c  ..function references..
      real abs,sqrt
c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
      subroutine fpgrdi(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v,mv,z,mz,dz,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,spu,spv,
     * right,q,au,av1,av2,bu,bv,aa,bb,cc,cosi,nru,nrv)
c  ..
c  ..scalar arguments..
      real p,sq,fp
      integer ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mz,iop0,iop1,nu,nv,nc,
     * mm,mvnu
c  ..array arguments..
      real u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),
     * spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),av1(nv,6),
     * av2(nv,4),aa(2,mv),bb(2,nv),cc(nv),cosi(2,nv),bu(nu,5),bv(nv,5)
      integer nru(mu),nrv(mv)
c  ..local scalars..
      real arg,co,dz1,dz2,dz3,fac,fac0,pinv,piv,si,term,one,three,half
      integer i,ic,ii,ij,ik,iq,irot,it,iz,i0,i1,i2,i3,j,jj,jk,jper,
     * j0,j1,k,k1,k2,l,l0,l1,l2,mvv,ncof,nrold,nroldu,nroldv,number,
     * numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
c  ..local arrays..
      real h(5),h1(5),h2(4)
c  ..function references..
      integer min0
      real cos,sin
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota
c  ..
c  let
c               |   (spu)    |            |   (spv)    |
c        (au) = | ---------- |     (av) = | ---------- |
c               | (1/p) (bu) |            | (1/p) (bv) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
c                coefficients.
c       z      : the mu x mv matrix which contains the function values.
c       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
c                according to the least-squares problems in the u-,resp.
c                v-direction.
c       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
c                containing the discontinuity jumps of the derivatives
c                of the b-splines in the u-,resp.v-variable at the knots
c  the b-spline coefficients of the smoothing spline are then calculated
c  as the least-squares solution of the following over-determined linear
c  system of equations
c
c    (1)  (av) c (au)' = q
c
c  subject to the constraints
c
c    (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
c
c    (3)  if iop0 = 0  c(1,j) = dz(1)
c            iop0 = 1  c(1,j) = dz(1)
c                      c(2,j) = dz(1)+(dz(2)*cosi(1,j)+dz(3)*cosi(2,j))*
c                               tu(5)/3. = cc(j) , j=1,2,...nv-4
c
c    (4)  if iop1 = 1  c(nu-4,j) = 0, j=1,2,...,nv-4.
c
c  set constants
      one = 1
      three = 3
      half = 0.5
c  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-1
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv and iop0 and
c  on the value of p whether the matrices (spu), (spv), (bu), (bv) and
c  (cosi) still must be determined.
      if(ifsu.ne.0) go to 30
c  calculate the non-zero elements of the matrix (spu) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 25 it=1,mu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 15
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  15    call fpbspl(tu,nu,3,arg,l,h)
        do 20 i=1,4
          spu(it,i) = h(i)
  20    continue
        nru(it) = number
  25  continue
      ifsu = 1
c  calculate the non-zero elements of the matrix (spv) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the v-direction.
  30  if(ifsv.ne.0) go to 85
      l = 4
      l1 = 5
      number = 0
      do 50 it=1,mv
        arg = v(it)
  35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40
        l = l1
        l1 = l+1
        number = number+1
        go to 35
  40    call fpbspl(tv,nv,3,arg,l,h)
        do 45 i=1,4
          spv(it,i) = h(i)
  45    continue
        nrv(it) = number
  50  continue
      ifsv = 1
      if(iop0.eq.0) go to 85
c  calculate the coefficients of the interpolating splines for cos(v)
c  and sin(v).
      do 55 i=1,nv4
         cosi(1,i) = 0.
         cosi(2,i) = 0.
  55  continue
      if(nv7.lt.4) go to 85
      do 65 i=1,nv7
         l = i+3
         arg = tv(l)
         call fpbspl(tv,nv,3,arg,l,h)
         do 60 j=1,3
            av1(i,j) = h(j)
  60     continue
         cosi(1,i) = cos(arg)
         cosi(2,i) = sin(arg)
  65  continue
      call fpcyt1(av1,nv7,nv)
      do 80 j=1,2
         do 70 i=1,nv7
            right(i) = cosi(j,i)
  70     continue
         call fpcyt2(av1,nv7,right,right,nv)
         do 75 i=1,nv7
            cosi(j,i+1) = right(i)
  75     continue
         cosi(j,1) = cosi(j,nv7+1)
         cosi(j,nv7+2) = cosi(j,2)
         cosi(j,nv4) = cosi(j,3)
  80  continue
  85  if(p.le.0.) go to  150
c  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 90
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
c  calculate the non-zero elements of the matrix (bv).
  90  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
c  substituting (2),(3) and (4) into (1), we obtain the overdetermined
c  system
c         (5)  (avv) (cr) (auu)' = (qq)
c  from which the nuu*nv7 remaining coefficients
c         c(i,j) , i=2+iop0,3+iop0,...,nu-4-iop1 ; j=1,2,...,nv-7 ,
c  the elements of (cr), are then determined in the least-squares sense.
c  simultaneously, we compute the resulting sum of squared residuals sq.
 150  dz1 = dz(1)
      do 155 i=1,mv
         aa(1,i) = dz1
 155  continue
      if(nv8.eq.0 .or. p.le.0.) go to 165
      do 160 i=1,nv8
         bb(1,i) = 0.
 160  continue
 165  mvv = mv
      if(iop0.eq.0) go to 220
      fac = tu(5)/three
      dz2 = dz(2)*fac
      dz3 = dz(3)*fac
      do 170 i=1,nv4
         cc(i) = dz1+dz2*cosi(1,i)+dz3*cosi(2,i)
 170  continue
      do 190 i=1,mv
         number = nrv(i)
         fac = 0.
         do 180 j=1,4
            number = number+1
            fac = fac+cc(number)*spv(i,j)
 180     continue
         aa(2,i) = fac
 190  continue
      if(nv8.eq.0 .or. p.le.0.) go to 220
      do 210 i=1,nv8
         number = i
         fac = 0.
         do 200 j=1,5
            fac = fac+cc(number)*bv(i,j)
            number = number+1
 200     continue
         bb(2,i) = fac*pinv
 210  continue
      mvv = mvv+nv8
c  we first determine the matrices (auu) and (qq). then we reduce the
c  matrix (auu) to upper triangular form (ru) using givens rotations.
c  we apply the same transformations to the rows of matrix qq to obtain
c  the (mv+nv8) x nuu matrix g.
c  we store matrix (ru) into au and g into q.
 220  l = mvv*nuu
c  initialization.
      sq = 0.
      do 230 i=1,l
        q(i) = 0.
 230  continue
      do 240 i=1,nuu
        do 240 j=1,5
          au(i,j) = 0.
 240  continue
      l = 0
      nrold = 0
      n1 = nrold+1
      do 420 it=1,mu
        number = nru(it)
c  find the appropriate column of q.
 250    do 260 j=1,mvv
           right(j) = 0.
 260    continue
        if(nrold.eq.number) go to 280
        if(p.le.0.) go to 410
c  fetch a new row of matrix (bu).
        do 270 j=1,5
          h(j) = bu(n1,j)*pinv
 270    continue
        i0 = 1
        i1 = 5
        go to 310
c  fetch a new row of matrix (spu).
 280    do 290 j=1,4
          h(j) = spu(it,j)
 290    continue
c  find the appropriate column of q.
        do 300 j=1,mv
          l = l+1
          right(j) = z(l)
 300    continue
        i0 = 1
        i1 = 4
 310    if(nu7-number .eq. iop1) i1 = i1-1
        j0 = n1
c  take into account that we eliminate the constraints (3)
 320     if(j0-1.gt.iop0) go to 360
         fac0 = h(i0)
         do 330 j=1,mv
            right(j) = right(j)-fac0*aa(j0,j)
 330     continue
         if(mv.eq.mvv) go to 350
         j = mv
         do 340 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac0*bb(j0,jj)
 340     continue
 350     j0 = j0+1
         i0 = i0+1
         go to 320
 360     irot = nrold-iop0-1
         if(irot.lt.0) irot = 0
c  rotate the new row of matrix (auu) into triangle.
        do 390 i=i0,i1
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,au(irot,1),co,si)
c  apply that transformation to the rows of matrix (qq).
          iq = (irot-1)*mvv
          do 370 j=1,mvv
            iq = iq+1
            call fprota(co,si,right(j),q(iq))
 370      continue
c  apply that transformation to the columns of (auu).
          if(i.eq.i1) go to 390
          i2 = 1
          i3 = i+1
          do 380 j=i3,i1
            i2 = i2+1
            call fprota(co,si,h(j),au(irot,i2))
 380      continue
 390    continue
c we update the sum of squared residuals
        do 395 j=1,mvv
          sq = sq+right(j)**2
 395    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = n1
        n1 = n1+1
        go to 250
 420  continue
c  we determine the matrix (avv) and then we reduce her to
c  upper triangular form (rv) using givens rotations.
c  we apply the same transformations to the columns of matrix
c  g to obtain the (nv-7) x (nu-5-iop0-iop1) matrix h.
c  we store matrix (rv) into av1 and av2, h into c.
c  the nv7 x nv7 upper triangular matrix (rv) has the form
c              | av1 '     |
c       (rv) = |     ' av2 |
c              |  0  '     |
c  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 upper
c  triangular matrix of bandwidth 5.
      ncof = nuu*nv7
c  initialization.
      do 430 i=1,ncof
        c(i) = 0.
 430  continue
      do 440 i=1,nv4
        av1(i,5) = 0.
        do 440 j=1,4
          av1(i,j) = 0.
          av2(i,j) = 0.
 440  continue
      jper = 0
      nrold = 0
      do 770 it=1,mv
        number = nrv(it)
 450    if(nrold.eq.number) go to 480
        if(p.le.0.) go to 760
c  fetch a new row of matrix (bv).
        n1 = nrold+1
        do 460 j=1,5
          h(j) = bv(n1,j)*pinv
 460    continue
c  find the appropiate row of g.
        do 465 j=1,nuu
          right(j) = 0.
 465    continue
        if(mv.eq.mvv) go to 510
        l = mv+n1
        do 470 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 470    continue
        go to 510
c  fetch a new row of matrix (spv)
 480    h(5) = 0.
        do 490 j=1,4
          h(j) = spv(it,j)
 490    continue
c  find the appropiate row of g.
        l = it
        do 500 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 500    continue
c  test whether there are non-zero values in the new row of (avv)
c  corresponding to the b-splines n(j,v),j=nv7+1,...,nv4.
 510     if(nrold.lt.nv11) go to 710
         if(jper.ne.0) go to 550
c  initialize the matrix (av2).
         jk = nv11+1
         do 540 i=1,4
            ik = jk
            do 520 j=1,5
               if(ik.le.0) go to 530
               av2(ik,i) = av1(ik,j)
               ik = ik-1
 520        continue
 530        jk = jk+1
 540     continue
         jper = 1
c  if one of the non-zero elements of the new row corresponds to one of
c  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
c  (2) for setting up this row of (avv). the row is stored in h1( the
c  part with respect to av1) and h2 (the part with respect to av2).
 550     do 560 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 560     continue
         h1(5) = 0.
         j = nrold-nv11
         do 600 i=1,5
            j = j+1
            l0 = j
 570        l1 = l0-4
            if(l1.le.0) go to 590
            if(l1.le.nv11) go to 580
            l0 = l1-nv11
            go to 570
 580        h1(l1) = h(i)
            go to 600
 590        h2(l0) = h2(l0) + h(i)
 600     continue
c  rotate the new row of (avv) into triangle.
         if(nv11.le.0) go to 670
c  rotations with the rows 1,2,...,nv11 of (avv).
         do 660 j=1,nv11
            piv = h1(1)
            i2 = min0(nv11-j,4)
            if(piv.eq.0.) go to 640
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(j,1),co,si)
c  apply that transformation to the columns of matrix g.
            ic = j
            do 610 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 610        continue
c  apply that transformation to the rows of (avv) with respect to av2.
            do 620 i=1,4
               call fprota(co,si,h2(i),av2(j,i))
 620        continue
c  apply that transformation to the rows of (avv) with respect to av1.
            if(i2.eq.0) go to 670
            do 630 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),av1(j,i1))
 630        continue
 640        do 650 i=1,i2
               h1(i) = h1(i+1)
 650        continue
            h1(i2+1) = 0.
 660     continue
c  rotations with the rows nv11+1,...,nv7 of avv.
 670     do 700 j=1,4
            ij = nv11+j
            if(ij.le.0) go to 700
            piv = h2(j)
            if(piv.eq.0.) go to 700
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av2(ij,j),co,si)
c  apply that transformation to the columns of matrix g.
            ic = ij
            do 680 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 680        continue
            if(j.eq.4) go to 700
c  apply that transformation to the rows of (avv) with respect to av2.
            j1 = j+1
            do 690 i=j1,4
               call fprota(co,si,h2(i),av2(ij,i))
 690        continue
 700     continue
c we update the sum of squared residuals
         do 705 i=1,nuu
           sq = sq+right(i)**2
 705     continue
         go to 750
c  rotation into triangle of the new row of (avv), in case the elements
c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
 710     irot =nrold
         do 740 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 740
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(irot,1),co,si)
c  apply that transformation to the columns of matrix g.
            ic = irot
            do 720 j=1,nuu
               call fprota(co,si,right(j),c(ic))
               ic = ic+nv7
 720        continue
c  apply that transformation to the rows of (avv).
            if(i.eq.5) go to 740
            i2 = 1
            i3 = i+1
            do 730 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),av1(irot,i2))
 730        continue
 740     continue
c we update the sum of squared residuals
         do 745 i=1,nuu
           sq = sq+right(i)**2
 745     continue
 750     if(nrold.eq.number) go to 770
 760     nrold = nrold+1
         go to 450
 770  continue
c  test whether the b-spline coefficients must be determined.
      if(iback.ne.0) return
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (rv) (cr) (ru)' = h.
c  first step: solve the system  (rv) (c1) = h.
      k = 1
      do 780 i=1,nuu
         call fpbacp(av1,av2,c(k),nv7,4,c(k),5,nv)
         k = k+nv7
 780  continue
c  second step: solve the system  (cr) (ru)' = (c1).
      k = 0
      do 800 j=1,nv7
        k = k+1
        l = k
        do 790 i=1,nuu
          right(i) = c(l)
          l = l+nv7
 790    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 795 i=1,nuu
           c(l) = right(i)
           l = l+nv7
 795    continue
 800  continue
c  calculate from the conditions (2)-(3)-(4), the remaining b-spline
c  coefficients.
      ncof = nu4*nv4
      i = nv4
      j = 0
      do 805 l=1,nv4
         q(l) = dz1
 805  continue
      if(iop0.eq.0) go to 815
      do 810 l=1,nv4
         i = i+1
         q(i) = cc(l)
 810  continue
 815  if(nuu.eq.0) go to 850
      do 840 l=1,nuu
         ii = i
         do 820 k=1,nv7
            i = i+1
            j = j+1
            q(i) = c(j)
 820     continue
         do 830 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 830     continue
 840  continue
 850  if(iop1.eq.0) go to 870
      do 860 l=1,nv4
         i = i+1
         q(i) = 0.
 860  continue
 870  do 880 i=1,ncof
         c(i) = q(i)
 880  continue
c  calculate the quantities
c    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
c    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
c    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
c                  tu(r+3) <= u(i) <= tu(r+4)
c    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
c                  tv(r+3) <= v(j) <= tv(r+4)
      fp = 0.
      do 890 i=1,nu
        fpu(i) = 0.
 890  continue
      do 900 i=1,nv
        fpv(i) = 0.
 900  continue
      iz = 0
      nroldu = 0
c  main loop for the different grid points.
      do 950 i1=1,mu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        do 940 i2=1,mv
          numv = nrv(i2)
          numv1 = numv+1
          iz = iz+1
c  evaluate s(u,v) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (u,v), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numu*nv4+numv
          do 920 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 910 l2=1,4
              k2 = k2+1
              term = term+fac*spv(i2,l2)*c(k2)
 910        continue
            k1 = k1+nv4
 920      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 930
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 930      nroldv = numv
          if(numu.eq.nroldu) go to 940
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 940    continue
        nroldu = numu
 950  continue
      return
      end
      subroutine fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,
     * tu,nu,tv,nv,p,c,nc,fp,fpu,fpv,mm,mvnu,spu,spv,right,q,au,au1,
     * av,av1,bu,bv,nru,nrv)
c  ..
c  ..scalar arguments..
      real p,fp
      integer ifsu,ifsv,ifbu,ifbv,idim,mu,mv,mz,nu,nv,nc,mm,mvnu
c  ..array arguments..
      real u(mu),v(mv),z(mz*idim),tu(nu),tv(nv),c(nc*idim),fpu(nu),
     * fpv(nv),spu(mu,4),spv(mv,4),right(mm*idim),q(mvnu),au(nu,5),
     * au1(nu,4),av(nv,5),av1(nv,4),bu(nu,5),bv(nv,5)
      integer ipar(2),nru(mu),nrv(mv)
c  ..local scalars..
      real arg,fac,term,one,half,value
      integer i,id,ii,it,iz,i1,i2,j,jz,k,k1,k2,l,l1,l2,mvv,k0,muu,
     * ncof,nroldu,nroldv,number,nmd,numu,numu1,numv,numv1,nuu,nvv,
     * nu4,nu7,nu8,nv4,nv7,nv8
c  ..local arrays..
      real h(5)
c  ..subroutine references..
c    fpback,fpbspl,fpdisc,fpbacp,fptrnp,fptrpe
c  ..
c  let
c               |   (spu)    |            |   (spv)    |
c        (au) = | ---------- |     (av) = | ---------- |
c               | (1/p) (bu) |            | (1/p) (bv) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
c                coefficients.
c       z      : the mu x mv matrix which contains the function values.
c       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
c                according to the least-squares problems in the u-,resp.
c                v-direction.
c       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
c                containing the discontinuity jumps of the derivatives
c                of the b-splines in the u-,resp.v-variable at the knots
c  the b-spline coefficients of the smoothing spline are then calculated
c  as the least-squares solution of the following over-determined linear
c  system of equations
c
c    (1)  (av) c (au)' = q
c
c  subject to the constraints
c
c    (2)  c(nu-3+i,j) = c(i,j), i=1,2,3 ; j=1,2,...,nv-4
c            if(ipar(1).ne.0)
c
c    (3)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
c            if(ipar(2).ne.0)
c
c  set constants
      one = 1
      half = 0.5
c  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      muu = mu
      if(ipar(1).ne.0) muu = mu-1
      mvv = mv
      if(ipar(2).ne.0) mvv = mv-1
c  it depends on the value of the flags ifsu,ifsv,ifbu  and ibvand
c  on the value of p whether the matrices (spu), (spv), (bu) and (bv)
c  still must be determined.
      if(ifsu.ne.0) go to 50
c  calculate the non-zero elements of the matrix (spu) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 40 it=1,muu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tu,nu,3,arg,l,h)
        do 30 i=1,4
          spu(it,i) = h(i)
  30    continue
        nru(it) = number
  40  continue
      ifsu = 1
c  calculate the non-zero elements of the matrix (spv) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the v-direction.
  50  if(ifsv.ne.0) go to 100
      l = 4
      l1 = 5
      number = 0
      do 90 it=1,mvv
        arg = v(it)
  60    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(tv,nv,3,arg,l,h)
        do 80 i=1,4
          spv(it,i) = h(i)
  80    continue
        nrv(it) = number
  90  continue
      ifsv = 1
 100  if(p.le.0.) go to  150
c  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 110
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
c  calculate the non-zero elements of the matrix (bv).
 110  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
c  substituting (2)  and (3) into (1), we obtain the overdetermined
c  system
c         (4)  (avv) (cr) (auu)' = (qq)
c  from which the nuu*nvv remaining coefficients
c         c(i,j) , i=1,...,nu-4-3*ipar(1) ; j=1,...,nv-4-3*ipar(2) ,
c  the elements of (cr), are then determined in the least-squares sense.
c  we first determine the matrices (auu) and (qq). then we reduce the
c  matrix (auu) to upper triangular form (ru) using givens rotations.
c  we apply the same transformations to the rows of matrix qq to obtain
c  the (mv) x nuu matrix g.
c  we store matrix (ru) into au (and au1 if ipar(1)=1) and g into q.
 150  if(ipar(1).ne.0) go to 160
      nuu = nu4
      call fptrnp(mu,mv,idim,nu,nru,spu,p,bu,z,au,q,right)
      go to 180
 160  nuu = nu7
      call fptrpe(mu,mv,idim,nu,nru,spu,p,bu,z,au,au1,q,right)
c  we determine the matrix (avv) and then we reduce this matrix to
c  upper triangular form (rv) using givens rotations.
c  we apply the same transformations to the columns of matrix
c  g to obtain the (nvv) x (nuu) matrix h.
c  we store matrix (rv) into av (and av1 if ipar(2)=1) and h into c.
 180  if(ipar(2).ne.0) go to 190
      nvv = nv4
      call fptrnp(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,c,right)
      go to 200
 190  nvv = nv7
      call fptrpe(mv,nuu,idim,nv,nrv,spv,p,bv,q,av,av1,c,right)
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (rv) (cr) (ru)' = h.
c  first step: solve the system  (rv) (c1) = h.
 200  ncof = nuu*nvv
      k = 1
      if(ipar(2).ne.0) go to 240
      do 220 ii=1,idim
      do 220 i=1,nuu
         call fpback(av,c(k),nvv,5,c(k),nv)
         k = k+nvv
 220  continue
      go to 300
 240  do 260 ii=1,idim
      do 260 i=1,nuu
         call fpbacp(av,av1,c(k),nvv,4,c(k),5,nv)
         k = k+nvv
 260  continue
c  second step: solve the system  (cr) (ru)' = (c1).
 300  if(ipar(1).ne.0) go to 400
      do 360 ii=1,idim
      k = (ii-1)*ncof
      do 360 j=1,nvv
        k = k+1
        l = k
        do 320 i=1,nuu
          right(i) = c(l)
          l = l+nvv
 320    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 340 i=1,nuu
           c(l) = right(i)
           l = l+nvv
 340    continue
 360  continue
      go to 500
 400  do 460 ii=1,idim
      k = (ii-1)*ncof
      do 460 j=1,nvv
        k = k+1
        l = k
        do 420 i=1,nuu
          right(i) = c(l)
          l = l+nvv
 420    continue
        call fpbacp(au,au1,right,nuu,4,right,5,nu)
        l = k
        do 440 i=1,nuu
           c(l) = right(i)
           l = l+nvv
 440    continue
 460  continue
c  calculate from the conditions (2)-(3), the remaining b-spline
c  coefficients.
 500  if(ipar(2).eq.0) go to 600
      i = 0
      j = 0
      do 560 id=1,idim
      do 560 l=1,nuu
         ii = i
         do 520 k=1,nvv
            i = i+1
            j = j+1
            q(i) = c(j)
 520     continue
         do 540 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 540     continue
 560  continue
      ncof = nv4*nuu
      nmd = ncof*idim
      do 580 i=1,nmd
         c(i) = q(i)
 580  continue
 600  if(ipar(1).eq.0) go to 700
      i = 0
      j = 0
      n33 = 3*nv4
      do 660 id=1,idim
         ii = i
         do 620 k=1,ncof
            i = i+1
            j = j+1
            q(i) = c(j)
 620     continue
         do 640 k=1,n33
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 640     continue
 660  continue
      ncof = nv4*nu4
      nmd = ncof*idim
      do 680 i=1,nmd
         c(i) = q(i)
 680  continue
c  calculate the quantities
c    res(i,j) = (z(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
c    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
c    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
c                  tu(r+3) <= u(i) <= tu(r+4)
c    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
c                  tv(r+3) <= v(j) <= tv(r+4)
 700  fp = 0.
      do 720 i=1,nu
        fpu(i) = 0.
 720  continue
      do 740 i=1,nv
        fpv(i) = 0.
 740  continue
      nroldu = 0
c  main loop for the different grid points.
      do 860 i1=1,muu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        iz = (i1-1)*mv
        do 840 i2=1,mvv
          numv = nrv(i2)
          numv1 = numv+1
          iz = iz+1
c  evaluate s(u,v) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (u,v), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k0 = numu*nv4+numv
          jz = iz
          do 800 id=1,idim
          k1 = k0
          value = 0.
          do 780 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 760 l2=1,4
              k2 = k2+1
              value = value+fac*spv(i2,l2)*c(k2)
 760        continue
            k1 = k1+nv4
 780      continue
c  calculate the squared residual at the current grid point.
          term = term+(z(jz)-value)**2
          jz = jz+mz
          k0 = k0+ncof
 800      continue
c  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 820
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 820      nroldv = numv
          if(numu.eq.nroldu) go to 840
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 840    continue
        nroldu = numu
 860  continue
      return
      end
      subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,
     * ty,ny,p,c,nc,fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,
     * ax,ay,bx,by,nrx,nry)
c  ..
c  ..scalar arguments..
      real p,fp
      integer ifsx,ifsy,ifbx,ifby,mx,my,mz,kx,ky,nx,ny,nc,mm,mynx,
     * kx1,kx2,ky1,ky2
c  ..array arguments..
      real x(mx),y(my),z(mz),tx(nx),ty(ny),c(nc),spx(mx,kx1),spy(my,ky1)
     * ,right(mm),q(mynx),ax(nx,kx2),bx(nx,kx2),ay(ny,ky2),by(ny,ky2),
     * fpx(nx),fpy(ny)
      integer nrx(mx),nry(my)
c  ..local scalars..
      real arg,cos,fac,pinv,piv,sin,term,one,half
      integer i,ibandx,ibandy,ic,iq,irot,it,iz,i1,i2,i3,j,k,k1,k2,l,
     * l1,l2,ncof,nk1x,nk1y,nrold,nroldx,nroldy,number,numx,numx1,
     * numy,numy1,n1
c  ..local arrays..
      real h(7)
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fprota
c  ..
c  the b-spline coefficients of the smoothing spline are calculated as
c  the least-squares solution of the over-determined linear system of
c  equations  (ay) c (ax)' = q       where
c
c               |   (spx)    |            |   (spy)    |
c        (ax) = | ---------- |     (ay) = | ---------- |
c               | (1/p) (bx) |            | (1/p) (by) |
c
c                                | z  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
c                b-spline coefficients.
c       z      : the my x mx matrix which contains the function values.
c       spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
c                matrices according to the least-squares problems in
c                the x- and y-direction.
c       bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
c                matrices which contain the discontinuity jumps of the
c                derivatives of the b-splines in the x- and y-direction.
      one = 1
      half = 0.5
      nk1x = nx-kx1
      nk1y = ny-ky1
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
c  the value of p whether the matrices (spx),(spy),(bx) and (by) still
c  must be determined.
      if(ifsx.ne.0) go to 50
c  calculate the non-zero elements of the matrix (spx) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the x-direction.
      l = kx1
      l1 = kx2
      number = 0
      do 40 it=1,mx
        arg = x(it)
  10    if(arg.lt.tx(l1) .or. l.eq.nk1x) go to 20
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        do 30 i=1,kx1
          spx(it,i) = h(i)
  30    continue
        nrx(it) = number
  40  continue
      ifsx = 1
  50  if(ifsy.ne.0) go to 100
c  calculate the non-zero elements of the matrix (spy) which is the
c  observation matrix according to the least-squares spline approximat-
c  ion problem in the y-direction.
      l = ky1
      l1 = ky2
      number = 0
      do 90 it=1,my
        arg = y(it)
  60    if(arg.lt.ty(l1) .or. l.eq.nk1y) go to 70
        l = l1
        l1 = l+1
        number = number+1
        go to 60
  70    call fpbspl(ty,ny,ky,arg,l,h)
        do 80 i=1,ky1
          spy(it,i) = h(i)
  80    continue
        nry(it) = number
  90  continue
      ifsy = 1
 100  if(p.le.0.) go to 120
c  calculate the non-zero elements of the matrix (bx).
      if(ifbx.ne.0 .or. nx.eq.2*kx1) go to 110
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
c  calculate the non-zero elements of the matrix (by).
 110  if(ifby.ne.0 .or. ny.eq.2*ky1) go to 120
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
c  reduce the matrix (ax) to upper triangular form (rx) using givens
c  rotations. apply the same transformations to the rows of matrix q
c  to obtain the my x (nx-kx-1) matrix g.
c  store matrix (rx) into (ax) and g into q.
 120  l = my*nk1x
c  initialization.
      do 130 i=1,l
        q(i) = 0.
 130  continue
      do 140 i=1,nk1x
        do 140 j=1,kx2
          ax(i,j) = 0.
 140  continue
      l = 0
      nrold = 0
c  ibandx denotes the bandwidth of the matrices (ax) and (rx).
      ibandx = kx1
      do 270 it=1,mx
        number = nrx(it)
 150    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 260
        ibandx = kx2
c  fetch a new row of matrix (bx).
        n1 = nrold+1
        do 160 j=1,kx2
          h(j) = bx(n1,j)*pinv
 160    continue
c  find the appropriate column of q.
        do 170 j=1,my
          right(j) = 0.
 170    continue
        irot = nrold
        go to 210
c  fetch a new row of matrix (spx).
 180    h(ibandx) = 0.
        do 190 j=1,kx1
          h(j) = spx(it,j)
 190    continue
c  find the appropriate column of q.
        do 200 j=1,my
          l = l+1
          right(j) = z(l)
 200    continue
        irot = number
c  rotate the new row of matrix (ax) into triangle.
 210    do 240 i=1,ibandx
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 240
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ax(irot,1),cos,sin)
c  apply that transformation to the rows of matrix q.
          iq = (irot-1)*my
          do 220 j=1,my
            iq = iq+1
            call fprota(cos,sin,right(j),q(iq))
 220      continue
c  apply that transformation to the columns of (ax).
          if(i.eq.ibandx) go to 250
          i2 = 1
          i3 = i+1
          do 230 j=i3,ibandx
            i2 = i2+1
            call fprota(cos,sin,h(j),ax(irot,i2))
 230      continue
 240    continue
 250    if(nrold.eq.number) go to 270
 260    nrold = nrold+1
        go to 150
 270  continue
c  reduce the matrix (ay) to upper triangular form (ry) using givens
c  rotations. apply the same transformations to the columns of matrix g
c  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
c  store matrix (ry) into (ay) and h into c.
      ncof = nk1x*nk1y
c  initialization.
      do 280 i=1,ncof
        c(i) = 0.
 280  continue
      do 290 i=1,nk1y
        do 290 j=1,ky2
          ay(i,j) = 0.
 290  continue
      nrold = 0
c  ibandy denotes the bandwidth of the matrices (ay) and (ry).
      ibandy = ky1
      do 420 it=1,my
        number = nry(it)
 300    if(nrold.eq.number) go to 330
        if(p.le.0.) go to 410
        ibandy = ky2
c  fetch a new row of matrix (by).
        n1 = nrold+1
        do 310 j=1,ky2
          h(j) = by(n1,j)*pinv
 310    continue
c  find the appropiate row of g.
        do 320 j=1,nk1x
          right(j) = 0.
 320    continue
        irot = nrold
        go to 360
c  fetch a new row of matrix (spy)
 330    h(ibandy) = 0.
        do 340 j=1,ky1
          h(j) = spy(it,j)
 340    continue
c  find the appropiate row of g.
        l = it
        do 350 j=1,nk1x
          right(j) = q(l)
          l = l+my
 350    continue
        irot = number
c  rotate the new row of matrix (ay) into triangle.
 360    do 390 i=1,ibandy
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 390
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,ay(irot,1),cos,sin)
c  apply that transformation to the colums of matrix g.
          ic = irot
          do 370 j=1,nk1x
            call fprota(cos,sin,right(j),c(ic))
            ic = ic+nk1y
 370      continue
c  apply that transformation to the columns of matrix (ay).
          if(i.eq.ibandy) go to 400
          i2 = 1
          i3 = i+1
          do 380 j=i3,ibandy
            i2 = i2+1
            call fprota(cos,sin,h(j),ay(irot,i2))
 380      continue
 390    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = nrold+1
        go to 300
 420  continue
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (ry) c (rx)' = h.
c  first step: solve the system  (ry) (c1) = h.
      k = 1
      do 450 i=1,nk1x
        call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
        k = k+nk1y
 450  continue
c  second step: solve the system  c (rx)' = (c1).
      k = 0
      do 480 j=1,nk1y
        k = k+1
        l = k
        do 460 i=1,nk1x
          right(i) = c(l)
          l = l+nk1y
 460    continue
        call fpback(ax,right,nk1x,ibandx,right,nx)
        l = k
        do 470 i=1,nk1x
          c(l) = right(i)
          l = l+nk1y
 470    continue
 480  continue
c  calculate the quantities
c    res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
c    fp = sumi=1,mx(sumj=1,my(res(i,j)))
c    fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
c                  tx(r+kx) <= x(i) <= tx(r+kx+1)
c    fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
c                  ty(r+ky) <= y(j) <= ty(r+ky+1)
      fp = 0.
      do 490 i=1,nx
        fpx(i) = 0.
 490  continue
      do 500 i=1,ny
        fpy(i) = 0.
 500  continue
      nk1y = ny-ky1
      iz = 0
      nroldx = 0
c  main loop for the different grid points.
      do 550 i1=1,mx
        numx = nrx(i1)
        numx1 = numx+1
        nroldy = 0
        do 540 i2=1,my
          numy = nry(i2)
          numy1 = numy+1
          iz = iz+1
c  evaluate s(x,y) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (x,y), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numx*nk1y+numy
          do 520 l1=1,kx1
            k2 = k1
            fac = spx(i1,l1)
            do 510 l2=1,ky1
              k2 = k2+1
              term = term+fac*spy(i2,l2)*c(k2)
 510        continue
            k1 = k1+nk1y
 520      continue
c  calculate the squared residual at the current grid point.
          term = (z(iz)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpx(numx1) = fpx(numx1)+term
          fpy(numy1) = fpy(numy1)+term
          fac = term*half
          if(numy.eq.nroldy) go to 530
          fpy(numy1) = fpy(numy1)-fac
          fpy(numy) = fpy(numy)+fac
 530      nroldy = numy
          if(numx.eq.nroldx) go to 540
          fpx(numx1) = fpx(numx1)-fac
          fpx(numx) = fpx(numx)+fac
 540    continue
        nroldx = numx
 550  continue
      return
      end

      subroutine fpgrsp(ifsu,ifsv,ifbu,ifbv,iback,u,mu,v,mv,r,mr,dr,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,spu,spv,
     * right,q,au,av1,av2,bu,bv,a0,a1,b0,b1,c0,c1,cosi,nru,nrv)
c  ..
c  ..scalar arguments..
      real p,sq,fp
      integer ifsu,ifsv,ifbu,ifbv,iback,mu,mv,mr,iop0,iop1,nu,nv,nc,
     * mm,mvnu
c  ..array arguments..
      real u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),
     * spu(mu,4),spv(mv,4),right(mm),q(mvnu),au(nu,5),av1(nv,6),c0(nv),
     * av2(nv,4),a0(2,mv),b0(2,nv),cosi(2,nv),bu(nu,5),bv(nv,5),c1(nv),
     * a1(2,mv),b1(2,nv)
      integer nru(mu),nrv(mv)
c  ..local scalars..
      real arg,co,dr01,dr02,dr03,dr11,dr12,dr13,fac,fac0,fac1,pinv,piv,
     * si,term,one,three,half
      integer i,ic,ii,ij,ik,iq,irot,it,ir,i0,i1,i2,i3,j,jj,jk,jper,
     * j0,j1,k,k1,k2,l,l0,l1,l2,mvv,ncof,nrold,nroldu,nroldv,number,
     * numu,numu1,numv,numv1,nuu,nu4,nu7,nu8,nu9,nv11,nv4,nv7,nv8,n1
c  ..local arrays..
      real h(5),h1(5),h2(4)
c  ..function references..
      integer min0
      real cos,sin
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpcyt1,fpcyt2,fpdisc,fpbacp,fprota
c  ..
c  let
c               |     (spu)      |            |     (spv)      |
c        (au) = | -------------- |     (av) = | -------------- |
c               | sqrt(1/p) (bu) |            | sqrt(1/p) (bv) |
c
c                                | r  ' 0 |
c                            q = | ------ |
c                                | 0  ' 0 |
c
c  with c      : the (nu-4) x (nv-4) matrix which contains the b-spline
c                coefficients.
c       r      : the mu x mv matrix which contains the function values.
c       spu,spv: the mu x (nu-4), resp. mv x (nv-4) observation matrices
c                according to the least-squares problems in the u-,resp.
c                v-direction.
c       bu,bv  : the (nu-7) x (nu-4),resp. (nv-7) x (nv-4) matrices
c                containing the discontinuity jumps of the derivatives
c                of the b-splines in the u-,resp.v-variable at the knots
c  the b-spline coefficients of the smoothing spline are then calculated
c  as the least-squares solution of the following over-determined linear
c  system of equations
c
c  (1)  (av) c (au)' = q
c
c  subject to the constraints
c
c  (2)  c(i,nv-3+j) = c(i,j), j=1,2,3 ; i=1,2,...,nu-4
c
c  (3)  if iop0 = 0  c(1,j) = dr(1)
c          iop0 = 1  c(1,j) = dr(1)
c                    c(2,j) = dr(1)+(dr(2)*cosi(1,j)+dr(3)*cosi(2,j))*
c                            tu(5)/3. = c0(j) , j=1,2,...nv-4
c
c  (4)  if iop1 = 0  c(nu-4,j) = dr(4)
c          iop1 = 1  c(nu-4,j) = dr(4)
c                    c(nu-5,j) = dr(4)+(dr(5)*cosi(1,j)+dr(6)*cosi(2,j))
c                                *(tu(nu-4)-tu(nu-3))/3. = c1(j)
c
c  set constants
      one = 1
      three = 3
      half = 0.5
c  initialization
      nu4 = nu-4
      nu7 = nu-7
      nu8 = nu-8
      nu9 = nu-9
      nv4 = nv-4
      nv7 = nv-7
      nv8 = nv-8
      nv11 = nv-11
      nuu = nu4-iop0-iop1-2
      if(p.gt.0.) pinv = one/p
c  it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1
c  and on the value of p whether the matrices (spu), (spv), (bu), (bv),
c  (cosi) still must be determined.
      if(ifsu.ne.0) go to 30
c  calculate the non-zero elements of the matrix (spu) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the u-direction.
      l = 4
      l1 = 5
      number = 0
      do 25 it=1,mu
        arg = u(it)
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 15
        l = l1
        l1 = l+1
        number = number+1
        go to 10
  15    call fpbspl(tu,nu,3,arg,l,h)
        do 20 i=1,4
          spu(it,i) = h(i)
  20    continue
        nru(it) = number
  25  continue
      ifsu = 1
c  calculate the non-zero elements of the matrix (spv) which is the ob-
c  servation matrix according to the least-squares spline approximation
c  problem in the v-direction.
  30  if(ifsv.ne.0) go to 85
      l = 4
      l1 = 5
      number = 0
      do 50 it=1,mv
        arg = v(it)
  35    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 40
        l = l1
        l1 = l+1
        number = number+1
        go to 35
  40    call fpbspl(tv,nv,3,arg,l,h)
        do 45 i=1,4
          spv(it,i) = h(i)
  45    continue
        nrv(it) = number
  50  continue
      ifsv = 1
      if(iop0.eq.0 .and. iop1.eq.0) go to 85
c  calculate the coefficients of the interpolating splines for cos(v)
c  and sin(v).
      do 55 i=1,nv4
         cosi(1,i) = 0.
         cosi(2,i) = 0.
  55  continue
      if(nv7.lt.4) go to 85
      do 65 i=1,nv7
         l = i+3
         arg = tv(l)
         call fpbspl(tv,nv,3,arg,l,h)
         do 60 j=1,3
            av1(i,j) = h(j)
  60     continue
         cosi(1,i) = cos(arg)
         cosi(2,i) = sin(arg)
  65  continue
      call fpcyt1(av1,nv7,nv)
      do 80 j=1,2
         do 70 i=1,nv7
            right(i) = cosi(j,i)
  70     continue
         call fpcyt2(av1,nv7,right,right,nv)
         do 75 i=1,nv7
            cosi(j,i+1) = right(i)
  75     continue
         cosi(j,1) = cosi(j,nv7+1)
         cosi(j,nv7+2) = cosi(j,2)
         cosi(j,nv4) = cosi(j,3)
  80  continue
  85  if(p.le.0.) go to  150
c  calculate the non-zero elements of the matrix (bu).
      if(ifbu.ne.0 .or. nu8.eq.0) go to 90
      call fpdisc(tu,nu,5,bu,nu)
      ifbu = 1
c  calculate the non-zero elements of the matrix (bv).
  90  if(ifbv.ne.0 .or. nv8.eq.0) go to 150
      call fpdisc(tv,nv,5,bv,nv)
      ifbv = 1
c  substituting (2),(3) and (4) into (1), we obtain the overdetermined
c  system
c         (5)  (avv) (cc) (auu)' = (qq)
c  from which the nuu*nv7 remaining coefficients
c         c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
c  the elements of (cc), are then determined in the least-squares sense.
c  simultaneously, we compute the resulting sum of squared residuals sq.
 150  dr01 = dr(1)
      dr11 = dr(4)
      do 155 i=1,mv
         a0(1,i) = dr01
         a1(1,i) = dr11
 155  continue
      if(nv8.eq.0 .or. p.le.0.) go to 165
      do 160 i=1,nv8
         b0(1,i) = 0.
         b1(1,i) = 0.
 160  continue
 165  mvv = mv
      if(iop0.eq.0) go to 195
      fac = (tu(5)-tu(4))/three
      dr02 = dr(2)*fac
      dr03 = dr(3)*fac
      do 170 i=1,nv4
         c0(i) = dr01+dr02*cosi(1,i)+dr03*cosi(2,i)
 170  continue
      do 180 i=1,mv
         number = nrv(i)
         fac = 0.
         do 175 j=1,4
            number = number+1
            fac = fac+c0(number)*spv(i,j)
 175     continue
         a0(2,i) = fac
 180  continue
      if(nv8.eq.0 .or. p.le.0.) go to 195
      do 190 i=1,nv8
         number = i
         fac = 0.
         do 185 j=1,5
            fac = fac+c0(number)*bv(i,j)
            number = number+1
 185     continue
         b0(2,i) = fac*pinv
 190  continue
      mvv = mv+nv8
 195  if(iop1.eq.0) go to 225
      fac = (tu(nu4)-tu(nu4+1))/three
      dr12 = dr(5)*fac
      dr13 = dr(6)*fac
      do 200 i=1,nv4
         c1(i) = dr11+dr12*cosi(1,i)+dr13*cosi(2,i)
 200  continue
      do 210 i=1,mv
         number = nrv(i)
         fac = 0.
         do 205 j=1,4
            number = number+1
            fac = fac+c1(number)*spv(i,j)
 205     continue
         a1(2,i) = fac
 210  continue
      if(nv8.eq.0 .or. p.le.0.) go to 225
      do 220 i=1,nv8
         number = i
         fac = 0.
         do 215 j=1,5
            fac = fac+c1(number)*bv(i,j)
            number = number+1
 215     continue
         b1(2,i) = fac*pinv
 220  continue
      mvv = mv+nv8
c  we first determine the matrices (auu) and (qq). then we reduce the
c  matrix (auu) to an unit upper triangular form (ru) using givens
c  rotations without square roots. we apply the same transformations to
c  the rows of matrix qq to obtain the mv x nuu matrix g.
c  we store matrix (ru) into au and g into q.
 225  l = mvv*nuu
c  initialization.
      sq = 0.
      if(l.eq.0) go to 245
      do 230 i=1,l
        q(i) = 0.
 230  continue
      do 240 i=1,nuu
        do 240 j=1,5
          au(i,j) = 0.
 240  continue
      l = 0
 245  nrold = 0
      n1 = nrold+1
      do 420 it=1,mu
        number = nru(it)
c  find the appropriate column of q.
 250    do 260 j=1,mvv
           right(j) = 0.
 260    continue
        if(nrold.eq.number) go to 280
        if(p.le.0.) go to 410
c  fetch a new row of matrix (bu).
        do 270 j=1,5
          h(j) = bu(n1,j)*pinv
 270    continue
        i0 = 1
        i1 = 5
        go to 310
c  fetch a new row of matrix (spu).
 280    do 290 j=1,4
          h(j) = spu(it,j)
 290    continue
c  find the appropriate column of q.
        do 300 j=1,mv
          l = l+1
          right(j) = r(l)
 300    continue
        i0 = 1
        i1 = 4
 310    j0 = n1
        j1 = nu7-number
c  take into account that we eliminate the constraints (3)
 315     if(j0-1.gt.iop0) go to 335
         fac0 = h(i0)
         do 320 j=1,mv
            right(j) = right(j)-fac0*a0(j0,j)
 320     continue
         if(mv.eq.mvv) go to 330
         j = mv
         do 325 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac0*b0(j0,jj)
 325     continue
 330     j0 = j0+1
         i0 = i0+1
         go to 315
c  take into account that we eliminate the constraints (4)
 335     if(j1-1.gt.iop1) go to 360
         fac1 = h(i1)
         do 340 j=1,mv
            right(j) = right(j)-fac1*a1(j1,j)
 340     continue
         if(mv.eq.mvv) go to 350
         j = mv
         do 345 jj=1,nv8
            j = j+1
            right(j) = right(j)-fac1*b1(j1,jj)
 345     continue
 350     j1 = j1+1
         i1 = i1-1
         go to 335
 360     irot = nrold-iop0-1
         if(irot.lt.0) irot = 0
c  rotate the new row of matrix (auu) into triangle.
        if(i0.gt.i1) go to 390
        do 385 i=i0,i1
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 385
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,au(irot,1),co,si)
c  apply that transformation to the rows of matrix (qq).
          iq = (irot-1)*mvv
          do 370 j=1,mvv
            iq = iq+1
            call fprota(co,si,right(j),q(iq))
 370      continue
c  apply that transformation to the columns of (auu).
          if(i.eq.i1) go to 385
          i2 = 1
          i3 = i+1
          do 380 j=i3,i1
            i2 = i2+1
            call fprota(co,si,h(j),au(irot,i2))
 380      continue
 385    continue
c  we update the sum of squared residuals.
 390    do 395 j=1,mvv
          sq = sq+right(j)**2
 395    continue
 400    if(nrold.eq.number) go to 420
 410    nrold = n1
        n1 = n1+1
        go to 250
 420  continue
      if(nuu.eq.0) go to 800
c  we determine the matrix (avv) and then we reduce her to an unit
c  upper triangular form (rv) using givens rotations without square
c  roots. we apply the same transformations to the columns of matrix
c  g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h.
c  we store matrix (rv) into av1 and av2, h into c.
c  the nv7 x nv7 triangular unit upper matrix (rv) has the form
c              | av1 '     |
c       (rv) = |     ' av2 |
c              |  0  '     |
c  with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper
c  triangular matrix of bandwidth 5.
      ncof = nuu*nv7
c  initialization.
      do 430 i=1,ncof
        c(i) = 0.
 430  continue
      do 440 i=1,nv4
        av1(i,5) = 0.
        do 440 j=1,4
          av1(i,j) = 0.
          av2(i,j) = 0.
 440  continue
      jper = 0
      nrold = 0
      do 770 it=1,mv
        number = nrv(it)
 450    if(nrold.eq.number) go to 480
        if(p.le.0.) go to 760
c  fetch a new row of matrix (bv).
        n1 = nrold+1
        do 460 j=1,5
          h(j) = bv(n1,j)*pinv
 460    continue
c  find the appropiate row of g.
        do 465 j=1,nuu
          right(j) = 0.
 465    continue
        if(mv.eq.mvv) go to 510
        l = mv+n1
        do 470 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 470    continue
        go to 510
c  fetch a new row of matrix (spv)
 480    h(5) = 0.
        do 490 j=1,4
          h(j) = spv(it,j)
 490    continue
c  find the appropiate row of g.
        l = it
        do 500 j=1,nuu
          right(j) = q(l)
          l = l+mvv
 500    continue
c  test whether there are non-zero values in the new row of (avv)
c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
 510     if(nrold.lt.nv11) go to 710
         if(jper.ne.0) go to 550
c  initialize the matrix (av2).
         jk = nv11+1
         do 540 i=1,4
            ik = jk
            do 520 j=1,5
               if(ik.le.0) go to 530
               av2(ik,i) = av1(ik,j)
               ik = ik-1
 520        continue
 530        jk = jk+1
 540     continue
         jper = 1
c  if one of the non-zero elements of the new row corresponds to one of
c  the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
c  (2) for setting up this row of (avv). the row is stored in h1( the
c  part with respect to av1) and h2 (the part with respect to av2).
 550     do 560 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 560     continue
         h1(5) = 0.
         j = nrold-nv11
         do 600 i=1,5
            j = j+1
            l0 = j
 570        l1 = l0-4
            if(l1.le.0) go to 590
            if(l1.le.nv11) go to 580
            l0 = l1-nv11
            go to 570
 580        h1(l1) = h(i)
            go to 600
 590        h2(l0) = h2(l0) + h(i)
 600     continue
c  rotate the new row of (avv) into triangle.
         if(nv11.le.0) go to 670
c  rotations with the rows 1,2,...,nv11 of (avv).
         do 660 j=1,nv11
            piv = h1(1)
            i2 = min0(nv11-j,4)
            if(piv.eq.0.) go to 640
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(j,1),co,si)
c  apply that transformation to the columns of matrix g.
            ic = j
            do 610 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 610        continue
c  apply that transformation to the rows of (avv) with respect to av2.
            do 620 i=1,4
               call fprota(co,si,h2(i),av2(j,i))
 620        continue
c  apply that transformation to the rows of (avv) with respect to av1.
            if(i2.eq.0) go to 670
            do 630 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),av1(j,i1))
 630        continue
 640        do 650 i=1,i2
               h1(i) = h1(i+1)
 650        continue
            h1(i2+1) = 0.
 660     continue
c  rotations with the rows nv11+1,...,nv7 of avv.
 670     do 700 j=1,4
            ij = nv11+j
            if(ij.le.0) go to 700
            piv = h2(j)
            if(piv.eq.0.) go to 700
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av2(ij,j),co,si)
c  apply that transformation to the columns of matrix g.
            ic = ij
            do 680 i=1,nuu
               call fprota(co,si,right(i),c(ic))
               ic = ic+nv7
 680        continue
            if(j.eq.4) go to 700
c  apply that transformation to the rows of (avv) with respect to av2.
            j1 = j+1
            do 690 i=j1,4
               call fprota(co,si,h2(i),av2(ij,i))
 690        continue
 700     continue
c  we update the sum of squared residuals.
         do 705 i=1,nuu
           sq = sq+right(i)**2
 705     continue
         go to 750
c  rotation into triangle of the new row of (avv), in case the elements
c  corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
 710     irot =nrold
         do 740 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 740
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,av1(irot,1),co,si)
c  apply that transformation to the columns of matrix g.
            ic = irot
            do 720 j=1,nuu
               call fprota(co,si,right(j),c(ic))
               ic = ic+nv7
 720        continue
c  apply that transformation to the rows of (avv).
            if(i.eq.5) go to 740
            i2 = 1
            i3 = i+1
            do 730 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),av1(irot,i2))
 730        continue
 740     continue
c  we update the sum of squared residuals.
         do 745 i=1,nuu
           sq = sq+right(i)**2
 745     continue
 750     if(nrold.eq.number) go to 770
 760     nrold = nrold+1
         go to 450
 770  continue
c  test whether the b-spline coefficients must be determined.
      if(iback.ne.0) return
c  backward substitution to obtain the b-spline coefficients as the
c  solution of the linear system    (rv) (cr) (ru)' = h.
c  first step: solve the system  (rv) (c1) = h.
      k = 1
      do 780 i=1,nuu
         call fpbacp(av1,av2,c(k),nv7,4,c(k),5,nv)
         k = k+nv7
 780  continue
c  second step: solve the system  (cr) (ru)' = (c1).
      k = 0
      do 795 j=1,nv7
        k = k+1
        l = k
        do 785 i=1,nuu
          right(i) = c(l)
          l = l+nv7
 785    continue
        call fpback(au,right,nuu,5,right,nu)
        l = k
        do 790 i=1,nuu
           c(l) = right(i)
           l = l+nv7
 790    continue
 795  continue
c  calculate from the conditions (2)-(3)-(4), the remaining b-spline
c  coefficients.
 800  ncof = nu4*nv4
      j = ncof
      do 805 l=1,nv4
         q(l) = dr01
         q(j) = dr11
         j = j-1
 805  continue
      i = nv4
      j = 0
      if(iop0.eq.0) go to 815
      do 810 l=1,nv4
         i = i+1
         q(i) = c0(l)
 810  continue
 815  if(nuu.eq.0) go to 835
      do 830 l=1,nuu
         ii = i
         do 820 k=1,nv7
            i = i+1
            j = j+1
            q(i) = c(j)
 820     continue
         do 825 k=1,3
            ii = ii+1
            i = i+1
            q(i) = q(ii)
 825     continue
 830  continue
 835  if(iop1.eq.0) go to 845
      do 840 l=1,nv4
         i = i+1
         q(i) = c1(l)
 840  continue
 845  do 850 i=1,ncof
         c(i) = q(i)
 850  continue
c  calculate the quantities
c    res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
c    fp = sumi=1,mu(sumj=1,mv(res(i,j)))
c    fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
c                  tu(r+3) <= u(i) <= tu(r+4)
c    fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
c                  tv(r+3) <= v(j) <= tv(r+4)
      fp = 0.
      do 890 i=1,nu
        fpu(i) = 0.
 890  continue
      do 900 i=1,nv
        fpv(i) = 0.
 900  continue
      ir = 0
      nroldu = 0
c  main loop for the different grid points.
      do 950 i1=1,mu
        numu = nru(i1)
        numu1 = numu+1
        nroldv = 0
        do 940 i2=1,mv
          numv = nrv(i2)
          numv1 = numv+1
          ir = ir+1
c  evaluate s(u,v) at the current grid point by making the sum of the
c  cross products of the non-zero b-splines at (u,v), multiplied with
c  the appropiate b-spline coefficients.
          term = 0.
          k1 = numu*nv4+numv
          do 920 l1=1,4
            k2 = k1
            fac = spu(i1,l1)
            do 910 l2=1,4
              k2 = k2+1
              term = term+fac*spv(i2,l2)*c(k2)
 910        continue
            k1 = k1+nv4
 920      continue
c  calculate the squared residual at the current grid point.
          term = (r(ir)-term)**2
c  adjust the different parameters.
          fp = fp+term
          fpu(numu1) = fpu(numu1)+term
          fpv(numv1) = fpv(numv1)+term
          fac = term*half
          if(numv.eq.nroldv) go to 930
          fpv(numv1) = fpv(numv1)-fac
          fpv(numv) = fpv(numv)+fac
 930      nroldv = numv
          if(numu.eq.nroldu) go to 940
          fpu(numu1) = fpu(numu1)-fac
          fpu(numu) = fpu(numu)+fac
 940    continue
        nroldu = numu
 950  continue
      return
      end
      subroutine fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)
c  given the b-spline representation (knots t(j),j=1,2,...,n, b-spline
c  coefficients c(j),j=1,2,...,n-k-1) of a spline of degree k, fpinst
c  calculates the b-spline representation (knots tt(j),j=1,2,...,nn,
c  b-spline coefficients cc(j),j=1,2,...,nn-k-1) of the same spline if
c  an additional knot is inserted at the point x situated in the inter-
c  val t(l)<=x<t(l+1). iopt denotes whether (iopt.ne.0) or not (iopt=0)
c  the given spline is periodic. in case of a periodic spline at least
c  one of the following conditions must be fulfilled: l>2*k or l<n-2*k.
c
c  ..scalar arguments..
      integer k,n,l,nn,iopt,nest
      real x
c  ..array arguments..
      real t(nest),c(nest),tt(nest),cc(nest)
c  ..local scalars..
      real fac,per,one
      integer i,i1,j,k1,m,mk,nk,nk1,nl,ll
c  ..
      one = 0.1e+01
      k1 = k+1
      nk1 = n-k1
c  the new knots
      ll = l+1
      i = n
      do 10 j=ll,n
         tt(i+1) = t(i)
         i = i-1
  10  continue
      tt(ll) = x
      do 20 j=1,l
         tt(j) = t(j)
  20  continue
c  the new b-spline coefficients
      i = nk1
      do 30 j=l,nk1
         cc(i+1) = c(i)
         i = i-1
  30  continue
      i = l
      do 40 j=1,k
         m = i+k1
         fac = (x-tt(i))/(tt(m)-tt(i))
         i1 = i-1
         cc(i) = fac*c(i)+(one-fac)*c(i1)
         i = i1
  40  continue
      do 50 j=1,i
         cc(j) = c(j)
  50  continue
      nn = n+1
      if(iopt.eq.0) return
c   incorporate the boundary conditions for a periodic spline.
      nk = nn-k
      nl = nk-k1
      per = tt(nk)-tt(k1)
      i = k1
      j = nk
      if(ll.le.nl) go to 70
      do 60 m=1,k
         mk = m+nl
         cc(m) = cc(mk)
         i = i-1
         j = j-1
         tt(i) = tt(j)-per
  60  continue
      return
  70  if(ll.gt.(k1+k)) return
      do 80 m=1,k
         mk = m+nl
         cc(mk) = cc(m)
         i = i+1
         j = j+1
         tt(j) = tt(i)+per
  80  continue
      return
      end
      subroutine fpintb(t,n,bint,nk1,x,y)
c  subroutine fpintb calculates integrals of the normalized b-splines
c  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
c  it makes use of the formulae of gaffney for the calculation of
c  indefinite integrals of b-splines.
c
c  calling sequence:
c     call fpintb(t,n,bint,nk1,x,y)
c
c  input parameters:
c    t    : real array,length n, containing the position of the knots.
c    n    : integer value, giving the number of knots.
c    nk1  : integer value, giving the number of b-splines of degree k,
c           defined on the set of knots ,i.e. nk1 = n-k-1.
c    x,y  : real values, containing the end points of the integration
c           interval.
c  output parameter:
c    bint : array,length nk1, containing the integrals of the b-splines.
c  ..
c  ..scalars arguments..
      integer n,nk1
      real x,y
c  ..array arguments..
      real t(n),bint(nk1)
c  ..local scalars..
      integer i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0,min
      real a,ak,arg,b,f,one
c  ..local arrays..
      real aint(6),h(6),h1(6)
c  initialization.
      one = 0.1e+01
      k1 = n-nk1
      ak = k1
      k = k1-1
      do 10 i=1,nk1
        bint(i) = 0.
  10  continue
c  the integration limits are arranged in increasing order.
      a = x
      b = y
      min = 0
      if(a-b) 30,160,20
  20  a = y
      b = x
      min = 1
  30  if(a.lt.t(k1)) a = t(k1)
      if(b.gt.t(nk1+1)) b = t(nk1+1)
c  using the expression of gaffney for the indefinite integral of a
c  b-spline we find that
c  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
c    where for t(l) <= x < t(l+1)
c    res(j,x) = 0, j=1,2,...,l-k-1
c             = 1, j=l+1,l+2,...,nk1
c             = aint(j+k-l+1), j=l-k,l-k+1,...,l
c               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
c                 i=0,1,...,k
      l = k1
      l0 = l+1
c  set arg = a.
      arg = a
      do 90 it=1,2
c  search for the knot interval t(l) <= arg < t(l+1).
  40    if(arg.lt.t(l0) .or. l.eq.nk1) go to 50
        l = l0
        l0 = l+1
        go to 40
c  calculation of aint(j), j=1,2,...,k+1.
c  initialization.
  50    do 55 j=1,k1
          aint(j) = 0.
  55    continue
        aint(1) = (arg-t(l))/(t(l+1)-t(l))
        h1(1) = one
        do 70 j=1,k
c  evaluation of the non-zero b-splines of degree j at arg,i.e.
c    h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
          h(1) = 0.
          do 60 i=1,j
            li = l+i
            lj = li-j
            f = h1(i)/(t(li)-t(lj))
            h(i) = h(i)+f*(t(li)-arg)
            h(i+1) = f*(arg-t(lj))
  60      continue
c  updating of the integrals aint.
          j1 = j+1
          do 70 i=1,j1
            li = l+i
            lj = li-j1
            aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
            h1(i) = h(i)
  70    continue
        if(it.eq.2) go to 100
c  updating of the integrals bint
        lk = l-k
        ia = lk
        do 80 i=1,k1
          bint(lk) = -aint(i)
          lk = lk+1
  80    continue
c  set arg = b.
        arg = b
  90  continue
c  updating of the integrals bint.
 100  lk = l-k
      ib = lk-1
      do 110 i=1,k1
        bint(lk) = bint(lk)+aint(i)
        lk = lk+1
 110  continue
      if(ib.lt.ia) go to 130
      do 120 i=ia,ib
        bint(i) = bint(i)+one
 120  continue
c  the scaling factors are taken into account.
 130  f = one/ak
      do 140 i=1,nk1
        j = i+k1
        bint(i) = bint(i)*(t(j)-t(i))*f
 140  continue
c  the order of the integration limits is taken into account.
      if(min.eq.0) go to 160
      do 150 i=1,nk1
        bint(i) = -bint(i)
 150  continue
 160  return
      end
      subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
c  subroutine fpknot locates an additional knot for a spline of degree
c  k and adjusts the corresponding parameters,i.e.
c    t     : the position of the knots.
c    n     : the number of knots.
c    nrint : the number of knotintervals.
c    fpint : the sum of squares of residual right hand sides
c            for each knot interval.
c    nrdata: the number of data points inside each knot interval.
c  istart indicates that the smallest data point at which the new knot
c  may be added is x(istart+1)
c  ..
c  ..scalar arguments..
      integer m,n,nrint,nest,istart
c  ..array arguments..
      real x(m),t(nest),fpint(nest)
      integer nrdata(nest)
c  ..local scalars..
      real an,am,fpmax
      integer ihalf,j,jbegin,jj,jk,jpoint,k,maxbeg,maxpt,
     * next,nrx,number
c  ..
      k = (n-nrint-1)/2
c  search for knot interval t(number+k) <= x <= t(number+k+1) where
c  fpint(number) is maximal on the condition that nrdata(number)
c  not equals zero.
      fpmax = 0.
      jbegin = istart
      do 20 j=1,nrint
        jpoint = nrdata(j)
        if(fpmax.ge.fpint(j) .or. jpoint.eq.0) go to 10
        fpmax = fpint(j)
        number = j
        maxpt = jpoint
        maxbeg = jbegin
  10    jbegin = jbegin+jpoint+1
  20  continue
c  let coincide the new knot t(number+k+1) with a data point x(nrx)
c  inside the old knot interval t(number+k) <= x <= t(number+k+1).
      ihalf = maxpt/2+1
      nrx = maxbeg+ihalf
      next = number+1
      if(next.gt.nrint) go to 40
c  adjust the different parameters.
      do 30 j=next,nrint
        jj = next+nrint-j
        fpint(jj+1) = fpint(jj)
        nrdata(jj+1) = nrdata(jj)
        jk = jj+k
        t(jk+1) = t(jk)
  30  continue
  40  nrdata(number) = ihalf-1
      nrdata(next) = maxpt-ihalf
      am = maxpt
      an = nrdata(number)
      fpint(number) = fpmax*an/am
      an = nrdata(next)
      fpint(next) = fpmax*an/am
      jk = next+k
      t(jk) = x(nrx)
      n = n+1
      nrint = nrint+1
      return
      end
      subroutine fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,
     * iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpu,fpv,
     * nru,nrv,wrk,lwrk)
c  given the set of function values z(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopdi determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) will be periodic in the variable v and will satisfy
c  the following constraints
c
c     s(tu(1),v) = dz(1) , tv(4) <=v<= tv(nv-3)
c
c  and (if iopt(2) = 1)
c
c     d s(tu(1),v)
c     ------------ =  dz(2)*cos(v)+dz(3)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  and (if iopt(3) = 1)
c
c     s(tu(nu),v)  =  0   tv(4) <=v<= tv(nv-3)
c
c  where the parameters dz(i) correspond to the derivative values g(i,j)
c  as defined in subroutine pogrid.
c
c  the b-spline coefficients of sp(u,v) are determined as the least-
c  squares solution  of an overdetermined linear system which depends
c  on the value of p and on the values dz(i),i=1,2,3. the correspond-
c  ing sum of squared residuals sq is a simple quadratic function in
c  the variables dz(i). these may or may not be provided. the values
c  dz(i) which are not given will be determined so as to minimize the
c  resulting sum of squared residuals sq. in that case the user must
c  provide some initial guess dz(i) and some estimate (dz(i)-step,
c  dz(i)+step) of the range of possible values for these latter.
c
c  sp(u,v) also depends on the parameter p (p>0) in such a way that
c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
c      with given knots, satisfying the constraints.
c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
c      satisfying the constraints.
c    - the function  f(p)=sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)
c      is continuous and strictly decreasing for p>0.
c
c  ..scalar arguments..
      integer ifsu,ifsv,ifbu,ifbv,mu,mv,mz,nu,nv,nuest,nvest,
     * nc,lwrk
      real z0,p,step,fp
c  ..array arguments..
      integer ider(2),nru(mu),nrv(mv),iopt(3)
      real u(mu),v(mv),z(mz),dz(3),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),
     * wrk(lwrk)
c  ..local scalars..
      real res,sq,sqq,step1,step2,three
      integer i,id0,iop0,iop1,i1,j,l,laa,lau,lav1,lav2,lbb,lbu,lbv,
     * lcc,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number
c  ..local arrays..
      integer nr(3)
      real delta(3),dzz(3),sum(3),a(6,6),g(6)
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpgrdi,fpsysy
c  ..
c  set constant
      three = 3
c  we partition the working space
      lsu = 1
      lsv = lsu+4*mu
      lri = lsv+4*mv
      mm = max0(nuest,mv+nvest)
      lq = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu = lav2+4*nvest
      lbv = lbu+5*nuest
      laa = lbv+5*nvest
      lbb = laa+2*mv
      lcc = lbb+2*nvest
      lcs = lcc+nvest
c  we calculate the smoothing spline sp(u,v) according to the input
c  values dz(i),i=1,2,3.
      iop0 = iopt(2)
      iop1 = iopt(3)
      call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     * wrk(lcc),wrk(lcs),nru,nrv)
      id0 = ider(1)
      if(id0.ne.0) go to 5
      res = (z0-dz(1))**2
      fp = fp+res
      sq = sq+res
c in case all derivative values dz(i) are given (step<=0) or in case
c we have spline interpolation, we accept this spline as a solution.
  5   if(step.le.0. .or. sq.le.0.) return
      dzz(1) = dz(1)
      dzz(2) = dz(2)
      dzz(3) = dz(3)
c number denotes the number of derivative values dz(i) that still must
c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 10
      number = 1
      nr(1) = 1
      delta(1) = step
  10  if(iop0.eq.0) go to 20
      if(ider(2).ne.0) go to 20
      step2 = step*three/tu(5)
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  20  if(number.eq.0) return
c the sum of squared residuals sq is a quadratic polynomial in the
c parameters g(j). we determine the unknown coefficients of this
c polymomial by calculating (number+1)*(number+2)/2 different splines
c according to specific values for g(j).
      do 30 i=1,number
         l = nr(i)
         step1 = delta(i)
         dzz(l) = dz(l)+step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *    wrk(lcc),wrk(lcs),nru,nrv)
         if(id0.eq.0) sum(i) = sum(i)+(z0-dzz(1))**2
         dzz(l) = dz(l)-step1
         call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *    wrk(lcc),wrk(lcs),nru,nrv)
         if(id0.eq.0) sqq = sqq+(z0-dzz(1))**2
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if(a(i,i).le.0.) go to 80
         g(i) = (sqq-sum(i))/(step1+step1)
         dzz(l) = dz(l)
  30  continue
      if(number.eq.1) go to 60
      do 50 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         dzz(l1) = dz(l1)+step1
         i1 = i-1
         do 40 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            dzz(l2) = dz(l2)+step2
            call fpgrdi(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,z,mz,dzz,
     *       iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *       wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *       wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     *       wrk(lcc),wrk(lcs),nru,nrv)
            if(id0.eq.0) sqq = sqq+(z0-dzz(1))**2
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            dzz(l2) = dz(l2)
  40     continue
         dzz(l1) = dz(l1)
  50  continue
c the optimal values g(j) are found as the solution of the system
c d (sq) / d (g(j)) = 0 , j=1,...,number.
  60  call fpsysy(a,number,g)
      do 70 i=1,number
         l = nr(i)
         dz(l) = dz(l)+g(i)
  70  continue
c we determine the spline sp(u,v) according to the optimal values g(j).
  80  call fpgrdi(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,z,mz,dz,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(laa),wrk(lbb),
     * wrk(lcc),wrk(lcs),nru,nrv)
      if(id0.eq.0) fp = fp+(z0-dz(1))**2
      return
      end
      subroutine fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,
     * iopt,ider,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpu,fpv,
     * nru,nrv,wrk,lwrk)
c  given the set of function values r(i,j) defined on the rectangular
c  grid (u(i),v(j)),i=1,2,...,mu;j=1,2,...,mv, fpopsp determines a
c  smooth bicubic spline approximation with given knots tu(i),i=1,..,nu
c  in the u-direction and tv(j),j=1,2,...,nv in the v-direction. this
c  spline sp(u,v) will be periodic in the variable v and will satisfy
c  the following constraints
c
c     s(tu(1),v) = dr(1) , tv(4) <=v<= tv(nv-3)
c
c     s(tu(nu),v) = dr(4) , tv(4) <=v<= tv(nv-3)
c
c  and (if iopt(2) = 1)
c
c     d s(tu(1),v)
c     ------------ =  dr(2)*cos(v)+dr(3)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  and (if iopt(3) = 1)
c
c     d s(tu(nu),v)
c     ------------- =  dr(5)*cos(v)+dr(6)*sin(v) , tv(4) <=v<= tv(nv-3)
c     d u
c
c  where the parameters dr(i) correspond to the derivative values at the
c  poles as defined in subroutine spgrid.
c
c  the b-spline coefficients of sp(u,v) are determined as the least-
c  squares solution  of an overdetermined linear system which depends
c  on the value of p and on the values dr(i),i=1,...,6. the correspond-
c  ing sum of squared residuals sq is a simple quadratic function in
c  the variables dr(i). these may or may not be provided. the values
c  dr(i) which are not given will be determined so as to minimize the
c  resulting sum of squared residuals sq. in that case the user must
c  provide some initial guess dr(i) and some estimate (dr(i)-step,
c  dr(i)+step) of the range of possible values for these latter.
c
c  sp(u,v) also depends on the parameter p (p>0) in such a way that
c    - if p tends to infinity, sp(u,v) becomes the least-squares spline
c      with given knots, satisfying the constraints.
c    - if p tends to zero, sp(u,v) becomes the least-squares polynomial,
c      satisfying the constraints.
c    - the function  f(p)=sumi=1,mu(sumj=1,mv((r(i,j)-sp(u(i),v(j)))**2)
c      is continuous and strictly decreasing for p>0.
c
c  ..scalar arguments..
      integer ifsu,ifsv,ifbu,ifbv,mu,mv,mr,nu,nv,nuest,nvest,
     * nc,lwrk
      real r0,r1,p,fp
c  ..array arguments..
      integer ider(4),nru(mu),nrv(mv),iopt(3)
      real u(mu),v(mv),r(mr),dr(6),tu(nu),tv(nv),c(nc),fpu(nu),fpv(nv),
     * wrk(lwrk),step(2)
c  ..local scalars..
      real res,sq,sqq,sq0,sq1,step1,step2,three
      integer i,id0,iop0,iop1,i1,j,l,lau,lav1,lav2,la0,la1,lbu,lbv,lb0,
     * lb1,lc0,lc1,lcs,lq,lri,lsu,lsv,l1,l2,mm,mvnu,number
c  ..local arrays..
      integer nr(6)
      real delta(6),drr(6),sum(6),a(6,6),g(6)
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpgrsp,fpsysy
c  ..
c  set constant
      three = 3
c  we partition the working space
      lsu = 1
      lsv = lsu+4*mu
      lri = lsv+4*mv
      mm = max0(nuest,mv+nvest)
      lq = lri+mm
      mvnu = nuest*(mv+nvest-8)
      lau = lq+mvnu
      lav1 = lau+5*nuest
      lav2 = lav1+6*nvest
      lbu = lav2+4*nvest
      lbv = lbu+5*nuest
      la0 = lbv+5*nvest
      la1 = la0+2*mv
      lb0 = la1+2*mv
      lb1 = lb0+2*nvest
      lc0 = lb1+2*nvest
      lc1 = lc0+nvest
      lcs = lc1+nvest
c  we calculate the smoothing spline sp(u,v) according to the input
c  values dr(i),i=1,...,6.
      iop0 = iopt(2)
      iop1 = iopt(3)
      id0 = ider(1)
      id1 = ider(3)
      call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     * wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      sq0 = 0.
      sq1 = 0.
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
c in case all derivative values dr(i) are given (step<=0) or in case
c we have spline interpolation, we accept this spline as a solution.
      if(sq.le.0.) return
      if(step(1).le.0. .and. step(2).le.0.) return
      do 10 i=1,6
        drr(i) = dr(i)
  10  continue
c number denotes the number of derivative values dr(i) that still must
c be optimized. let us denote these parameters by g(j),j=1,...,number.
      number = 0
      if(id0.gt.0) go to 20
      number = 1
      nr(1) = 1
      delta(1) = step(1)
  20  if(iop0.eq.0) go to 30
      if(ider(2).ne.0) go to 30
      step2 = step(1)*three/(tu(5)-tu(4))
      nr(number+1) = 2
      nr(number+2) = 3
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  30  if(id1.gt.0) go to 40
      number = number+1
      nr(number) = 4
      delta(number) = step(2)
  40  if(iop1.eq.0) go to 50
      if(ider(4).ne.0) go to 50
      step2 = step(2)*three/(tu(nu)-tu(nu-4))
      nr(number+1) = 5
      nr(number+2) = 6
      delta(number+1) = step2
      delta(number+2) = step2
      number = number+2
  50  if(number.eq.0) return
c the sum of squared residulas sq is a quadratic polynomial in the
c parameters g(j). we determine the unknown coefficients of this
c polymomial by calculating (number+1)*(number+2)/2 different splines
c according to specific values for g(j).
      do 60 i=1,number
         l = nr(i)
         step1 = delta(i)
         drr(l) = dr(l)+step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sum(i),fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *    wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sum(i) = sum(i)+sq0+sq1
         drr(l) = dr(l)-step1
         call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *    iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *    wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *    wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *    wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
         if(id0.eq.0) sq0 = (r0-drr(1))**2
         if(id1.eq.0) sq1 = (r1-drr(4))**2
         sqq = sqq+sq0+sq1
         drr(l) = dr(l)
         a(i,i) = (sum(i)+sqq-sq-sq)/step1**2
         if(a(i,i).le.0.) go to 110
         g(i) = (sqq-sum(i))/(step1+step1)
  60  continue
      if(number.eq.1) go to 90
      do 80 i=2,number
         l1 = nr(i)
         step1 = delta(i)
         drr(l1) = dr(l1)+step1
         i1 = i-1
         do 70 j=1,i1
            l2 = nr(j)
            step2 = delta(j)
            drr(l2) = dr(l2)+step2
            call fpgrsp(ifsu,ifsv,ifbu,ifbv,1,u,mu,v,mv,r,mr,drr,
     *       iop0,iop1,tu,nu,tv,nv,p,c,nc,sqq,fp,fpu,fpv,mm,mvnu,
     *       wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     *       wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     *       wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
            if(id0.eq.0) sq0 = (r0-drr(1))**2
            if(id1.eq.0) sq1 = (r1-drr(4))**2
            sqq = sqq+sq0+sq1
            a(i,j) = (sq+sqq-sum(i)-sum(j))/(step1*step2)
            drr(l2) = dr(l2)
  70     continue
         drr(l1) = dr(l1)
  80  continue
c the optimal values g(j) are found as the solution of the system
c d (sq) / d (g(j)) = 0 , j=1,...,number.
  90  call fpsysy(a,number,g)
      do 100 i=1,number
         l = nr(i)
         dr(l) = dr(l)+g(i)
 100  continue
c we determine the spline sp(u,v) according to the optimal values g(j).
 110  call fpgrsp(ifsu,ifsv,ifbu,ifbv,0,u,mu,v,mv,r,mr,dr,
     * iop0,iop1,tu,nu,tv,nv,p,c,nc,sq,fp,fpu,fpv,mm,mvnu,
     * wrk(lsu),wrk(lsv),wrk(lri),wrk(lq),wrk(lau),wrk(lav1),
     * wrk(lav2),wrk(lbu),wrk(lbv),wrk(la0),wrk(la1),wrk(lb0),
     * wrk(lb1),wrk(lc0),wrk(lc1),wrk(lcs),nru,nrv)
      if(id0.eq.0) sq0 = (r0-dr(1))**2
      if(id1.eq.0) sq1 = (r1-dr(4))**2
      sq = sq+sq0+sq1
      return
      end
      subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
c  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
c  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
c  to. for each panel a stack is constructed  containing the numbers
c  of data points lying inside; index(j),j=1,2,...,nreg points to the
c  first data point in the jth panel while nummer(i),i=1,2,...,m gives
c  the number of the next data point in the panel.
c  ..
c  ..scalar arguments..
      integer m,kx,ky,nx,ny,nreg
c  ..array arguments..
      real x(m),y(m),tx(nx),ty(ny)
      integer nummer(m),index(nreg)
c  ..local scalars..
      real xi,yi
      integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
c  ..
      kx1 = kx+1
      ky1 = ky+1
      nk1x = nx-kx1
      nk1y = ny-ky1
      nyy = nk1y-ky
      do 10 i=1,nreg
        index(i) = 0
  10  continue
      do 60 im=1,m
        xi = x(im)
        yi = y(im)
        l = kx1
        l1 = l+1
  20    if(xi.lt.tx(l1) .or. l.eq.nk1x) go to 30
        l = l1
        l1 = l+1
        go to 20
  30    k = ky1
        k1 = k+1
  40    if(yi.lt.ty(k1) .or. k.eq.nk1y) go to 50
        k = k1
        k1 = k+1
        go to 40
  50    num = (l-kx1)*nyy+k-ky
        nummer(im) = index(num)
        index(num) = im
  60  continue
      return
      end

      subroutine fppara(iopt,idim,m,u,mx,x,w,ub,ue,k,s,nest,tol,maxit,
     * k1,k2,n,t,nc,c,fp,fpint,z,a,b,g,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real ub,ue,s,tol,fp
      integer iopt,idim,m,mx,k,nest,maxit,k1,k2,n,nc,ier
c  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),fpint(nest),
     * z(nc),a(nest,k1),b(nest,k2),g(nest,k2),q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,con1,con4,con9,cos,fac,fpart,fpms,fpold,fp0,f1,f2,f3,
     * half,one,p,pinv,piv,p1,p2,p3,rn,sin,store,term,ui,wi
      integer i,ich1,ich3,it,iter,i1,i2,i3,j,jj,j1,j2,k3,l,l0,
     * mk1,new,nk1,nmax,nmin,nplus,npl1,nrint,n8
c  ..local arrays..
      real h(7),xi(10)
c  ..function references
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares curve sinf(u),    c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1 sinf(u) is the requested curve.                          c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares curve until finally fp<=s.         c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+k+1.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial curve of   c
c      degree k; n = nmin = 2*k+2                                      c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the polynomial curve of degree k.           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine nmin, the number of knots for polynomial approximation.
      nmin = 2*k1
      if(iopt.lt.0) go to 60
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for spline interpolation.
      nmax = m+k1
      if(s.gt.0.) go to 45
c  if s=0, s(u) is an interpolating curve.
c  test whether the required storage space exceeds the available one.
      n = nmax
      if(nmax.gt.nest) go to 420
c  find the position of the interior knots in case of interpolation.
  10  mk1 = m-k1
      if(mk1.eq.0) go to 60
      k3 = k/2
      i = k2
      j = k3+2
      if(k3*2.eq.k) go to 30
      do 20 l=1,mk1
        t(i) = u(j)
        i = i+1
        j = j+1
  20  continue
      go to 60
  30  do 40 l=1,mk1
        t(i) = (u(j)+u(j-1))*half
        i = i+1
        j = j+1
  40  continue
      go to 60
c  if s>0 our initial choice of knots depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial curve which is a spline curve without interior knots.
c  if iopt=1 and fp0>s we start computing the least squares spline curve
c  according to the set of knots found at the last call of the routine.
  45  if(iopt.eq.0) go to 50
      if(n.eq.nmin) go to 50
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 60
  50  n = nmin
      fpold = 0.
      nplus = 0
      nrdata(1) = m-2
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  60  do 200 iter = 1,m
        if(n.eq.nmin) ier = -2
c  find nrint, tne number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(u).
        nk1 = n-k1
        i = n
        do 70 j=1,k1
          t(j) = ub
          t(i) = ue
          i = i-1
  70    continue
c  compute the b-spline coefficients of the least-squares spline curve
c  sinf(u). the observation matrix a is built up row by row and
c  reduced to upper triangular form by givens transformations.
c  at the same time fp=f(p=inf) is computed.
        fp = 0.
c  initialize the b-spline coefficients and the observation matrix a.
        do 75 i=1,nc
          z(i) = 0.
  75    continue
        do 80 i=1,nk1
          do 80 j=1,k1
            a(i,j) = 0.
  80    continue
        l = k1
        jj = 0
        do 130 it=1,m
c  fetch the current data point u(it),x(it).
          ui = u(it)
          wi = w(it)
          do 83 j=1,idim
             jj = jj+1
             xi(j) = x(jj)*wi
  83      continue
c  search for knot interval t(l) <= ui < t(l+1).
  85      if(ui.lt.t(l+1) .or. l.eq.nk1) go to 90
          l = l+1
          go to 85
c  evaluate the (k+1) non-zero b-splines at ui and store them in q.
  90      call fpbspl(t,n,k,ui,l,h)
          do 95 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  95      continue
c  rotate the new row of the observation matrix into triangle.
          j = l-k1
          do 110 i=1,k1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 110
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cos,sin)
c  transformations to right hand side.
            j1 = j
            do 97 j2 =1,idim
               call fprota(cos,sin,xi(j2),z(j1))
               j1 = j1+n
  97        continue
            if(i.eq.k1) go to 120
            i2 = 1
            i3 = i+1
            do 100 i1 = i3,k1
              i2 = i2+1
c  transformations to left hand side.
              call fprota(cos,sin,h(i1),a(j,i2))
 100        continue
 110      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 120      do 125 j2=1,idim
             fp  = fp+xi(j2)**2
 125      continue
 130    continue
        if(ier.eq.(-2)) fp0 = fp
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients.
        j1 = 1
        do 135 j2=1,idim
           call fpback(a,z(j1),nk1,k1,c(j1),nest)
           j1 = j1+n
 135    continue
c  test whether the approximation sinf(u) is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 250
c  if n = nmax, sinf(u) is an interpolating spline curve.
        if(n.eq.nmax) go to 430
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of
c  the storage capacity limitation.
        if(n.eq.nest) go to 420
c  determine the number of knots nplus we are going to add.
        if(ier.eq.0) go to 140
        nplus = 1
        ier = 0
        go to 150
 140    npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
 150    fpold = fp
c  compute the sum of squared residuals for each knot interval
c  t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k2
        new = 0
        jj = 0
        do 180 it=1,m
          if(u(it).lt.t(l) .or. l.gt.nk1) go to 160
          new = 1
          l = l+1
 160      term = 0.
          l0 = l-k2
          do 175 j2=1,idim
            fac = 0.
            j1 = l0
            do 170 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 170        continue
            jj = jj+1
            term = term+(w(it)*(fac-x(jj)))**2
            l0 = l0+n
 175      continue
          fpart = fpart+term
          if(new.eq.0) go to 180
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 180    continue
        fpint(nrint) = fpart
        do 190 l=1,nplus
c  add a new knot.
          call fpknot(u,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation
          if(n.eq.nmax) go to 10
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 200
 190    continue
c  restart the computations with the new set of knots.
 200  continue
c  test whether the least-squares kth degree polynomial curve is a
c  solution of our approximation problem.
 250  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing spline curve sp(u).          c
c  **********************************************************          c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing curve     c
c  sp(u). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(u) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/p.            c
c  iteratively we then have to determine the value of p such that f(p),c
c  the sum of squared residuals be = s. we already know that the least c
c  squares kth degree polynomial curve corresponds to p=0, and that    c
c  the least-squares spline curve corresponds to p=infinity. the       c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
      call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 252 i=1,nk1
         p = p+a(i,1)
 252  continue
      rn = nk1
      p = rn/p
      ich1 = 0
      ich3 = 0
      n8 = n-nmin
c  iteration process to find the root of f(p) = s.
      do 360 iter=1,maxit
c  the rows of matrix b with weight 1/p are rotated into the
c  triangularised observation matrix a which is stored in g.
        pinv = one/p
        do 255 i=1,nc
          c(i) = z(i)
 255    continue
        do 260 i=1,nk1
          g(i,k2) = 0.
          do 260 j=1,k1
            g(i,j) = a(i,j)
 260    continue
        do 300 it=1,n8
c  the row of matrix b is rotated into triangle by givens transformation
          do 270 i=1,k2
            h(i) = b(it,i)*pinv
 270      continue
          do 275 j=1,idim
            xi(j) = 0.
 275      continue
          do 290 j=it,nk1
            piv = h(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g(j,1),cos,sin)
c  transformations to right hand side.
            j1 = j
            do 277 j2=1,idim
              call fprota(cos,sin,xi(j2),c(j1))
              j1 = j1+n
 277        continue
            if(j.eq.nk1) go to 300
            i2 = k1
            if(j.gt.n8) i2 = nk1-j
            do 280 i=1,i2
c  transformations to left hand side.
              i1 = i+1
              call fprota(cos,sin,h(i1),g(j,i1))
              h(i) = h(i1)
 280        continue
            h(i2+1) = 0.
 290      continue
 300    continue
c  backward substitution to obtain the b-spline coefficients.
        j1 = 1
        do 305 j2=1,idim
          call fpback(g,c(j1),nk1,k2,c(j1),nest)
          j1 =j1+n
 305    continue
c  computation of f(p).
        fp = 0.
        l = k2
        jj = 0
        do 330 it=1,m
          if(u(it).lt.t(l) .or. l.gt.nk1) go to 310
          l = l+1
 310      l0 = l-k2
          term = 0.
          do 325 j2=1,idim
            fac = 0.
            j1 = l0
            do 320 j=1,k1
              j1 = j1+1
              fac = fac+c(j1)*q(it,j)
 320        continue
            jj = jj+1
            term = term+(fac-x(jj))**2
            l0 = l0+n
 325      continue
          fp = fp+term*w(it)**2
 330    continue
c  test whether the approximation sp(u) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 340
        if((f2-f3).gt.acc) go to 335
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p=p1*con9 + p2*con1
        go to 360
 335    if(f2.lt.0.) ich3=1
 340    if(ich1.ne.0) go to 350
        if((f1-f2).gt.acc) go to 345
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 360
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 360
 345    if(f2.gt.0.) ich1=1
c  test whether the iteration process proceeds as theoretically
c  expected.
 350    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 360  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
 440  return
      end
      subroutine fppasu(iopt,ipar,idim,u,mu,v,mv,z,mz,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,
     * fpintv,lastdi,nplusu,nplusv,nru,nrv,nrdatu,nrdatv,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      real s,tol,fp,fp0,fpold,reducu,reducv
      integer iopt,idim,mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi,
     * nplusu,nplusv,lwrk,ier
c  ..array arguments..
      real u(mu),v(mv),z(mz*idim),tu(nuest),tv(nvest),c(nc*idim),
     * fpintu(nuest),fpintv(nvest),wrk(lwrk)
      integer ipar(2),nrdatu(nuest),nrdatv(nvest),nru(mu),nrv(mv)
c  ..local scalars
      real acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,con1,con9,con4,
     * peru,perv,ub,ue,vb,ve
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,iter,j,lau1,lav1,laa,
     * l,lau,lav,lbu,lbv,lq,lri,lsu,lsv,l1,l2,l3,l4,mm,mpm,mvnu,ncof,
     * nk1u,nk1v,nmaxu,nmaxv,nminu,nminv,nplu,nplv,npl1,nrintu,
     * nrintv,nue,nuk,nve,nuu,nvv
c  ..function references..
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpgrpa,fpknot
c  ..
c   set constants
      one = 1
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c  set boundaries of the approximation domain
      ub = u(1)
      ue = u(mu)
      vb = v(1)
      ve = v(mv)
c  we partition the working space.
      lsu = 1
      lsv = lsu+mu*4
      lri = lsv+mv*4
      mm = max0(nuest,mv)
      lq = lri+mm*idim
      mvnu = nuest*mv*idim
      lau = lq+mvnu
      nuk = nuest*5
      lbu = lau+nuk
      lav = lbu+nuk
      nuk = nvest*5
      lbv = lav+nuk
      laa = lbv+nuk
      lau1 = lau
      if(ipar(1).eq.0) go to 10
      peru = ue-ub
      lau1 = laa
      laa = laa+4*nuest
  10  lav1 = lav
      if(ipar(2).eq.0) go to 20
      perv = ve-vb
      lav1 = laa
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(u,v), c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1  sinf(u,v) is the requested approximation.               c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmaxu = mu+4+2*ipar(1) and  nmaxv = mv+4+2*ipar(2)   c
c    if s>0 and                                                        c
c     *iopt=0 we first compute the least-squares polynomial            c
c          nu=nminu=8 and nv=nminv=8                                   c
c     *iopt=1 we start with the knots found at the last call of the    c
c      routine, except for the case that s > fp0; then we can compute  c
c      the least-squares polynomial directly.                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine the number of knots for polynomial approximation.
  20  nminu = 8
      nminv = 8
      if(iopt.lt.0) go to 100
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  find nmaxu and nmaxv which denote the number of knots in u- and v-
c  direction in case of spline interpolation.
      nmaxu = mu+4+2*ipar(1)
      nmaxv = mv+4+2*ipar(2)
c  find nue and nve which denote the maximum number of knots
c  allowed in each direction
      nue = min0(nmaxu,nuest)
      nve = min0(nmaxv,nvest)
      if(s.gt.0.) go to 60
c  if s = 0, s(u,v) is an interpolating spline.
      nu = nmaxu
      nv = nmaxv
c  test whether the required storage space exceeds the available one.
      if(nv.gt.nvest .or. nu.gt.nuest) go to 420
c  find the position of the interior knots in case of interpolation.
c  the knots in the u-direction.
      nuu = nu-8
      if(nuu.eq.0) go to 40
      i = 5
      j = 3-ipar(1)
      do 30 l=1,nuu
        tu(i) = u(j)
        i = i+1
        j = j+1
  30  continue
c  the knots in the v-direction.
  40  nvv = nv-8
      if(nvv.eq.0) go to 60
      i = 5
      j = 3-ipar(2)
      do 50 l=1,nvv
        tv(i) = v(j)
        i = i+1
        j = j+1
  50  continue
      go to 100
c  if s > 0 our initial choice of knots depends on the value of iopt.
  60  if(iopt.eq.0) go to 90
      if(fp0.le.s) go to 90
c  if iopt=1 and fp0 > s we start computing the least- squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates u(i) inside each knot
c  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mpm = mu-1
      do 70 i=2,mpm
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 70
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
  70  continue
c  we determine the number of grid coordinates v(i) inside each knot
c  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      mpm = mv-1
      do 80 i=2,mpm
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 80
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
  80  continue
      go to 100
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial (which is a spline without interior knots).
  90  nu = nminu
      nv = nminv
      nrdatu(1) = mu-2
      nrdatv(1) = mv-2
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
 100  mpm = mu+mv
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
c  main loop for the different sets of knots.mpm=mu+mv is a save upper
c  bound for the number of trials.
      do 250 iter=1,mpm
        if(nu.eq.nminu .and. nv.eq.nminv) ier = -2
c  find nrintu (nrintv) which is the number of knot intervals in the
c  u-direction (v-direction).
        nrintu = nu-nminu+1
        nrintv = nv-nminv+1
c  find ncof, the number of b-spline coefficients for the current set
c  of knots.
        nk1u = nu-4
        nk1v = nv-4
        ncof = nk1u*nk1v
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
        if(ipar(1).ne.0) go to 110
        i = nu
        do 105 j=1,4
          tu(j) = ub
          tu(i) = ue
          i = i-1
 105    continue
        go to 120
 110    l1 = 4
        l2 = l1
        l3 = nu-3
        l4 = l3
        tu(l2) = ub
        tu(l3) = ue
        do 115 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tu(l2) = tu(l4)-peru
          tu(l3) = tu(l1)+peru
 115    continue
 120    if(ipar(2).ne.0) go to 130
        i = nv
        do 125 j=1,4
          tv(j) = vb
          tv(i) = ve
          i = i-1
 125    continue
        go to 140
 130    l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 135 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-perv
          tv(l3) = tv(l1)+perv
 135    continue
c  find the least-squares spline sinf(u,v) and calculate for each knot
c  interval tu(j+3)<=u<=tu(j+4) (tv(j+3)<=v<=tv(j+4)) the sum
c  of squared residuals fpintu(j),j=1,2,...,nu-7 (fpintv(j),j=1,2,...
c  ,nv-7) for the data points having their absciss (ordinate)-value
c  belonging to that interval.
c  fp gives the total sum of squared residuals.
 140    call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu,
     *  nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv),
     *  wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1),
     *  wrk(lbu),wrk(lbv),nru,nrv)
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nu=nmaxu and nv=nmaxv, sinf(u,v) is an interpolating spline.
        if(nu.eq.nmaxu .and. nv.eq.nmaxv) go to 430
c  increase the number of knots.
c  if nu=nue and nv=nve we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        ier = 0
c  adjust the parameter reducu or reducv according to the direction
c  in which the last added knots were located.
        if(lastdi) 150,170,160
 150    reducu = fpold-fp
        go to 170
 160    reducv = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.nminu) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 1
        if(nv.eq.nminv) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
 190    if(nplu-nplv) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
c  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        do 220 l=1,nplusu
c  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,1)
c  test whether we cannot further increase the number of knots in the
c  u-direction.
          if(nu.eq.nue) go to 250
 220    continue
        go to 250
 230    if(nv.eq.nve) go to 210
c  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
c  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
c  test whether we cannot further increase the number of knots in the
c  v-direction.
          if(nv.eq.nve) go to 250
 240    continue
c  restart the computations with the new set of knots.
 250  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(u,v)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
c  this smoothing spline varies with the parameter p in such a way thatc
c  f(p)=suml=1,idim(sumi=1,mu(sumj=1,mv((z(i,j,l)-sp(u(i),v(j),l))**2) c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. iteratively we then have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(u,v) and the corresponding sum of
c  squared residuals fp.
        call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu,
     *  nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv),
     *  wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1),
     *  wrk(lbu),wrk(lbv),nru,nrv)
c  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
      subroutine fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,
     * fp,fpint,z,a1,a2,b,g1,g2,q,nrdata,ier)
c  ..
c  ..scalar arguments..
      real s,tol,fp
      integer iopt,m,k,nest,maxit,k1,k2,n,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),fpint(nest),z(nest),
     * a1(nest,k1),a2(nest,k),b(nest,k2),g1(nest,k2),g2(nest,k1),
     * q(m,k1)
      integer nrdata(nest)
c  ..local scalars..
      real acc,cos,c1,d1,fpart,fpms,fpold,fp0,f1,f2,f3,p,per,pinv,piv,
     * p1,p2,p3,sin,store,term,wi,xi,yi,rn,one,con1,con4,con9,half
      integer i,ich1,ich3,ij,ik,it,iter,i1,i2,i3,j,jk,jper,j1,j2,kk,
     * kk1,k3,l,l0,l1,l5,mm,m1,new,nk1,nk2,nmax,nmin,nplus,npl1,
     * nrint,n10,n11,n7,n8
c  ..local arrays..
      real h(6),h1(7),h2(6)
c  ..function references..
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpbacp,fpbspl,fpgivs,fpdisc,fpknot,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1: determination of the number of knots and their position     c
c  **************************************************************      c
c  given a set of knots we compute the least-squares periodic spline   c
c  sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    c
c  the initial choice of knots depends on the value of s and iopt.     c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmax = m+2*k.                                        c
c    if s > 0 and                                                      c
c      iopt=0 we first compute the least-squares polynomial of         c
c      degree k; n = nmin = 2*k+2. since s(x) must be periodic we      c
c      find that s(x) is a constant function.                          c
c      iopt=1 we start with the set of knots found at the last         c
c      call of the routine, except for the case that s > fp0; then     c
c      we compute directly the least-squares periodic polynomial.      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      m1 = m-1
      kk = k
      kk1 = k1
      k3 = 3*k+1
      nmin = 2*k1
c  determine the length of the period of s(x).
      per = x(m)-x(1)
      if(iopt.lt.0) go to 50
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  determine nmax, the number of knots for periodic spline interpolation
      nmax = m+2*k
      if(s.gt.0. .or. nmax.eq.nmin) go to 30
c  if s=0, s(x) is an interpolating spline.
      n = nmax
c  test whether the required storage space exceeds the available one.
      if(n.gt.nest) go to 620
c  find the position of the interior knots in case of interpolation.
   5  if((k/2)*2 .eq. k) go to 20
      do 10 i=2,m1
        j = i+k
        t(j) = x(i)
  10  continue
      if(s.gt.0.) go to 50
      kk = k-1
      kk1 = k
      if(kk.gt.0) go to 50
      t(1) = t(m)-per
      t(2) = x(1)
      t(m+1) = x(m)
      t(m+2) = t(3)+per
      do 15 i=1,m1
        c(i) = y(i)
  15  continue
      c(m) = c(1)
      fp = 0.
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
      go to 630
  20  do 25 i=2,m1
        j = i+k
        t(j) = (x(i)+x(i-1))*half
  25  continue
      go to 50
c  if s > 0 our initial choice depends on the value of iopt.
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  periodic polynomial. (i.e. a constant function).
c  if iopt=1 and fp0>s we start computing the least-squares periodic
c  spline according the set of knots found at the last call of the
c  routine.
  30  if(iopt.eq.0) go to 35
      if(n.eq.nmin) go to 35
      fp0 = fpint(n)
      fpold = fpint(n-1)
      nplus = nrdata(n)
      if(fp0.gt.s) go to 50
c  the case that s(x) is a constant function is treated separetely.
c  find the least-squares constant c1 and compute fp0 at the same time.
  35  fp0 = 0.
      d1 = 0.
      c1 = 0.
      do 40 it=1,m1
        wi = w(it)
        yi = y(it)*wi
        call fpgivs(wi,d1,cos,sin)
        call fprota(cos,sin,yi,c1)
        fp0 = fp0+yi**2
  40  continue
      c1 = c1/d1
c  test whether that constant function is a solution of our problem.
      fpms = fp0-s
      if(fpms.lt.acc .or. nmax.eq.nmin) go to 640
      fpold = fp0
c  test whether the required storage space exceeds the available one.
      if(nmin.ge.nest) go to 620
c  start computing the least-squares periodic spline with one
c  interior knot.
      nplus = 1
      n = nmin+1
      mm = (m+1)/2
      t(k2) = x(mm)
      nrdata(1) = mm-2
      nrdata(2) = m1-mm
c  main loop for the different sets of knots. m is a save upper
c  bound for the number of trials.
  50  do 340 iter=1,m
c  find nrint, the number of knot intervals.
        nrint = n-nmin+1
c  find the position of the additional knots which are needed for
c  the b-spline representation of s(x). if we take
c      t(k+1) = x(1), t(n-k) = x(m)
c      t(k+1-j) = t(n-k-j) - per, j=1,2,...k
c      t(n-k+j) = t(k+1+j) + per, j=1,2,...k
c  then s(x) is a periodic spline with period per if the b-spline
c  coefficients satisfy the following conditions
c      c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
        t(k1) = x(1)
        nk1 = n-k1
        nk2 = nk1+1
        t(nk2) = x(m)
        do 60 j=1,k
          i1 = nk2+j
          i2 = nk2-j
          j1 = k1+j
          j2 = k1-j
          t(i1) = t(j1)+per
          t(j2) = t(i2)-per
  60    continue
c  compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
c  periodic spline sinf(x). the observation matrix a is built up row
c  by row while taking into account condition (**) and is reduced to
c  triangular form by givens transformations .
c  at the same time fp=f(p=inf) is computed.
c  the n7 x n7 triangularised upper matrix a has the form
c            ! a1 '    !
c        a = !    ' a2 !
c            ! 0  '    !
c  with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
c  matrix of bandwith k+1 ( n10 = n7-k).
c  initialization.
        do 70 i=1,nk1
          z(i) = 0.
          do 70 j=1,kk1
            a1(i,j) = 0.
  70    continue
        n7 = nk1-k
        n10 = n7-kk
        jper = 0
        fp = 0.
        l = k1
        do 290 it=1,m1
c  fetch the current data point x(it),y(it)
          xi = x(it)
          wi = w(it)
          yi = y(it)*wi
c  search for knot interval t(l) <= xi < t(l+1).
  80      if(xi.lt.t(l+1)) go to 85
          l = l+1
          go to 80
c  evaluate the (k+1) non-zero b-splines at xi and store them in q.
  85      call fpbspl(t,n,k,xi,l,h)
          do 90 i=1,k1
            q(it,i) = h(i)
            h(i) = h(i)*wi
  90      continue
          l5 = l-k1
c  test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
          if(l5.lt.n10) go to 285
          if(jper.ne.0) go to 160
c  initialize the matrix a2.
          do 95 i=1,n7
          do 95 j=1,kk
              a2(i,j) = 0.
  95      continue
          jk = n10+1
          do 110 i=1,kk
            ik = jk
            do 100 j=1,kk1
              if(ik.le.0) go to 105
              a2(ik,i) = a1(ik,j)
              ik = ik-1
 100        continue
 105        jk = jk+1
 110      continue
          jper = 1
c  if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
c  we take account of condition (**) for setting up the new row
c  of the observation matrix a. this row is stored in the arrays h1
c  (the part with respect to a1) and h2 (the part with
c  respect to a2).
 160      do 170 i=1,kk
            h1(i) = 0.
            h2(i) = 0.
 170      continue
          h1(kk1) = 0.
          j = l5-n10
          do 210 i=1,kk1
            j = j+1
            l0 = j
 180        l1 = l0-kk
            if(l1.le.0) go to 200
            if(l1.le.n10) go to 190
            l0 = l1-n10
            go to 180
 190        h1(l1) = h(i)
            go to 210
 200        h2(l0) = h2(l0)+h(i)
 210      continue
c  rotate the new row of the observation matrix into triangle
c  by givens transformations.
          if(n10.le.0) go to 250
c  rotation with the rows 1,2,...n10 of matrix a.
          do 240 j=1,n10
            piv = h1(1)
            if(piv.ne.0.) go to 214
            do 212 i=1,kk
              h1(i) = h1(i+1)
 212        continue
            h1(kk1) = 0.
            go to 240
c  calculate the parameters of the givens transformation.
 214        call fpgivs(piv,a1(j,1),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,z(j))
c  transformations to the left hand side with respect to a2.
            do 220 i=1,kk
              call fprota(cos,sin,h2(i),a2(j,i))
 220        continue
            if(j.eq.n10) go to 250
            i2 = min0(n10-j,kk)
c  transformations to the left hand side with respect to a1.
            do 230 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),a1(j,i1))
              h1(i) = h1(i1)
 230        continue
            h1(i1) = 0.
 240      continue
c  rotation with the rows n10+1,...n7 of matrix a.
 250      do 270 j=1,kk
            ij = n10+j
            if(ij.le.0) go to 270
            piv = h2(j)
            if(piv.eq.0.) go to 270
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a2(ij,j),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(ij))
            if(j.eq.kk) go to 280
            j1 = j+1
c  transformations to left hand side.
            do 260 i=j1,kk
              call fprota(cos,sin,h2(i),a2(ij,i))
 260        continue
 270      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 280      fp = fp+yi**2
          go to 290
c  rotation of the new row of the observation matrix into
c  triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero
c  at xi.
 285      j = l5
          do 140 i=1,kk1
            j = j+1
            piv = h(i)
            if(piv.eq.0.) go to 140
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a1(j,1),cos,sin)
c  transformations to right hand side.
            call fprota(cos,sin,yi,z(j))
            if(i.eq.kk1) go to 150
            i2 = 1
            i3 = i+1
c  transformations to left hand side.
            do 130 i1=i3,kk1
              i2 = i2+1
              call fprota(cos,sin,h(i1),a1(j,i2))
 130        continue
 140      continue
c  add contribution of this row to the sum of squares of residual
c  right hand sides.
 150      fp = fp+yi**2
 290    continue
        fpint(n) = fp0
        fpint(n-1) = fpold
        nrdata(n) = nplus
c  backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        call fpbacp(a1,a2,z,n7,kk,c,kk1,nest)
c  calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        do 295 i=1,k
          j = i+n7
          c(j) = c(i)
 295    continue
        if(iopt.lt.0) go to 660
c  test whether the approximation sinf(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  if f(p=inf) < s accept the choice of knots.
        if(fpms.lt.0.) go to 350
c  if n=nmax, sinf(x) is an interpolating spline.
        if(n.eq.nmax) go to 630
c  increase the number of knots.
c  if n=nest we cannot increase the number of knots because of the
c  storage capacity limitation.
        if(n.eq.nest) go to 620
c  determine the number of knots nplus we are going to add.
        npl1 = nplus*2
        rn = nplus
        if(fpold-fp.gt.acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min0(nplus*2,max0(npl1,nplus/2,1))
        fpold = fp
c  compute the sum(wi*(yi-s(xi))**2) for each knot interval
c  t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.
        i = 1
        l = k1
        do 320 it=1,m1
          if(x(it).lt.t(l)) go to 300
          new = 1
          l = l+1
 300      term = 0.
          l0 = l-k2
          do 310 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 310      continue
          term = (w(it)*(term-y(it)))**2
          fpart = fpart+term
          if(new.eq.0) go to 320
          if(l.gt.k2) go to 315
          fpint(nrint) = term
          new = 0
          go to 320
 315      store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
 320    continue
        fpint(nrint) = fpint(nrint)+fpart
        do 330 l=1,nplus
c  add a new knot
          call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
c  if n=nmax we locate the knots as for interpolation.
          if(n.eq.nmax) go to 5
c  test whether we cannot further increase the number of knots.
          if(n.eq.nest) go to 340
 330    continue
c  restart the computations with the new set of knots.
 340  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 2: determination of the smoothing periodic spline sp(x).       c
c  *************************************************************       c
c  we have determined the number of knots and their position.          c
c  we now compute the b-spline coefficients of the smoothing spline    c
c  sp(x). the observation matrix a is extended by the rows of matrix   c
c  b expressing that the kth derivative discontinuities of sp(x) at    c
c  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
c  ponding weights of these additional rows are set to 1/sqrt(p).      c
c  iteratively we then have to determine the value of p such that      c
c  f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      c
c  the least-squares constant function corresponds to p=0, and that    c
c  the least-squares periodic spline corresponds to p=infinity. the    c
c  iteration process which is proposed here, makes use of rational     c
c  interpolation. since f(p) is a convex and strictly decreasing       c
c  function of p, it can be approximated by a rational function        c
c  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
c  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
c  to calculate the new value of p such that r(p)=s. convergence is    c
c  guaranteed by taking f1>0 and f3<0.                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jump of the kth derivative of the
c  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
 350  call fpdisc(t,n,k2,b,nest)
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      n11 = n10-1
      n8 = n7-1
      p = 0.
      l = n7
      do 352 i=1,k
         j = k+1-i
         p = p+a2(l,j)
         l = l-1
         if(l.eq.0) go to 356
 352  continue
      do 354 i=1,n10
         p = p+a1(i,1)
 354  continue
 356  rn = n7
      p = rn/p
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p) = s.
      do 595 iter=1,maxit
c  form the matrix g  as the matrix a extended by the rows of matrix b.
c  the rows of matrix b with weight 1/p are rotated into
c  the triangularised observation matrix a.
c  after triangularisation our n7 x n7 matrix g takes the form
c            ! g1 '    !
c        g = !    ' g2 !
c            ! 0  '    !
c  with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
c  matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one/p
c  store matrix a into g
        do 360 i=1,n7
          c(i) = z(i)
          g1(i,k1) = a1(i,k1)
          g1(i,k2) = 0.
          g2(i,1) = 0.
          do 360 j=1,k
            g1(i,j) = a1(i,j)
            g2(i,j+1) = a2(i,j)
 360    continue
        l = n10
        do 370 j=1,k1
          if(l.le.0) go to 375
          g2(l,1) = a1(l,j)
          l = l-1
 370    continue
 375    do 540 it=1,n8
c  fetch a new row of matrix b and store it in the arrays h1 (the part
c  with respect to g1) and h2 (the part with respect to g2).
          yi = 0.
          do 380 i=1,k1
            h1(i) = 0.
            h2(i) = 0.
 380      continue
          h1(k2) = 0.
          if(it.gt.n11) go to 420
          l = it
          l0 = it
          do 390 j=1,k2
            if(l0.eq.n10) go to 400
            h1(j) = b(it,j)*pinv
            l0 = l0+1
 390      continue
          go to 470
 400      l0 = 1
          do 410 l1=j,k2
            h2(l0) = b(it,l1)*pinv
            l0 = l0+1
 410      continue
          go to 470
 420      l = 1
          i = it-n10
          do 460 j=1,k2
            i = i+1
            l0 = i
 430        l1 = l0-k1
            if(l1.le.0) go to 450
            if(l1.le.n11) go to 440
            l0 = l1-n11
            go to 430
 440        h1(l1) = b(it,j)*pinv
            go to 460
 450        h2(l0) = h2(l0)+b(it,j)*pinv
 460      continue
          if(n11.le.0) go to 510
c  rotate this row into triangle by givens transformations without
c  square roots.
c  rotation with the rows l,l+1,...n11.
 470      do 500 j=l,n11
            piv = h1(1)
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,g1(j,1),cos,sin)
c  transformation to right hand side.
            call fprota(cos,sin,yi,c(j))
c  transformation to the left hand side with respect to g2.
            do 480 i=1,k1
              call fprota(cos,sin,h2(i),g2(j,i))
 480        continue
            if(j.eq.n11) go to 510
            i2 = min0(n11-j,k1)
c  transformation to the left hand side with respect to g1.
            do 490 i=1,i2
              i1 = i+1
              call fprota(cos,sin,h1(i1),g1(j,i1))
              h1(i) = h1(i1)
 490        continue
            h1(i1) = 0.
 500      continue
c  rotation with the rows n11+1,...n7
 510      do 530 j=1,k1
            ij = n11+j
            if(ij.le.0) go to 530
            piv = h2(j)
c  calculate the parameters of the givens transformation
            call fpgivs(piv,g2(ij,j),cos,sin)
c  transformation to the right hand side.
            call fprota(cos,sin,yi,c(ij))
            if(j.eq.k1) go to 540
            j1 = j+1
c  transformation to the left hand side.
            do 520 i=j1,k1
              call fprota(cos,sin,h2(i),g2(ij,i))
 520        continue
 530      continue
 540    continue
c  backward substitution to obtain the b-spline coefficients
c  c(j),j=1,2,...n7 of sp(x).
        call fpbacp(g1,g2,c,n7,k1,c,k2,nest)
c  calculate from condition (**) the b-spline coefficients c(n7+j),j=1,.
        do 545 i=1,k
          j = i+n7
          c(j) = c(i)
 545    continue
c  computation of f(p).
        fp = 0.
        l = k1
        do 570 it=1,m1
          if(x(it).lt.t(l)) go to 550
          l = l+1
 550      l0 = l-k2
          term = 0.
          do 560 j=1,k1
            l0 = l0+1
            term = term+c(l0)*q(it,j)
 560      continue
          fp = fp+(w(it)*(term-y(it)))**2
 570    continue
c  test whether the approximation sp(x) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 660
c  test whether the maximal number of iterations is reached.
        if(iter.eq.maxit) go to 600
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 580
        if((f2-f3) .gt. acc) go to 575
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 +p2*con1
        go to 595
 575    if(f2.lt.0.) ich3 = 1
 580    if(ich1.ne.0) go to 590
        if((f1-f2) .gt. acc) go to 585
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 595
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 595
 585    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 590    if(f2.ge.f1 .or. f2.le.f3) go to 610
c  find the new value for p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 595  continue
c  error codes and messages.
 600  ier = 3
      go to 660
 610  ier = 2
      go to 660
 620  ier = 1
      go to 660
 630  ier = -1
      go to 660
 640  ier = -2
c  the least-squares constant function c1 is a solution of our problem.
c  a constant function is a spline of degree k with all b-spline
c  coefficients equal to that constant c1.
      do 650 i=1,k1
        rn = k1-i
        t(i) = x(1)-rn*per
        c(i) = c1
        j = i+k1
        rn = i-1
        t(j) = x(m)+rn*per
 650  continue
      n = nmin
      fp = fp0
      fpint(n) = fp0
      fpint(n-1) = 0.
      nrdata(n) = 0
 660  return
      end
      subroutine fppocu(idim,k,a,b,ib,db,nb,ie,de,ne,cp,np)
c  subroutine fppocu finds a idim-dimensional polynomial curve p(u) =
c  (p1(u),p2(u),...,pidim(u)) of degree k, satisfying certain derivative
c  constraints at the end points a and b, i.e.
c                  (l)
c    if ib > 0 : pj   (a) = db(idim*l+j), l=0,1,...,ib-1
c                  (l)
c    if ie > 0 : pj   (b) = de(idim*l+j), l=0,1,...,ie-1
c
c  the polynomial curve is returned in its b-spline representation
c  ( coefficients cp(j), j=1,2,...,np )
c  ..
c  ..scalar arguments..
      integer idim,k,ib,nb,ie,ne,np
      real a,b
c  ..array arguments..
      real db(nb),de(ne),cp(np)
c  ..local scalars..
      real ab,aki
      integer i,id,j,jj,l,ll,k1,k2
c  ..local array..
      real work(6,6)
c  ..
      k1 = k+1
      k2 = 2*k1
      ab = b-a
      do 110 id=1,idim
        do 10 j=1,k1
          work(j,1) = 0.
  10    continue
        if(ib.eq.0) go to 50
        l = id
        do 20 i=1,ib
          work(1,i) = db(l)
          l = l+idim
  20    continue
        if(ib.eq.1) go to 50
        ll = ib
        do 40 j=2,ib
          ll =  ll-1
          do 30 i=1,ll
            aki = k1-i
            work(j,i) = ab*work(j-1,i+1)/aki + work(j-1,i)
  30      continue
  40    continue
  50    if(ie.eq.0) go to 90
        l = id
        j = k1
        do 60 i=1,ie
          work(j,i) = de(l)
          l = l+idim
          j = j-1
  60    continue
        if(ie.eq.1) go to 90
        ll = ie
        do 80 jj=2,ie
          ll =  ll-1
          j = k1+1-jj
          do 70 i=1,ll
            aki = k1-i
            work(j,i) = work(j+1,i) - ab*work(j,i+1)/aki
            j = j-1
  70      continue
  80    continue
  90    l = (id-1)*k2
        do 100 j=1,k1
          l = l+1
          cp(l) = work(j,1)
 100    continue
 110  continue
      return
      end
      subroutine fppogr(iopt,ider,u,mu,v,mv,z,mz,z0,r,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,
     * fpintv,dz,step,lastdi,nplusu,nplusv,lasttu,nru,nrv,nrdatu,
     * nrdatv,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      integer mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi,nplusu,nplusv,
     * lasttu,lwrk,ier
      real z0,r,s,tol,fp,fp0,fpold,reducu,reducv,step
c  ..array arguments..
      integer iopt(3),ider(2),nrdatu(nuest),nrdatv(nvest),nru(mu),
     * nrv(mv)
      real u(mu),v(mv),z(mz),tu(nuest),tv(nvest),c(nc),fpintu(nuest),
     * fpintv(nvest),dz(3),wrk(lwrk)
c  ..local scalars..
      real acc,fpms,f1,f2,f3,p,per,pi,p1,p2,p3,vb,ve,zmax,zmin,rn,one,
     * con1,con4,con9
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju,
     * ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu,
     * nrintv,nue,numax,nve,nvmax
c  ..local arrays..
      integer idd(2)
      real dzz(3)
c  ..function references..
      real abs,atan2,fprati
      integer max0,min0
c  ..subroutine references..
c    fpknot,fpopdi
c  ..
c   set constants
      one = 1
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      mumin = 4-iopt(3)
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      pi = atan2(0.,-one)
      per = pi+pi
      vb = v(1)
      ve = vb+per
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(u,v)  c
c  and the corresponding sum of squared residuals fp = f(p=inf).       c
c  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c
c  if iopt(1)>=0  we check whether we can accept the knots:            c
c    if fp <= s we will continue with the current set of knots.        c
c    if fp >  s we will increase the number of knots and compute the   c
c       corresponding least-squares spline until finally fp <= s.      c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c     knots in the u-direction equals nu=numax=mu+5+iopt(2)+iopt(3)    c
c     and in the v-direction nv=nvmax=mv+7.                            c
c    if s>0 and                                                        c
c      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
c       spline without interior knots : nu=8 ; nv=8.                   c
c      iopt(1)=1 we start with the set of knots found at the last call c
c       of the routine, except for the case that s > fp0; then we      c
c       compute the least-squares polynomial directly.                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iopt(1).lt.0) go to 120
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+5+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)
      if(s.gt.0.) go to 100
c  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax
c  test whether the required storage space exceeds the available one.
      if(nu.gt.nuest .or. nv.gt.nvest) go to 420
c  find the position of the knots in the v-direction.
      do 10 l=1,mv
        tv(l+3) = v(l)
  10  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 20 i=1,3
         tv(i) = v(l1)-per
         tv(l2) = v(i+1)+per
         l1 = l1+1
         l2 = l2+1
  20  continue
c  if not all the derivative values g(i,j) are given, we will first
c  estimate these values by computing a least-squares spline
      idd(1) = ider(1)
      if(idd(1).eq.0) idd(1) = 1
      if(idd(1).gt.0) dz(1) = z0
      idd(2) = ider(2)
      if(ider(1).lt.0) go to 30
      if(iopt(2).eq.0 .or. ider(2).ne.0) go to 70
c we set up the knots in the u-direction for computing the least-squares
c spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do 40 i=1,mu
         if(i1.gt.i2) go to 50
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  40  continue
  50  do 60 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = r
  60  continue
c we compute the least-squares spline for estimating the derivatives.
      call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      ifsu = 0
c if all the derivatives at the origin are known, we compute the
c interpolating spline.
c we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn.eq.0) go to 95
      ju = 2-iopt(2)
      do 80 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  80  continue
      nu = numax
      l = nu
      do 90 i=1,4
         tu(i) = 0.
         tu(l) = r
         l = l-1
  90  continue
c we compute the interpolating spline.
  95  call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      go to 430
c  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = 0
      if(iopt(1).eq.0) go to 115
      step = -step
      if(fp0.le.s) go to 115
c  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates u(i) inside each knot
c  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mu0 = 2-iopt(2)
      mu1 = mu-2+iopt(3)
      do 105 i=mu0,mu1
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 105
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
 105  continue
c  we determine the number of grid coordinates v(i) inside each knot
c  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      do 110 i=2,mv
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 110
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
 110  continue
      idd(1) = ider(1)
      idd(2) = ider(2)
      go to 120
c  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
c  squares polynomial (which is a spline without interior knots).
 115  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      nu = 8
      nv = 8
      nrdatu(1) = mu-3+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
c  main loop for the different sets of knots.mpm=mu+mv is a save upper
c  bound for the number of trials.
 120  mpm = mu+mv
      do 270 iter=1,mpm
c  find nrintu (nrintv) which is the number of knot intervals in the
c  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
        i = nu
        do 130 j=1,4
          tu(j) = 0.
          tu(i) = r
          i = i-1
 130    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 140 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-per
          tv(l3) = tv(l1)+per
 140    continue
c  find an estimate of the range of possible values for the optimal
c  derivatives at the origin.
        ktu = nrdatu(1)+2-iopt(2)
        if(nrintu.eq.1) ktu = mu
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lasttu) go to 150
         zmin = z0
         zmax = z0
         l = mv*ktu
         do 145 i=1,l
            if(z(i).lt.zmin) zmin = z(i)
            if(z(i).gt.zmax) zmax = z(i)
 145     continue
         step = zmax-zmin
         lasttu = ktu
c  find the least-squares spline sinf(u,v).
 150    call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dz,iopt,idd,
     *   tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *   wrk,lwrk)
        if(step.lt.0.) step = -step
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt(1).lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if(nu.eq.numax .and. nv.eq.nvmax) go to 430
c  increase the number of knots.
c  if nu=nue and nv=nve we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        if(ider(1).eq.0) fpintu(1) = fpintu(1)+(z0-c(1))**2
        ier = 0
c  adjust the parameter reducu or reducv according to the direction
c  in which the last added knots were located.
        if(lastdi) 160,155,170
 155     nplv = 3
         idd(2) = ider(2)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.8) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv.eq.8) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
c  test whether we are going to add knots in the u- or v-direction.
 190    if(nplu-nplv) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
c  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        istart = 0
        if(iopt(2).eq.0) istart = 1
        do 220 l=1,nplusu
c  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)
c  test whether we cannot further increase the number of knots in the
c  u-direction.
          if(nu.eq.nue) go to 270
 220    continue
        go to 270
 230    if(nv.eq.nve) go to 210
c  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
c  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
c  test whether we cannot further increase the number of knots in the
c  v-direction.
          if(nv.eq.nve) go to 270
 240    continue
c  restart the computations with the new set of knots.
 270  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(u,v)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
c  this smoothing spline depends on the parameter p in such a way that c
c    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. then iteratively we have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      dzz(1) = dz(1)
      dzz(2) = dz(2)
      dzz(3) = dz(3)
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(u,v) and the corresponding sum f(p).
        call fpopdi(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,z,mz,z0,dzz,iopt,idd,
     *   tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *   wrk,lwrk)
c  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
      subroutine fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s,nuest,nvest,
     * eta,tol,maxit,ib1,ib3,nc,ncc,intest,nrest,nu,tu,nv,tv,c,fp,sup,
     * fpint,coord,f,ff,row,cs,cosi,a,q,bu,bv,spu,spv,h,index,nummer,
     * wrk,lwrk,ier)
c  ..scalar arguments..
      integer iopt1,iopt2,iopt3,m,nuest,nvest,maxit,ib1,ib3,nc,ncc,
     * intest,nrest,nu,nv,lwrk,ier
      real s,eta,tol,fp,sup
c  ..array arguments..
      integer index(nrest),nummer(m)
      real u(m),v(m),z(m),w(m),tu(nuest),tv(nvest),c(nc),fpint(intest),
     * coord(intest),f(ncc),ff(nc),row(nvest),cs(nvest),cosi(5,nvest),
     * a(ncc,ib1),q(ncc,ib3),bu(nuest,5),bv(nvest,5),spu(m,4),spv(m,4),
     * h(ib3),wrk(lwrk)
c  ..user supplied function..
      real rad
c  ..local scalars..
      real acc,arg,co,c1,c2,c3,c4,dmax,eps,fac,fac1,fac2,fpmax,fpms,
     * f1,f2,f3,hui,huj,p,pi,pinv,piv,pi2,p1,p2,p3,r,ratio,si,sigma,
     * sq,store,uu,u2,u3,wi,zi,rn,one,two,three,con1,con4,con9,half,ten
      integer i,iband,iband3,iband4,ich1,ich3,ii,il,in,ipar,ipar1,irot,
     * iter,i1,i2,i3,j,jl,jrot,j1,j2,k,l,la,lf,lh,ll,lu,lv,lwest,l1,l2,
     * l3,l4,ncof,ncoff,nvv,nv4,nreg,nrint,nrr,nr1,nuu,nu4,num,num1,
     * numin,nvmin,rank,iband1
c  ..local arrays..
      real hu(4),hv(4)
c  ..function references..
      real abs,atan,cos,fprati,sin,sqrt
      integer min0
c  ..subroutine references..
c    fporde,fpbspl,fpback,fpgivs,fprota,fprank,fpdisc,fprppo
c  ..
c  set constants
      one = 1
      two = 2
      three = 3
      ten = 10
      half = 0.5e0
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      pi = atan(one)*4
      pi2 = pi+pi
      ipar = iopt2*(iopt2+3)/2
      ipar1 = ipar+1
      eps = sqrt(eta)
      if(iopt1.lt.0) go to 90
      numin = 9
      nvmin = 9+iopt2*(iopt2+1)
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt1.eq.0) go to 10
      if(s.lt.sup) if(nv-nvmin) 70,90,90
c  if iopt1 = 0 we begin by computing the weighted least-squares
c  polymomial of the form
c     s(u,v) = f(1)*(1-u**3)+f(2)*u**3+f(3)*(u**2-u**3)+f(4)*(u-u**3)
c  where f(4) = 0 if iopt2> 0 , f(3) = 0 if iopt2 > 1 and
c        f(2) = 0 if iopt3> 0.
c  the corresponding weighted sum of squared residuals gives the upper
c  bound sup for the smoothing factor s.
  10  sup = 0.
      do 20 i=1,4
         f(i) = 0.
         do 20 j=1,4
            a(i,j) = 0.
 20   continue
      do 50 i=1,m
         wi = w(i)
         zi = z(i)*wi
         uu = u(i)
         u2 = uu*uu
         u3 = uu*u2
         h(1) = (one-u3)*wi
         h(2) = u3*wi
         h(3) = u2*(one-uu)*wi
         h(4) = uu*(one-u2)*wi
         if(iopt3.ne.0) h(2) = 0.
         if(iopt2.gt.1) h(3) = 0.
         if(iopt2.gt.0) h(4) = 0.
         do 40 j=1,4
            piv = h(j)
            if(piv.eq.0.) go to 40
            call fpgivs(piv,a(j,1),co,si)
            call fprota(co,si,zi,f(j))
            if(j.eq.4) go to 40
            j1 = j+1
            j2 = 1
            do 30 l=j1,4
               j2 = j2+1
               call fprota(co,si,h(l),a(j,j2))
  30        continue
  40     continue
         sup = sup+zi*zi
  50  continue
      if(a(4,1).ne.0.) f(4) = f(4)/a(4,1)
      if(a(3,1).ne.0.) f(3) = (f(3)-a(3,2)*f(4))/a(3,1)
      if(a(2,1).ne.0.) f(2) = (f(2)-a(2,2)*f(3)-a(2,3)*f(4))/a(2,1)
      if(a(1,1).ne.0.)
     * f(1) = (f(1)-a(1,2)*f(2)-a(1,3)*f(3)-a(1,4)*f(4))/a(1,1)
c  find the b-spline representation of this least-squares polynomial
      c1 = f(1)
      c4 = f(2)
      c2 = f(4)/three+c1
      c3 = (f(3)+two*f(4))/three+c1
      nu = 8
      nv = 8
      do 60 i=1,4
         c(i) = c1
         c(i+4) = c2
         c(i+8) = c3
         c(i+12) = c4
         tu(i) = 0.
         tu(i+4) = one
         rn = 2*i-9
         tv(i) = rn*pi
         rn = 2*i-1
         tv(i+4) = rn*pi
  60  continue
      fp = sup
c  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms.lt.acc) go to 960
c  test whether we cannot further increase the number of knots.
  70  if(nuest.lt.numin .or. nvest.lt.nvmin) go to 950
c  find the initial set of interior knots of the spline in case iopt1=0.
      nu = numin
      nv = nvmin
      tu(5) = half
      nvv = nv-8
      rn = nvv+1
      fac = pi2/rn
      do 80 i=1,nvv
         rn = i
         tv(i+4) = rn*fac-pi
  80  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1 : computation of least-squares bicubic splines.              c
c  ******************************************************              c
c  if iopt1<0 we compute the least-squares bicubic spline according    c
c  to the given set of knots.                                          c
c  if iopt1>=0 we compute least-squares bicubic splines with in-       c
c  creasing numbers of knots until the corresponding sum f(p=inf)<=s.  c
c  the initial set of knots then depends on the value of iopt1         c
c    if iopt1=0 we start with one interior knot in the u-direction     c
c              (0.5) and 1+iopt2*(iopt2+1) in the v-direction.         c
c    if iopt1>0 we start with the set of knots found at the last       c
c              call of the routine.                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  90  do 570 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
         l1 = 4
         l2 = l1
         l3 = nv-3
         l4 = l3
         tv(l2) = -pi
         tv(l3) = pi
         do 120 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tv(l2) = tv(l4)-pi2
            tv(l3) = tv(l1)+pi2
 120     continue
        l = nu
        do 130 i=1,4
          tu(i) = 0.
          tu(l) = one
          l = l-1
 130    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        nuu = nu-7
        nvv = nv-7
        nrr = nvv/2
        nr1 = nrr+1
        nrint = nuu+nvv
        nreg = nuu*nvv
c  arrange the data points according to the panel they belong to.
        call fporde(u,v,m,3,3,tu,nu,tv,nv,nummer,index,nreg)
        if(iopt2.eq.0) go to 195
c  find the b-spline coefficients cosi of the cubic spline
c  approximations for cr(v)=rad(v)*cos(v) and sr(v) = rad(v)*sin(v)
c  if iopt2=1, and additionally also for cr(v)**2,sr(v)**2 and
c  2*cr(v)*sr(v) if iopt2=2
        do 140 i=1,nvv
           do 135 j=1,ipar
              cosi(j,i) = 0.
 135       continue
           do 140 j=1,nvv
              a(i,j) = 0.
 140    continue
c  the coefficients cosi are obtained from interpolation conditions
c  at the knots tv(i),i=4,5,...nv-4.
        do 175 i=1,nvv
           l2 = i+3
           arg = tv(l2)
           call fpbspl(tv,nv,3,arg,l2,hv)
           do 145 j=1,nvv
              row(j) = 0.
 145       continue
           ll = i
           do 150 j=1,3
              if(ll.gt.nvv) ll= 1
              row(ll) = row(ll)+hv(j)
              ll = ll+1
 150       continue
           co = cos(arg)
           si = sin(arg)
           r = rad(arg)
           cs(1) = co*r
           cs(2) = si*r
           if(iopt2.eq.1) go to 155
           cs(3) = cs(1)*cs(1)
           cs(4) = cs(2)*cs(2)
           cs(5) = cs(1)*cs(2)
 155       do 170 j=1,nvv
              piv = row(j)
              if(piv.eq.0.) go to 170
              call fpgivs(piv,a(j,1),co,si)
              do 160 l=1,ipar
                 call fprota(co,si,cs(l),cosi(l,j))
 160          continue
              if(j.eq.nvv) go to 175
              j1 = j+1
              j2 = 1
              do 165 l=j1,nvv
                 j2 = j2+1
                 call fprota(co,si,row(l),a(j,j2))
 165          continue
 170       continue
 175    continue
         do 190 l=1,ipar
            do 180 j=1,nvv
               cs(j) = cosi(l,j)
 180        continue
            call fpback(a,cs,nvv,nvv,cs,ncc)
            do 185 j=1,nvv
               cosi(l,j) = cs(j)
 185        continue
 190     continue
c  find ncof, the dimension of the spline and ncoff, the number
c  of coefficients in the standard b-spline representation.
 195    nu4 = nu-4
        nv4 = nv-4
        ncoff = nu4*nv4
        ncof = ipar1+nvv*(nu4-1-iopt2-iopt3)
c  find the bandwidth of the observation matrix a.
        iband = 4*nvv
        if(nuu-iopt2-iopt3.le.1) iband = ncof
        iband1 = iband-1
c  initialize the observation matrix a.
        do 200 i=1,ncof
          f(i) = 0.
          do 200 j=1,iband
            a(i,j) = 0.
 200    continue
c  initialize the sum of squared residuals.
        fp = 0.
        ratio = one+tu(6)/tu(5)
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 380 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lu = num1/nvv
          l1 = lu+4
          lv = num1-lu*nvv+1
          l2 = lv+3
          jrot = 0
          if(lu.gt.iopt2) jrot = ipar1+(lu-iopt2-1)*nvv
          lu = lu+1
c  test whether there are still data points in the current panel.
          in = index(num)
 210      if(in.eq.0) go to 380
c  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
c  evaluate for the u-direction, the 4 non-zero b-splines at u(in)
          call fpbspl(tu,nu,3,u(in),l1,hu)
c  evaluate for the v-direction, the 4 non-zero b-splines at v(in)
          call fpbspl(tv,nv,3,v(in),l2,hv)
c  store the value of these b-splines in spu and spv resp.
          do 220 i=1,4
            spu(in,i) = hu(i)
            spv(in,i) = hv(i)
 220      continue
c  initialize the new row of observation matrix.
          do 240 i=1,iband
            h(i) = 0.
 240      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in u- and v-direction and
c  by taking into account the conditions of the splines.
          do 250 i=1,nvv
             row(i) = 0.
 250      continue
c  take into account the periodicity condition of the bicubic splines.
          ll = lv
          do 260 i=1,4
             if(ll.gt.nvv) ll=1
             row(ll) = row(ll)+hv(i)
             ll = ll+1
 260      continue
c  take into account the other conditions of the splines.
          if(iopt2.eq.0 .or. lu.gt.iopt2+1) go to 280
          do 270 l=1,ipar
             cs(l) = 0.
             do 270 i=1,nvv
                cs(l) = cs(l)+row(i)*cosi(l,i)
 270     continue
c  fill in the non-zero elements of the new row.
 280     j1 = 0
         do 330 j =1,4
            jlu = j+lu
            huj = hu(j)
            if(jlu.gt.iopt2+2) go to 320
            go to (290,290,300,310),jlu
 290        h(1) = huj
            j1 = 1
            go to 330
 300        h(1) = h(1)+huj
            h(2) = huj*cs(1)
            h(3) = huj*cs(2)
            j1 = 3
            go to 330
 310        h(1) = h(1)+huj
            h(2) = h(2)+huj*ratio*cs(1)
            h(3) = h(3)+huj*ratio*cs(2)
            h(4) = huj*cs(3)
            h(5) = huj*cs(4)
            h(6) = huj*cs(5)
            j1 = 6
            go to 330
 320        if(jlu.gt.nu4 .and. iopt3.ne.0) go to 330
            do 325 i=1,nvv
               j1 = j1+1
               h(j1) = row(i)*huj
 325        continue
 330      continue
          do 335 i=1,iband
            h(i) = h(i)*wi
 335      continue
c  rotate the row into triangle by givens transformations.
          irot = jrot
          do 350 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 350
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the right hand side.
            call fprota(co,si,zi,f(irot))
            if(i.eq.iband) go to 360
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 340 j=i3,iband
              i2 = i2+1
              call fprota(co,si,h(j),a(irot,i2))
 340        continue
 350      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 360      fp = fp+zi**2
c  find the number of the next data point in the panel.
 370      in = nummer(in)
          go to 210
 380    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 390 i=1,ncof
          if(a(i,1).le.dmax) go to 390
          dmax = a(i,1)
 390    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 400 i=1,ncof
          if(a(i,1).le.sigma) go to 410
 400    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,ncc)
        rank = ncof
        do 405 i=1,ncof
          q(i,1) = a(i,1)/dmax
 405    continue
        go to 430
c  in case of rank deficiency, find the minimum norm solution.
 410    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband
        do 420 i=1,ncof
          ff(i) = f(i)
          do 420 j=1,iband
            q(i,j) = a(i,j)
 420    continue
        call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
        do 425 i=1,ncof
          q(i,1) = q(i,1)/dmax
 425    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
c  find the coefficients in the standard b-spline representation of
c  the spline.
 430    call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff)
c  test whether the least-squares spline is an acceptable solution.
        if(iopt1.lt.0) if(fp) 970,970,980
        fpms = fp-s
        if(abs(fpms).le.acc) if(fp) 970,970,980
c  if f(p=inf) < s, accept the choice of knots.
        if(fpms.lt.0.) go to 580
c  test whether we cannot further increase the number of knots
        if(m.lt.ncof) go to 935
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
 440    do 450 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 450    continue
        do 490 num=1,nreg
          num1 = num-1
          lu = num1/nvv
          l1 = lu+1
          lv = num1-lu*nvv
          l2 = lv+1+nuu
          jrot = lu*nv4+lv
          in = index(num)
 460      if(in.eq.0) go to 490
          store = 0.
          i1 = jrot
          do 480 i=1,4
            hui = spu(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hui*spv(in,j)*c(j1)
 470        continue
            i1 = i1+nv4
 480      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*u(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*v(in)
          in = nummer(in)
          go to 460
 490    continue
c bring together the information concerning knot panels which are
c symmetric with respect to the origin.
        do 495 i=1,nrr
          l1 = nuu+i
          l2 = l1+nrr
          fpint(l1) = fpint(l1)+fpint(l2)
          coord(l1) = coord(l1)+coord(l2)-pi*fpint(l2)
 495    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
        l1 = 1
        l2 = nuu+nrr
        if(nuest.lt.nu+1) l1=nuu+1
        if(nvest.lt.nv+2) l2=nuu
c  test whether we cannot further increase the number of knots.
        if(l1.gt.l2) go to 950
 500    fpmax = 0.
        l = 0
        do 510 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l.eq.0) go to 930
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.nuu) go to 530
c  addition in the u-direction
        l4 = l+4
        fpint(l) = 0.
        fac1 = tu(l4)-arg
        fac2 = arg-tu(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        j = nu
        do 520 i=l4,nu
          tu(j+1) = tu(j)
          j = j-1
 520    continue
        tu(l4) = arg
        nu = nu+1
        go to 570
c  addition in the v-direction
 530    l4 = l+4-nuu
        fpint(l) = 0.
        fac1 = tv(l4)-arg
        fac2 = arg-tv(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tv(j+1) = tv(j)
          j = j-1
 550    continue
        tv(l4) = arg
        nv = nv+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tv(j) = tv(i)+pi
 560    continue
c  restart the computations with the new set of knots.
 570  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing bicubic spline.               c
c ******************************************************               c
c we have determined the number of knots and their position. we now    c
c compute the coefficients of the smoothing spline sp(u,v).            c
c the observation matrix a is extended by the rows of a matrix, expres-c
c sing that sp(u,v) must be a constant function in the variable        c
c v and a cubic polynomial in the variable u. the corresponding        c
c weights of these additional rows are set to 1/(p). iteratively       c
c we than have to determine the value of p such that f(p) = sum((w(i)* c
c (z(i)-sp(u(i),v(i))))**2)  be = s.                                   c
c we already know that the least-squares polynomial corresponds to p=0,c
c and that the least-squares bicubic spline corresponds to p=infin.    c
c the iteration process makes use of rational interpolation. since f(p)c
c is a convex and strictly decreasing function of p, it can be approx- c
c imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
c three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
c f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
c of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tu(l),l=5,...,nu-4.
 580  call fpdisc(tu,nu,5,bu,nuest)
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tv(l),l=5,...,nv-4.
      call fpdisc(tv,nv,5,bv,nvest)
c  initial value for p.
      p1 = 0.
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 590 i=1,ncof
        p = p+a(i,1)
 590  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband4 = iband+ipar1
      if(iband4.gt.ncof) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
      nuu = nu4-iopt3-1
c  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 630 i=1,ncof
          ff(i) = f(i)
          do 620 j=1,iband4
            q(i,j) = 0.
 620      continue
          do 630 j=1,iband
            q(i,j) = a(i,j)
 630    continue
c  extend the observation matrix with the rows of a matrix, expressing
c  that for u=constant sp(u,v) must be a constant function.
        do 720 i=5,nv4
          ii = i-4
          do 635 l=1,nvv
             row(l) = 0.
 635      continue
          ll = ii
          do 640  l=1,5
             if(ll.gt.nvv) ll=1
             row(ll) = row(ll)+bv(ii,l)
             ll = ll+1
 640      continue
          do 720 j=1,nuu
c  initialize the new row.
            do 645 l=1,iband
              h(l) = 0.
 645        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            if(j.gt.iopt2) go to 665
            if(j.eq.2) go to 655
            do 650 k=1,2
               cs(k) = 0.
               do 650 l=1,nvv
                  cs(k) = cs(k)+cosi(k,l)*row(l)
 650        continue
            h(1) = cs(1)
            h(2) = cs(2)
            jrot = 2
            go to 675
 655        do 660 k=3,5
               cs(k) = 0.
               do 660 l=1,nvv
                  cs(k) = cs(k)+cosi(k,l)*row(l)
 660        continue
            h(1) = cs(1)*ratio
            h(2) = cs(2)*ratio
            h(3) = cs(3)
            h(4) = cs(4)
            h(5) = cs(5)
            jrot = 2
            go to 675
 665        do 670 l=1,nvv
               h(l) = row(l)
 670        continue
            jrot = ipar1+1+(j-iopt2-1)*nvv
 675        do 677 l=1,iband
              h(l) = h(l)*pinv
 677        continue
            zi = 0.
c  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) if(i2) 720,720,690
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,zi,ff(irot))
              if(i2.eq.0) go to 720
c  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          do 700 l=1,i2
                h(l) = h(l+1)
 700          continue
              h(i2+1) = 0.
 710        continue
 720    continue
c  extend the observation matrix with the rows of a matrix expressing
c  that for v=constant. sp(u,v) must be a cubic polynomial.
        do 810 i=5,nu4
          ii = i-4
          do 810 j=1,nvv
c  initialize the new row
            do 730 l=1,iband4
              h(l) = 0.
 730        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l-1
               if(il.eq.nu4 .and. iopt3.ne.0) go to 760
               if(il.gt.iopt2+1) go to 750
               go to (735,740,745),il
 735           h(1) = bu(ii,l)
               j1 = j+1
               go to 760
 740           h(1) = h(1)+bu(ii,l)
               h(2) = bu(ii,l)*cosi(1,j)
               h(3) = bu(ii,l)*cosi(2,j)
               j1 = j+3
               go to 760
 745           h(1) = h(1)+bu(ii,l)
               h(2) = bu(ii,l)*cosi(1,j)*ratio
               h(3) = bu(ii,l)*cosi(2,j)*ratio
               h(4) = bu(ii,l)*cosi(3,j)
               h(5) = bu(ii,l)*cosi(4,j)
               h(6) = bu(ii,l)*cosi(5,j)
               j1 = j+6
               go to 760
 750           h(j1) = bu(ii,l)
               j1 = j1+nvv
 760        continue
            do 765 l=1,iband4
              h(l) = h(l)*pinv
 765        continue
            zi = 0.
            jrot = 1
            if(ii.gt.iopt2+1) jrot = ipar1+(ii-iopt2-2)*nvv+j
c  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) if(i2) 810,810,780
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,zi,ff(irot))
              if(i2.eq.0) go to 810
c  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          do 790 l=1,i2
                h(l) = h(l+1)
 790          continue
              h(i2+1) = 0.
 800        continue
 810    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
        dmax = 0.
        do 820 i=1,ncof
          if(q(i,1).le.dmax) go to 820
          dmax = q(i,1)
 820    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 830 i=1,ncof
          if(q(i,1).le.sigma) go to 840
 830    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,ncc)
        rank = ncof
        go to 845
c  in case of rank deficiency, find the minimum norm solution.
 840    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 845    do 850 i=1,ncof
           q(i,1) = q(i,1)/dmax
 850    continue
c  find the coefficients in the standard b-spline representation of
c  the polar spline.
        call fprppo(nu,nv,iopt2,iopt3,cosi,ratio,c,ff,ncoff)
c  compute f(p).
        fp = 0.
        do 890 num = 1,nreg
          num1 = num-1
          lu = num1/nvv
          lv = num1-lu*nvv
          jrot = lu*nv4+lv
          in = index(num)
 860      if(in.eq.0) go to 890
          store = 0.
          i1 = jrot
          do 880 i=1,4
            hui = spu(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hui*spv(in,j)*c(j1)
 870        continue
            i1 = i1+nv4
 880      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
c  test whether the approximation sp(u,v) is an acceptable solution
        fpms = fp-s
        if(abs(fpms).le.acc) go to 980
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 940
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 900
        if((f2-f3).gt.acc) go to 895
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2.lt.0.) ich3 = 1
 900    if(ich1.ne.0) go to 910
        if((f1-f2).gt.acc) go to 905
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 920
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 910    if(f2.ge.f1 .or. f2.le.f3) go to 945
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 920  continue
c  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = 3
      go to 990
 945  ier = 2
      go to 990
 950  ier = 1
      go to 990
 960  ier = -2
      go to 990
 970  ier = -1
      fp = 0.
 980  if(ncof.ne.rank) ier = -rank
 990  return
      end

      subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
c  subroutine fprank finds the minimum norm solution of a least-
c  squares problem in case of rank deficiency.
c
c  input parameters:
c    a : array, which contains the non-zero elements of the observation
c        matrix after triangularization by givens transformations.
c    f : array, which contains the transformed right hand side.
c    n : integer,wich contains the dimension of a.
c    m : integer, which denotes the bandwidth of a.
c  tol : real value, giving a threshold to determine the rank of a.
c
c  output parameters:
c    c : array, which contains the minimum norm solution.
c   sq : real value, giving the contribution of reducing the rank
c        to the sum of squared residuals.
c rank : integer, which contains the rank of matrix a.
c
c  ..scalar arguments..
      integer n,m,na,rank
      real tol,sq
c  ..array arguments..
      real a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
c  ..local scalars..
      integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
      real cos,fac,piv,sin,yi
      double precision store,stor1,stor2,stor3
c  ..function references..
      integer min0
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      m1 = m-1
c  the rank deficiency nl is considered to be the number of sufficient
c  small diagonal elements of a.
      nl = 0
      sq = 0.
      do 90 i=1,n
        if(a(i,1).gt.tol) go to 90
c  if a sufficient small diagonal element is found, we put it to
c  zero. the remainder of the row corresponding to that zero diagonal
c  element is then rotated into triangle by givens rotations .
c  the rank deficiency is increased by one.
        nl = nl+1
        if(i.eq.n) go to 90
        yi = f(i)
        do 10 j=1,m1
          h(j) = a(i,j+1)
  10    continue
        h(m) = 0.
        i1 = i+1
        do 60 ii=i1,n
          i2 = min0(n-ii,m1)
          piv = h(1)
          if(piv.eq.0.) go to 30
          call fpgivs(piv,a(ii,1),cos,sin)
          call fprota(cos,sin,yi,f(ii))
          if(i2.eq.0) go to 70
          do 20 j=1,i2
            j1 = j+1
            call fprota(cos,sin,h(j1),a(ii,j1))
            h(j) = h(j1)
  20      continue
          go to 50
  30      if(i2.eq.0) go to 70
          do 40 j=1,i2
            h(j) = h(j+1)
  40      continue
  50      h(i2+1) = 0.
  60    continue
c  add to the sum of squared residuals the contribution of deleting
c  the row with small diagonal element.
  70    sq = sq+yi**2
  90  continue
c  rank denotes the rank of a.
      rank = n-nl
c  let b denote the (rank*n) upper trapezoidal matrix which can be
c  obtained from the (n*n) upper triangular matrix a by deleting
c  the rows and interchanging the columns corresponding to a zero
c  diagonal element. if this matrix is factorized using givens
c  transformations as  b = (r) (u)  where
c    r is a (rank*rank) upper triangular matrix,
c    u is a (rank*n) orthonormal matrix
c  then the minimal least-squares solution c is given by c = b' v,
c  where v is the solution of the system  (r) (r)' v = g  and
c  g denotes the vector obtained from the old right hand side f, by
c  removing the elements corresponding to a zero diagonal element of a.
c  initialization.
      do 100 i=1,rank
        do 100 j=1,m
          aa(i,j) = 0.
 100  continue
c  form in aa the upper triangular matrix obtained from a by
c  removing rows and columns with zero diagonal elements. form in ff
c  the new right hand side by removing the elements of the old right
c  hand side corresponding to a deleted row.
      ii = 0
      do 120 i=1,n
        if(a(i,1).le.tol) go to 120
        ii = ii+1
        ff(ii) = f(i)
        aa(ii,1) = a(i,1)
        jj = ii
        kk = 1
        j = i
        j1 = min0(j-1,m1)
        if(j1.eq.0) go to 120
        do 110 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 110
          kk = kk+1
          jj = jj-1
          aa(jj,kk) = a(j,k+1)
 110    continue
 120  continue
c  form successively in h the columns of a with a zero diagonal element.
      ii = 0
      do 200 i=1,n
        ii = ii+1
        if(a(i,1).gt.tol) go to 200
        ii = ii-1
        if(ii.eq.0) go to 200
        jj = 1
        j = i
        j1 = min0(j-1,m1)
        do 130 k=1,j1
          j = j-1
          if(a(j,1).le.tol) go to 130
          h(jj) = a(j,k+1)
          jj = jj+1
 130    continue
        do 140 kk=jj,m
          h(kk) = 0.
 140    continue
c  rotate this column into aa by givens transformations.
        jj = ii
        do 190 i1=1,ii
          j1 = min0(jj-1,m1)
          piv = h(1)
          if(piv.ne.0.) go to 160
          if(j1.eq.0) go to 200
          do 150 j2=1,j1
            j3 = j2+1
            h(j2) = h(j3)
 150      continue
          go to 180
 160      call fpgivs(piv,aa(jj,1),cos,sin)
          if(j1.eq.0) go to 200
          kk = jj
          do 170 j2=1,j1
            j3 = j2+1
            kk = kk-1
            call fprota(cos,sin,h(j3),aa(kk,j3))
            h(j2) = h(j3)
 170      continue
 180      jj = jj-1
          h(j3) = 0.
 190    continue
 200  continue
c  solve the system (aa) (f1) = ff
      ff(rank) = ff(rank)/aa(rank,1)
      i = rank-1
      if(i.eq.0) go to 230
      do 220 j=2,rank
        store = ff(i)
        i1 = min0(j-1,m1)
        k = i
        do 210 ii=1,i1
          k = k+1
          stor1 = ff(k)
          stor2 = aa(i,ii+1)
          store = store-stor1*stor2
 210    continue
        stor1 = aa(i,1)
        ff(i) = store/stor1
        i = i-1
 220  continue
c  solve the system  (aa)' (f2) = f1
 230  ff(1) = ff(1)/aa(1,1)
      if(rank.eq.1) go to 260
      do 250 j=2,rank
        store = ff(j)
        i1 = min0(j-1,m1)
        k = j
        do 240 ii=1,i1
          k = k-1
          stor1 = ff(k)
          stor2 = aa(k,ii+1)
          store = store-stor1*stor2
 240    continue
        stor1 = aa(j,1)
        ff(j) = store/stor1
 250  continue
c  premultiply f2 by the transpoze of a.
 260  k = 0
      do 280 i=1,n
        store = 0.
        if(a(i,1).gt.tol) k = k+1
        j1 = min0(i,m)
        kk = k
        ij = i+1
        do 270 j=1,j1
          ij = ij-1
          if(a(ij,1).le.tol) go to 270
          stor1 = a(ij,j)
          stor2 = ff(kk)
          store = store+stor1*stor2
          kk = kk-1
 270    continue
        c(i) = store
 280  continue
c  add to the sum of squared residuals the contribution of putting
c  to zero the small diagonal elements of matrix (a).
      stor3 = 0.
      do 310 i=1,n
        if(a(i,1).gt.tol) go to 310
        store = f(i)
        i1 = min0(n-i,m1)
        if(i1.eq.0) go to 300
        do 290 j=1,i1
          ij = i+j
          stor1 = c(ij)
          stor2 = a(i,j+1)
          store = store-stor1*stor2
 290    continue
 300    fac = a(i,1)*c(i)
        stor1 = a(i,1)
        stor2 = c(i)
        stor1 = stor1*stor2
        stor3 = stor3+stor1*(stor1-store-store)
 310  continue
      fac = stor3
      sq = sq+fac
      return
      end

      real function fprati(p1,f1,p2,f2,p3,f3)
c  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
c  gives the value of p such that the rational interpolating function
c  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
c  ..
c  ..scalar arguments..
      real p1,f1,p2,f2,p3,f3
c  ..local scalars..
      real h1,h2,h3,p
c  ..
      if(p3.gt.0.) go to 10
c  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
c  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
c  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
      subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,
     * nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,
     * reducy,fpintx,fpinty,lastdi,nplusx,nplusy,nrx,nry,nrdatx,nrdaty,
     * wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      real xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy
      integer iopt,mx,my,mz,kx,ky,nxest,nyest,maxit,nc,nx,ny,lastdi,
     * nplusx,nplusy,lwrk,ier
c  ..array arguments..
      real x(mx),y(my),z(mz),tx(nxest),ty(nyest),c(nc),fpintx(nxest),
     * fpinty(nyest),wrk(lwrk)
      integer nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
c  ..local scalars
      real acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,half,con1,con9,con4
      integer i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2,
     * k3,l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,ncof,
     * nk1x,nk1y,nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx,
     * nrinty,nxe,nxk,nye
c  ..function references..
      real abs,fprati
      integer max0,min0
c  ..subroutine references..
c    fpgrre,fpknot
c  ..
c   set constants
      one = 1
      half = 0.5e0
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c  we partition the working space.
      kx1 = kx+1
      ky1 = ky+1
      kx2 = kx1+1
      ky2 = ky1+1
      lsx = 1
      lsy = lsx+mx*kx1
      lri = lsy+my*ky1
      mm = max0(nxest,my)
      lq = lri+mm
      mynx = nxest*my
      lax = lq+mynx
      nxk = nxest*kx2
      lbx = lax+nxk
      lay = lbx+nxk
      lby = lay+nyest*ky2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(x,y), c
c  and the corresponding sum of squared residuals fp=f(p=inf).         c
c  if iopt=-1  sinf(x,y) is the requested approximation.               c
c  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
c    if fp <=s we will continue with the current set of knots.         c
c    if fp > s we will increase the number of knots and compute the    c
c       corresponding least-squares spline until finally fp<=s.        c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
c    if s>0 and                                                        c
c     *iopt=0 we first compute the least-squares polynomial of degree  c
c      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
c     *iopt=1 we start with the knots found at the last call of the    c
c      routine, except for the case that s > fp0; then we can compute  c
c      the least-squares polynomial directly.                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  determine the number of knots for polynomial approximation.
      nminx = 2*kx1
      nminy = 2*ky1
      if(iopt.lt.0) go to 120
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  find nmaxx and nmaxy which denote the number of knots in x- and y-
c  direction in case of spline interpolation.
      nmaxx = mx+kx1
      nmaxy = my+ky1
c  find nxe and nye which denote the maximum number of knots
c  allowed in each direction
      nxe = min0(nmaxx,nxest)
      nye = min0(nmaxy,nyest)
      if(s.gt.0.) go to 100
c  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
c  test whether the required storage space exceeds the available one.
      if(ny.gt.nyest .or. nx.gt.nxest) go to 420
c  find the position of the interior knots in case of interpolation.
c  the knots in the x-direction.
      mk1 = mx-kx1
      if(mk1.eq.0) go to 60
      k3 = kx/2
      i = kx1+1
      j = k3+2
      if(k3*2.eq.kx) go to 40
      do 30 l=1,mk1
        tx(i) = x(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
  40  do 50 l=1,mk1
        tx(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  50  continue
c  the knots in the y-direction.
  60  mk1 = my-ky1
      if(mk1.eq.0) go to 120
      k3 = ky/2
      i = ky1+1
      j = k3+2
      if(k3*2.eq.ky) go to 80
      do 70 l=1,mk1
        ty(i) = y(j)
        i = i+1
        j = j+1
  70  continue
      go to 120
  80  do 90 l=1,mk1
        ty(i) = (y(j)+y(j-1))*half
        i = i+1
        j = j+1
  90  continue
      go to 120
c  if s > 0 our initial choice of knots depends on the value of iopt.
 100  if(iopt.eq.0) go to 115
      if(fp0.le.s) go to 115
c  if iopt=1 and fp0 > s we start computing the least- squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates x(i) inside each knot
c  interval (tx(l),tx(l+1)).
      l = kx2
      j = 1
      nrdatx(1) = 0
      mpm = mx-1
      do 105 i=2,mpm
        nrdatx(j) = nrdatx(j)+1
        if(x(i).lt.tx(l)) go to 105
        nrdatx(j) = nrdatx(j)-1
        l = l+1
        j = j+1
        nrdatx(j) = 0
 105  continue
c  we determine the number of grid coordinates y(i) inside each knot
c  interval (ty(l),ty(l+1)).
      l = ky2
      j = 1
      nrdaty(1) = 0
      mpm = my-1
      do 110 i=2,mpm
        nrdaty(j) = nrdaty(j)+1
        if(y(i).lt.ty(l)) go to 110
        nrdaty(j) = nrdaty(j)-1
        l = l+1
        j = j+1
        nrdaty(j) = 0
 110  continue
      go to 120
c  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
c  polynomial of degree kx in x and ky in y (which is a spline without
c  interior knots).
 115  nx = nminx
      ny = nminy
      nrdatx(1) = mx-2
      nrdaty(1) = my-2
      lastdi = 0
      nplusx = 0
      nplusy = 0
      fp0 = 0.
      fpold = 0.
      reducx = 0.
      reducy = 0.
 120  mpm = mx+my
      ifsx = 0
      ifsy = 0
      ifbx = 0
      ifby = 0
      p = -one
c  main loop for the different sets of knots.mpm=mx+my is a save upper
c  bound for the number of trials.
      do 250 iter=1,mpm
        if(nx.eq.nminx .and. ny.eq.nminy) ier = -2
c  find nrintx (nrinty) which is the number of knot intervals in the
c  x-direction (y-direction).
        nrintx = nx-nminx+1
        nrinty = ny-nminy+1
c  find ncof, the number of b-spline coefficients for the current set
c  of knots.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        i = nx
        do 130 j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
 130    continue
        i = ny
        do 140 j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
 140    continue
c  find the least-squares spline sinf(x,y) and calculate for each knot
c  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
c  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
c  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
c  value belonging to that interval.
c  fp gives the total sum of squared residuals.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,
     *  ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),
     *  wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),
     *  nrx,nry)
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430
c  increase the number of knots.
c  if nx=nxe and ny=nye we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nx.eq.nxe .and. ny.eq.nye) go to 420
        ier = 0
c  adjust the parameter reducx or reducy according to the direction
c  in which the last added knots were located.
        if(lastdi) 150,170,160
 150    reducx = fpold-fp
        go to 170
 160    reducy = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
c  find nplx, the number of knots we should add in the x-direction.
        nplx = 1
        if(nx.eq.nminx) go to 180
        npl1 = nplusx*2
        rn = nplusx
        if(reducx.gt.acc) npl1 = rn*fpms/reducx
        nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
c  find nply, the number of knots we should add in the y-direction.
 180    nply = 1
        if(ny.eq.nminy) go to 190
        npl1 = nplusy*2
        rn = nplusy
        if(reducy.gt.acc) npl1 = rn*fpms/reducy
        nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
 190    if(nplx-nply) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nx.eq.nxe) go to 230
c  addition in the x-direction.
        lastdi = -1
        nplusx = nplx
        ifsx = 0
        do 220 l=1,nplusx
c  add a new knot in the x-direction
          call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
c  test whether we cannot further increase the number of knots in the
c  x-direction.
          if(nx.eq.nxe) go to 250
 220    continue
        go to 250
 230    if(ny.eq.nye) go to 210
c  addition in the y-direction.
        lastdi = 1
        nplusy = nply
        ifsy = 0
        do 240 l=1,nplusy
c  add a new knot in the y-direction.
          call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
c  test whether we cannot further increase the number of knots in the
c  y-direction.
          if(ny.eq.nye) go to 250
 240    continue
c  restart the computations with the new set of knots.
 250  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
c  this smoothing spline varies with the parameter p in such a way thatc
c    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. iteratively we then have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(x,y) and the corresponding sum of
c  squared residuals fp.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,
     *  ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),
     *  wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),
     *  nrx,nry)
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end

      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      real cos,sin,a,b
c ..local scalars..
      real stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
      subroutine fprppo(nu,nv,if1,if2,cosi,ratio,c,f,ncoff)
c  given the coefficients of a constrained bicubic spline, as determined
c  in subroutine fppola, subroutine fprppo calculates the coefficients
c  in the standard b-spline representation of bicubic splines.
c  ..
c  ..scalar arguments..
      real ratio
      integer nu,nv,if1,if2,ncoff
c  ..array arguments
      real c(ncoff),f(ncoff),cosi(5,nv)
c  ..local scalars..
      integer i,iopt,ii,j,k,l,nu4,nvv
c  ..
      nu4 = nu-4
      nvv = nv-7
      iopt = if1+1
      do 10 i=1,ncoff
         f(i) = 0.
  10  continue
      i = 0
      do 120 l=1,nu4
         ii = i
         if(l.gt.iopt) go to 80
         go to (20,40,60),l
  20     do 30 k=1,nvv
            i = i+1
            f(i) = c(1)
  30     continue
         j = 1
         go to 100
  40     do 50 k=1,nvv
            i = i+1
            f(i) = c(1)+c(2)*cosi(1,k)+c(3)*cosi(2,k)
  50     continue
         j = 3
         go to 100
  60     do 70 k=1,nvv
            i = i+1
            f(i) = c(1)+ratio*(c(2)*cosi(1,k)+c(3)*cosi(2,k))+
     *             c(4)*cosi(3,k)+c(5)*cosi(4,k)+c(6)*cosi(5,k)
  70     continue
         j = 6
         go to 100
  80     if(l.eq.nu4 .and. if2.ne.0) go to 120
         do 90 k=1,nvv
            i = i+1
            j = j+1
            f(i) = c(j)
  90     continue
 100     do 110 k=1,3
            ii = ii+1
            i = i+1
            f(i) = f(ii)
 110     continue
 120  continue
      do 130 i=1,ncoff
         c(i) = f(i)
 130  continue
      return
      end

      subroutine fprpsp(nt,np,co,si,c,f,ncoff)
c  given the coefficients of a spherical spline function, subroutine
c  fprpsp calculates the coefficients in the standard b-spline re-
c  presentation of this bicubic spline.
c  ..
c  ..scalar arguments
      integer nt,np,ncoff
c  ..array arguments
      real co(np),si(np),c(ncoff),f(ncoff)
c  ..local scalars
      real cn,c1,c2,c3
      integer i,ii,j,k,l,ncof,npp,np4,nt4
c  ..
      nt4 = nt-4
      np4 = np-4
      npp = np4-3
      ncof = 6+npp*(nt4-4)
      c1 = c(1)
      cn = c(ncof)
      j = ncoff
      do 10 i=1,np4
         f(i) = c1
         f(j) = cn
         j = j-1
  10  continue
      i = np4
      j=1
      do 70 l=3,nt4
         ii = i
         if(l.eq.3 .or. l.eq.nt4) go to 30
         do 20 k=1,npp
            i = i+1
            j = j+1
            f(i) = c(j)
  20     continue
         go to 50
  30     if(l.eq.nt4) c1 = cn
         c2 = c(j+1)
         c3 = c(j+2)
         j = j+2
         do 40 k=1,npp
            i = i+1
            f(i) = c1+c2*co(k)+c3*si(k)
  40     continue
  50     do 60 k=1,3
            ii = ii+1
            i = i+1
            f(i) = f(ii)
  60     continue
  70  continue
      do 80 i=1,ncoff
         c(i) = f(i)
  80  continue
      return
      end
      subroutine fpseno(maxtr,up,left,right,info,merk,ibind,nbind)
c  subroutine fpseno fetches a branch of a triply linked tree the
c  information of which is kept in the arrays up,left,right and info.
c  the branch has a specified length nbind and is determined by the
c  parameter merk which points to its terminal node. the information
c  field of the nodes of this branch is stored in the array ibind. on
c  exit merk points to a new branch of length nbind or takes the value
c  1 if no such branch was found.
c  ..
c  ..scalar arguments..
      integer maxtr,merk,nbind
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),
     * ibind(nbind)
c  ..scalar arguments..
      integer i,j,k
c  ..
      k = merk
      j = nbind
      do 10 i=1,nbind
        ibind(j) = info(k)
        k = up(k)
        j = j-1
  10  continue
  20  k = right(merk)
      if(k.ne.0) go to 30
      merk = up(merk)
      if(merk-1) 40,40,20
  30  merk = k
      k = left(merk)
      if(k.ne.0) go to 30
  40  return
      end
      subroutine fpspgr(iopt,ider,u,mu,v,mv,r,mr,r0,r1,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,
     * fpintv,dr,step,lastdi,nplusu,nplusv,lastu0,lastu1,nru,nrv,
     * nrdatu,nrdatv,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      integer mu,mv,mr,nuest,nvest,maxit,nc,nu,nv,lastdi,nplusu,nplusv,
     * lastu0,lastu1,lwrk,ier
      real r0,r1,s,tol,fp,fp0,fpold,reducu,reducv
c  ..array arguments..
      integer iopt(3),ider(4),nrdatu(nuest),nrdatv(nvest),nru(mu),
     * nrv(mv)
      real u(mu),v(mv),r(mr),tu(nuest),tv(nvest),c(nc),fpintu(nuest),
     * fpintv(nvest),dr(6),wrk(lwrk),step(2)
c  ..local scalars..
      real acc,fpms,f1,f2,f3,p,per,pi,p1,p2,p3,vb,ve,rmax,rmin,rn,one,
     * con1,con4,con9
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,istart,iter,i1,i2,j,ju,
     * ktu,l,l1,l2,l3,l4,mpm,mumin,mu0,mu1,nn,nplu,nplv,npl1,nrintu,
     * nrintv,nue,numax,nve,nvmax
c  ..local arrays..
      integer idd(4)
      real drr(6)
c  ..function references..
      real abs,atan2,fprati
      integer max0,min0
c  ..subroutine references..
c    fpknot,fpopsp
c  ..
c   set constants
      one = 1
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
c   initialization
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
      mumin = 4
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      if(ider(3).ge.0) mumin = mumin-1
      if(iopt(3).eq.1 .and. ider(4).eq.1) mumin = mumin-1
      if(mumin.eq.0) mumin = 1
      pi = atan2(0.,-one)
      per = pi+pi
      vb = v(1)
      ve = vb+per
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c  given a set of knots we compute the least-squares spline sinf(u,v)  c
c  and the corresponding sum of squared residuals fp = f(p=inf).       c
c  if iopt(1)=-1  sinf(u,v) is the requested approximation.            c
c  if iopt(1)>=0  we check whether we can accept the knots:            c
c    if fp <= s we will continue with the current set of knots.        c
c    if fp >  s we will increase the number of knots and compute the   c
c       corresponding least-squares spline until finally fp <= s.      c
c    the initial choice of knots depends on the value of s and iopt.   c
c    if s=0 we have spline interpolation; in that case the number of   c
c     knots in the u-direction equals nu=numax=mu+6+iopt(2)+iopt(3)    c
c     and in the v-direction nv=nvmax=mv+7.                            c
c    if s>0 and                                                        c
c      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  c
c       spline without interior knots : nu=8 ; nv=8.                   c
c      iopt(1)=1 we start with the set of knots found at the last call c
c       of the routine, except for the case that s > fp0; then we      c
c       compute the least-squares polynomial directly.                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(iopt(1).lt.0) go to 120
c  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
c  numax and nvmax denote the number of knots needed for interpolation.
      numax = mu+6+iopt(2)+iopt(3)
      nvmax = mv+7
      nue = min0(numax,nuest)
      nve = min0(nvmax,nvest)
      if(s.gt.0.) go to 100
c  if s = 0, s(u,v) is an interpolating spline.
      nu = numax
      nv = nvmax
c  test whether the required storage space exceeds the available one.
      if(nu.gt.nuest .or. nv.gt.nvest) go to 420
c  find the position of the knots in the v-direction.
      do 10 l=1,mv
        tv(l+3) = v(l)
  10  continue
      tv(mv+4) = ve
      l1 = mv-2
      l2 = mv+5
      do 20 i=1,3
         tv(i) = v(l1)-per
         tv(l2) = v(i+1)+per
         l1 = l1+1
         l2 = l2+1
  20  continue
c  if not all the derivative values g(i,j) are given, we will first
c  estimate these values by computing a least-squares spline
      idd(1) = ider(1)
      if(idd(1).eq.0) idd(1) = 1
      if(idd(1).gt.0) dr(1) = r0
      idd(2) = ider(2)
      idd(3) = ider(3)
      if(idd(3).eq.0) idd(3) = 1
      if(idd(3).gt.0) dr(4) = r1
      idd(4) = ider(4)
      if(ider(1).lt.0 .or. ider(3).lt.0) go to 30
      if(iopt(2).ne.0 .and. ider(2).eq.0) go to 30
      if(iopt(3).eq.0 .or. ider(4).ne.0) go to 70
c we set up the knots in the u-direction for computing the least-squares
c spline.
  30  i1 = 3
      i2 = mu-2
      nu = 4
      do 40 i=1,mu
         if(i1.gt.i2) go to 50
         nu = nu+1
         tu(nu) = u(i1)
         i1 = i1+2
  40  continue
  50  do 60 i=1,4
         tu(i) = 0.
         nu = nu+1
         tu(nu) = pi
  60  continue
c we compute the least-squares spline for estimating the derivatives.
      call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      ifsu = 0
c if all the derivatives at the origin are known, we compute the
c interpolating spline.
c we set up the knots in the u-direction, needed for interpolation.
  70  nn = numax-8
      if(nn.eq.0) go to 95
      ju = 2-iopt(2)
      do 80 l=1,nn
        tu(l+4) = u(ju)
        ju = ju+1
  80  continue
      nu = numax
      l = nu
      do 90 i=1,4
         tu(i) = 0.
         tu(l) = pi
         l = l-1
  90  continue
c we compute the interpolating spline.
  95  call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,idd,
     *  tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,nrv,
     *  wrk,lwrk)
      go to 430
c  if s>0 our initial choice of knots depends on the value of iopt(1).
 100  ier = 0
      if(iopt(1).eq.0) go to 115
      step(1) = -step(1)
      step(2) = -step(2)
      if(fp0.le.s) go to 115
c  if iopt(1)=1 and fp0 > s we start computing the least-squares spline
c  according to the set of knots found at the last call of the routine.
c  we determine the number of grid coordinates u(i) inside each knot
c  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mu0 = 2-iopt(2)
      mu1 = mu-1+iopt(3)
      do 105 i=mu0,mu1
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 105
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
 105  continue
c  we determine the number of grid coordinates v(i) inside each knot
c  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      do 110 i=2,mv
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 110
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
 110  continue
      idd(1) = ider(1)
      idd(2) = ider(2)
      idd(3) = ider(3)
      idd(4) = ider(4)
      go to 120
c  if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
c  squares polynomial (which is a spline without interior knots).
 115  ier = -2
      idd(1) = ider(1)
      idd(2) = 1
      idd(3) = ider(3)
      idd(4) = 1
      nu = 8
      nv = 8
      nrdatu(1) = mu-2+iopt(2)+iopt(3)
      nrdatv(1) = mv-1
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
c  main loop for the different sets of knots.mpm=mu+mv is a save upper
c  bound for the number of trials.
 120  mpm = mu+mv
      do 270 iter=1,mpm
c  find nrintu (nrintv) which is the number of knot intervals in the
c  u-direction (v-direction).
        nrintu = nu-7
        nrintv = nv-7
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(u,v).
        i = nu
        do 125 j=1,4
          tu(j) = 0.
          tu(i) = pi
          i = i-1
 125    continue
        l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 130 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-per
          tv(l3) = tv(l1)+per
 130    continue
c  find an estimate of the range of possible values for the optimal
c  derivatives at the origin.
        ktu = nrdatu(1)+2-iopt(2)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu0) go to 140
         rmin = r0
         rmax = r0
         l = mv*ktu
         do 135 i=1,l
            if(r(i).lt.rmin) rmin = r(i)
            if(r(i).gt.rmax) rmax = r(i)
 135     continue
         step(1) = rmax-rmin
         lastu0 = ktu
 140    ktu = nrdatu(nrintu)+2-iopt(3)
        if(ktu.lt.mumin) ktu = mumin
        if(ktu.eq.lastu1) go to 150
         rmin = r1
         rmax = r1
         l = mv*ktu
         j = mr
         do 145 i=1,l
            if(r(j).lt.rmin) rmin = r(j)
            if(r(j).gt.rmax) rmax = r(j)
            j = j-1
 145     continue
         step(2) = rmax-rmin
         lastu1 = ktu
c  find the least-squares spline sinf(u,v).
 150    call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,dr,iopt,
     *   idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,
     *   nrv,wrk,lwrk)
        if(step(1).lt.0.) step(1) = -step(1)
        if(step(2).lt.0.) step(2) = -step(2)
        if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt(1).lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
c  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
c  if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if(nu.eq.numax .and. nv.eq.nvmax) go to 430
c  increase the number of knots.
c  if nu=nue and nv=nve we cannot further increase the number of knots
c  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        if(ider(1).eq.0) fpintu(1) = fpintu(1)+(r0-dr(1))**2
        if(ider(3).eq.0) fpintu(nrintu) = fpintu(nrintu)+(r1-dr(4))**2
        ier = 0
c  adjust the parameter reducu or reducv according to the direction
c  in which the last added knots were located.
        if(lastdi) 160,155,170
 155     nplv = 3
         idd(2) = ider(2)
         idd(4) = ider(4)
         fpold = fp
         go to 230
 160    reducu = fpold-fp
        go to 175
 170    reducv = fpold-fp
c  store the sum of squared residuals for the current set of knots.
 175    fpold = fp
c  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.8) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
c  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 3
        if(nv.eq.8) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
c  test whether we are going to add knots in the u- or v-direction.
 190    if(nplu-nplv) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
c  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        istart = 0
        if(iopt(2).eq.0) istart = 1
        do 220 l=1,nplusu
c  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,istart)
c  test whether we cannot further increase the number of knots in the
c  u-direction.
          if(nu.eq.nue) go to 270
 220    continue
        go to 270
 230    if(nv.eq.nve) go to 210
c  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
c  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
c  test whether we cannot further increase the number of knots in the
c  v-direction.
          if(nv.eq.nve) go to 270
 240    continue
c  restart the computations with the new set of knots.
 270  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 300  if(ier.eq.(-2)) go to 440
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(u,v)                c
c *****************************************************                c
c  we have determined the number of knots and their position. we now   c
c  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
c  this smoothing spline depends on the parameter p in such a way that c
c    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             c
c  is a continuous, strictly decreasing function of p. moreover the    c
c  least-squares polynomial corresponds to p=0 and the least-squares   c
c  spline to p=infinity. then iteratively we have to determine the     c
c  positive value of p such that f(p)=s. the process which is proposed c
c  here makes use of rational interpolation. f(p) is approximated by a c
c  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
c  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
c  are used to calculate the new value of p such that r(p)=s.          c
c  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      do 305 i=1,6
        drr(i) = dr(i)
 305  continue
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
c  find the smoothing spline sp(u,v) and the corresponding sum f(p).
        call fpopsp(ifsu,ifsv,ifbu,ifbv,u,mu,v,mv,r,mr,r0,r1,drr,iopt,
     *   idd,tu,nu,tv,nv,nuest,nvest,p,step,c,nc,fp,fpintu,fpintv,nru,
     *   nrv,wrk,lwrk)
c  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 400
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
c  test whether the iteration process proceeds as theoretically
c  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
c  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
      subroutine fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eta,tol,maxit,
     * ib1,ib3,nc,ncc,intest,nrest,nt,tt,np,tp,c,fp,sup,fpint,coord,f,
     * ff,row,coco,cosi,a,q,bt,bp,spt,spp,h,index,nummer,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      integer iopt,m,ntest,npest,maxit,ib1,ib3,nc,ncc,intest,nrest,
     * nt,np,lwrk,ier
      real s,eta,tol,fp,sup
c  ..array arguments..
      real teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),c(nc),
     * fpint(intest),coord(intest),f(ncc),ff(nc),row(npest),coco(npest),
     * cosi(npest),a(ncc,ib1),q(ncc,ib3),bt(ntest,5),bp(npest,5),
     * spt(m,4),spp(m,4),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
c  ..local scalars..
      real aa,acc,arg,cn,co,c1,dmax,d1,d2,eps,facc,facs,fac1,fac2,fn,
     * fpmax,fpms,f1,f2,f3,hti,htj,p,pi,pinv,piv,pi2,p1,p2,p3,ri,si,
     * sigma,sq,store,wi,rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ich1,ich3,ii,ij,il,in,irot,
     * iter,i1,i2,i3,j,jlt,jrot,j1,j2,l,la,lf,lh,ll,lp,lt,lwest,l1,l2,
     * l3,l4,ncof,ncoff,npp,np4,nreg,nrint,nrr,nr1,ntt,nt4,nt6,num,
     * num1,rank
c  ..local arrays..
      real ht(4),hp(4)
c  ..function references..
      real abs,atan,fprati,sqrt,cos,sin
      integer min0
c  ..subroutine references..
c   fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota,fprpsp
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
      pi = atan(one)*4
      pi2 = pi+pi
      eps = sqrt(eta)
      if(iopt.lt.0) go to 70
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(s.lt.sup) if(np-11) 60,70,70
c  if iopt=0 we begin by computing the weighted least-squares polynomial
c  of the form
c     s(teta,phi) = c1*f1(teta) + cn*fn(teta)
c  where f1(teta) and fn(teta) are the cubic polynomials satisfying
c     f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
c  the corresponding weighted sum of squared residuals gives the upper
c  bound sup for the smoothing factor s.
  10  sup = 0.
      d1 = 0.
      d2 = 0.
      c1 = 0.
      cn = 0.
      fac1 = pi*(one + half)
      fac2 = (one + one)/pi**3
      aa = 0.
      do 40 i=1,m
         wi = w(i)
         ri = r(i)*wi
         arg = teta(i)
         fn = fac2*arg*arg*(fac1-arg)
         f1 = (one-fn)*wi
         fn = fn*wi
         if(fn.eq.0.) go to 20
         call fpgivs(fn,d1,co,si)
         call fprota(co,si,f1,aa)
         call fprota(co,si,ri,cn)
 20      if(f1.eq.0.) go to 30
         call fpgivs(f1,d2,co,si)
         call fprota(co,si,ri,c1)
 30      sup = sup+ri*ri
 40   continue
      if(d2.ne.0.) c1 = c1/d2
      if(d1.ne.0.) cn = (cn-aa*c1)/d1
c  find the b-spline representation of this least-squares polynomial
      nt = 8
      np = 8
      do 50 i=1,4
         c(i) = c1
         c(i+4) = c1
         c(i+8) = cn
         c(i+12) = cn
         tt(i) = 0.
         tt(i+4) = pi
         tp(i) = 0.
         tp(i+4) = pi2
  50  continue
      fp = sup
c  test whether the least-squares polynomial is an acceptable solution
      fpms = sup-s
      if(fpms.lt.acc) go to 960
c  test whether we cannot further increase the number of knots.
  60  if(npest.lt.11 .or. ntest.lt.9) go to 950
c  find the initial set of interior knots of the spherical spline in
c  case iopt = 0.
      np = 11
      tp(5) = pi*half
      tp(6) = pi
      tp(7) = tp(5)+pi
      nt = 9
      tt(5) = tp(5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  part 1 : computation of least-squares spherical splines.            c
c  ********************************************************            c
c  if iopt < 0 we compute the least-squares spherical spline according c
c  to the given set of knots.                                          c
c  if iopt >=0 we compute least-squares spherical splines with increas-c
c  ing numbers of knots until the corresponding sum f(p=inf)<=s.       c
c  the initial set of knots then depends on the value of iopt:         c
c    if iopt=0 we start with one interior knot in the teta-direction   c
c              (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). c
c    if iopt>0 we start with the set of knots found at the last call   c
c              of the routine.                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  70  do 570 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(teta,phi).
         l1 = 4
         l2 = l1
         l3 = np-3
         l4 = l3
         tp(l2) = 0.
         tp(l3) = pi2
         do 80 i=1,3
            l1 = l1+1
            l2 = l2-1
            l3 = l3+1
            l4 = l4-1
            tp(l2) = tp(l4)-pi2
            tp(l3) = tp(l1)+pi2
  80     continue
        l = nt
        do 90 i=1,4
          tt(i) = 0.
          tt(l) = pi
          l = l-1
  90    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        ntt = nt-7
        npp = np-7
        nrr = npp/2
        nr1 = nrr+1
        nrint = ntt+npp
        nreg = ntt*npp
c  arrange the data points according to the panel they belong to.
        call fporde(teta,phi,m,3,3,tt,nt,tp,np,nummer,index,nreg)
c  find the b-spline coefficients coco and cosi of the cubic spline
c  approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
        do 100 i=1,npp
           coco(i) = 0.
           cosi(i) = 0.
           do 100 j=1,npp
              a(i,j) = 0.
 100    continue
c  the coefficients coco and cosi are obtained from the conditions
c  sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
        do 150 i=1,npp
           l2 = i+3
           arg = tp(l2)
           call fpbspl(tp,np,3,arg,l2,hp)
           do 110 j=1,npp
              row(j) = 0.
 110       continue
           ll = i
           do 120 j=1,3
              if(ll.gt.npp) ll= 1
              row(ll) = row(ll)+hp(j)
              ll = ll+1
 120       continue
           facc = cos(arg)
           facs = sin(arg)
           do 140 j=1,npp
              piv = row(j)
              if(piv.eq.0.) go to 140
              call fpgivs(piv,a(j,1),co,si)
              call fprota(co,si,facc,coco(j))
              call fprota(co,si,facs,cosi(j))
              if(j.eq.npp) go to 150
              j1 = j+1
              i2 = 1
              do 130 l=j1,npp
                 i2 = i2+1
                 call fprota(co,si,row(l),a(j,i2))
 130          continue
 140       continue
 150    continue
        call fpback(a,coco,npp,npp,coco,ncc)
        call fpback(a,cosi,npp,npp,cosi,ncc)
c  find ncof, the dimension of the spherical spline and ncoff, the
c  number of coefficients in the standard b-spline representation.
        nt4 = nt-4
        np4 = np-4
        ncoff = nt4*np4
        ncof = 6+npp*(ntt-1)
c  find the bandwidth of the observation matrix a.
        iband = 4*npp
        if(ntt.eq.4) iband = 3*(npp+1)
        if(ntt.lt.4) iband = ncof
        iband1 = iband-1
c  initialize the observation matrix a.
        do 160 i=1,ncof
          f(i) = 0.
          do 160 j=1,iband
            a(i,j) = 0.
 160    continue
c  initialize the sum of squared residuals.
        fp = 0.
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 340 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lt = num1/npp
          l1 = lt+4
          lp = num1-lt*npp+1
          l2 = lp+3
          lt = lt+1
          jrot = 0
          if(lt.gt.2) jrot = 3+(lt-3)*npp
c  test whether there are still data points in the current panel.
          in = index(num)
 170      if(in.eq.0) go to 340
c  fetch a new data point.
          wi = w(in)
          ri = r(in)*wi
c  evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
          call fpbspl(tt,nt,3,teta(in),l1,ht)
c  evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
          call fpbspl(tp,np,3,phi(in),l2,hp)
c  store the value of these b-splines in spt and spp resp.
          do 180 i=1,4
            spp(in,i) = hp(i)
            spt(in,i) = ht(i)
 180      continue
c  initialize the new row of observation matrix.
          do 190 i=1,iband
            h(i) = 0.
 190      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in teta- and phi-direction and
c  by taking into account the conditions of the spherical splines.
          do 200 i=1,npp
             row(i) = 0.
 200      continue
c  take into account the condition (3) of the spherical splines.
          ll = lp
          do 210 i=1,4
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+hp(i)
             ll = ll+1
 210      continue
c  take into account the other conditions of the spherical splines.
          if(lt.gt.2 .and. lt.lt.(ntt-1)) go to 230
          facc = 0.
          facs = 0.
          do 220 i=1,npp
             facc = facc+row(i)*coco(i)
             facs = facs+row(i)*cosi(i)
 220     continue
c  fill in the non-zero elements of the new row.
 230     j1 = 0
         do 280 j =1,4
            jlt = j+lt
            htj = ht(j)
            if(jlt.gt.2 .and. jlt.le.nt4) go to 240
            j1 = j1+1
            h(j1) = h(j1)+htj
            go to 280
 240        if(jlt.eq.3 .or. jlt.eq.nt4) go to 260
            do 250 i=1,npp
               j1 = j1+1
               h(j1) = row(i)*htj
 250        continue
            go to 280
 260        if(jlt.eq.3) go to 270
            h(j1+1) = facc*htj
            h(j1+2) = facs*htj
            h(j1+3) = htj
            j1 = j1+2
            go to 280
 270        h(1) = h(1)+htj
            h(2) = facc*htj
            h(3) = facs*htj
            j1 = 3
 280      continue
          do 290 i=1,iband
            h(i) = h(i)*wi
 290      continue
c  rotate the row into triangle by givens transformations.
          irot = jrot
          do 310 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 310
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the right hand side.
            call fprota(co,si,ri,f(irot))
            if(i.eq.iband) go to 320
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 300 j=i3,iband
              i2 = i2+1
              call fprota(co,si,h(j),a(irot,i2))
 300        continue
 310      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 320      fp = fp+ri**2
c  find the number of the next data point in the panel.
 330      in = nummer(in)
          go to 170
 340    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 350 i=1,ncof
          if(a(i,1).le.dmax) go to 350
          dmax = a(i,1)
 350    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 360 i=1,ncof
          if(a(i,1).le.sigma) go to 370
 360    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,ncc)
        rank = ncof
        do 365 i=1,ncof
          q(i,1) = a(i,1)/dmax
 365    continue
        go to 390
c  in case of rank deficiency, find the minimum norm solution.
 370    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband
        do 380 i=1,ncof
          ff(i) = f(i)
          do 380 j=1,iband
            q(i,j) = a(i,j)
 380    continue
        call fprank(q,ff,ncof,iband,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
        do 385 i=1,ncof
          q(i,1) = q(i,1)/dmax
 385    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
c  find the coefficients in the standard b-spline representation of
c  the spherical spline.
 390    call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) if(fp) 970,970,980
        fpms = fp-s
        if(abs(fpms).le.acc) if(fp) 970,970,980
c  if f(p=inf) < s, accept the choice of knots.
        if(fpms.lt.0.) go to 580
c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 935
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
 440    do 450 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 450    continue
        do 490 num=1,nreg
          num1 = num-1
          lt = num1/npp
          l1 = lt+1
          lp = num1-lt*npp
          l2 = lp+1+ntt
          jrot = lt*np4+lp
          in = index(num)
 460      if(in.eq.0) go to 490
          store = 0.
          i1 = jrot
          do 480 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 470 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 470        continue
            i1 = i1+np4
 480      continue
          store = (w(in)*(r(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*teta(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*phi(in)
          in = nummer(in)
          go to 460
 490    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
        l1 = 1
        l2 = nrint
        if(ntest.lt.nt+1) l1=ntt+1
        if(npest.lt.np+2) l2=ntt
c  test whether we cannot further increase the number of knots.
        if(l1.gt.l2) go to 950
 500    fpmax = 0.
        l = 0
        do 510 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 510
          l = i
          fpmax = fpint(i)
 510    continue
        if(l.eq.0) go to 930
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.ntt) go to 530
c  addition in the teta-direction
        l4 = l+4
        fpint(l) = 0.
        fac1 = tt(l4)-arg
        fac2 = arg-tt(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        j = nt
        do 520 i=l4,nt
          tt(j+1) = tt(j)
          j = j-1
 520    continue
        tt(l4) = arg
        nt = nt+1
        go to 570
c  addition in the phi-direction
 530    l4 = l+4-ntt
        if(arg.lt.pi) go to 540
        arg = arg-pi
        l4 = l4-nrr
 540    fpint(l) = 0.
        fac1 = tp(l4)-arg
        fac2 = arg-tp(l4-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 500
        ll = nrr+4
        j = ll
        do 550 i=l4,ll
          tp(j+1) = tp(j)
          j = j-1
 550    continue
        tp(l4) = arg
        np = np+2
        nrr = nrr+1
        do 560 i=5,ll
          j = i+nrr
          tp(j) = tp(i)+pi
 560    continue
c  restart the computations with the new set of knots.
 570  continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spherical spline.             c
c ********************************************************             c
c we have determined the number of knots and their position. we now    c
c compute the coefficients of the smoothing spline sp(teta,phi).       c
c the observation matrix a is extended by the rows of a matrix, expres-c
c sing that sp(teta,phi) must be a constant function in the variable   c
c phi and a cubic polynomial in the variable teta. the corresponding   c
c weights of these additional rows are set to 1/(p). iteratively       c
c we than have to determine the value of p such that f(p) = sum((w(i)* c
c (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              c
c we already know that the least-squares polynomial corresponds to p=0,c
c and that the least-squares spherical spline corresponds to p=infin.  c
c the iteration process makes use of rational interpolation. since f(p)c
c is a convex and strictly decreasing function of p, it can be approx- c
c imated by a rational function of the form r(p) = (u*p+v)/(p+w).      c
c three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  c
c f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   c
c of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tt(l),l=5,...,nt-4.
 580  call fpdisc(tt,nt,5,bt,ntest)
c  evaluate the discontinuity jumps of the 3-th order derivative of
c  the b-splines at the knots tp(l),l=5,...,np-4.
      call fpdisc(tp,np,5,bp,npest)
c  initial value for p.
      p1 = 0.
      f1 = sup-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 585 i=1,ncof
        p = p+a(i,1)
 585  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband4 = iband+3
      if(ntt.le.4) iband4 = ncof
      iband3 = iband4 -1
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 920 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 600 i=1,ncof
          ff(i) = f(i)
          do 590 j=1,iband4
            q(i,j) = 0.
 590      continue
          do 600 j=1,iband
            q(i,j) = a(i,j)
 600    continue
c  extend the observation matrix with the rows of a matrix, expressing
c  that for teta=cst. sp(teta,phi) must be a constant function.
        nt6 = nt-6
        do 720 i=5,np4
          ii = i-4
          do 610 l=1,npp
             row(l) = 0.
 610      continue
          ll = ii
          do 620  l=1,5
             if(ll.gt.npp) ll=1
             row(ll) = row(ll)+bp(ii,l)
             ll = ll+1
 620      continue
          facc = 0.
          facs = 0.
          do 630 l=1,npp
             facc = facc+row(l)*coco(l)
             facs = facs+row(l)*cosi(l)
 630      continue
          do 720 j=1,nt6
c  initialize the new row.
            do 640 l=1,iband
              h(l) = 0.
 640        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            jrot = 4+(j-2)*npp
            if(j.gt.1 .and. j.lt.nt6) go to 650
            h(1) = facc
            h(2) = facs
            if(j.eq.1) jrot = 2
            go to 670
 650        do 660 l=1,npp
               h(l)=row(l)
 660        continue
 670        do 675 l=1,iband
               h(l) = h(l)*pinv
 675        continue
            ri = 0.
c  rotate the new row into triangle by givens transformations.
            do 710 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) if(i2) 720,720,690
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 720
c  apply that givens transformation to the left hand side.
              do 680 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 680          continue
 690          do 700 l=1,i2
                h(l) = h(l+1)
 700          continue
              h(i2+1) = 0.
 710        continue
 720    continue
c  extend the observation matrix with the rows of a matrix expressing
c  that for phi=cst. sp(teta,phi) must be a cubic polynomial.
        do 810 i=5,nt4
          ii = i-4
          do 810 j=1,npp
c  initialize the new row
            do 730 l=1,iband4
              h(l) = 0.
 730        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 760 l=1,5
               il = ii+l
               ij = npp
               if(il.ne.3 .and. il.ne.nt4) go to 750
               j1 = j1+3-j
               j2 = j1-2
               ij = 0
               if(il.ne.3) go to 740
               j1 = 1
               j2 = 2
               ij = j+2
 740           h(j2) = bt(ii,l)*coco(j)
               h(j2+1) = bt(ii,l)*cosi(j)
 750           h(j1) = h(j1)+bt(ii,l)
               j1 = j1+ij
 760        continue
            do 765 l=1,iband4
               h(l) = h(l)*pinv
 765        continue
            ri = 0.
            jrot = 1
            if(ii.gt.2) jrot = 3+j+(ii-3)*npp
c  rotate the new row into triangle by givens transformations.
            do 800 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) if(i2) 810,810,780
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),co,si)
c  apply that givens transformation to the right hand side.
              call fprota(co,si,ri,ff(irot))
              if(i2.eq.0) go to 810
c  apply that givens transformation to the left hand side.
              do 770 l=1,i2
                l1 = l+1
                call fprota(co,si,h(l1),q(irot,l1))
 770          continue
 780          do 790 l=1,i2
                h(l) = h(l+1)
 790          continue
              h(i2+1) = 0.
 800        continue
 810    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
        dmax = 0.
        do 820 i=1,ncof
          if(q(i,1).le.dmax) go to 820
          dmax = q(i,1)
 820    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 830 i=1,ncof
          if(q(i,1).le.sigma) go to 840
 830    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,ncc)
        rank = ncof
        go to 845
c  in case of rank deficiency, find the minimum norm solution.
 840    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 925
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,ncc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 845    do 850 i=1,ncof
           q(i,1) = q(i,1)/dmax
 850    continue
c  find the coefficients in the standard b-spline representation of
c  the spherical spline.
        call fprpsp(nt,np,coco,cosi,c,ff,ncoff)
c  compute f(p).
        fp = 0.
        do 890 num = 1,nreg
          num1 = num-1
          lt = num1/npp
          lp = num1-lt*npp
          jrot = lt*np4+lp
          in = index(num)
 860      if(in.eq.0) go to 890
          store = 0.
          i1 = jrot
          do 880 i=1,4
            hti = spt(in,i)
            j1 = i1
            do 870 j=1,4
              j1 = j1+1
              store = store+hti*spp(in,j)*c(j1)
 870        continue
            i1 = i1+np4
 880      continue
          fp = fp+(w(in)*(r(in)-store))**2
          in = nummer(in)
          go to 860
 890    continue
c  test whether the approximation sp(teta,phi) is an acceptable solution
        fpms = fp-s
        if(abs(fpms).le.acc) go to 980
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 940
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 900
        if((f2-f3).gt.acc) go to 895
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 920
 895    if(f2.lt.0.) ich3 = 1
 900    if(ich1.ne.0) go to 910
        if((f1-f2).gt.acc) go to 905
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 920
        if(p.ge.p3) p = p2*con1 +p3*con9
        go to 920
 905    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 910    if(f2.ge.f1 .or. f2.le.f3) go to 945
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 920  continue
c  error codes and messages.
 925  ier = lwest
      go to 990
 930  ier = 5
      go to 990
 935  ier = 4
      go to 990
 940  ier = 3
      go to 990
 945  ier = 2
      go to 990
 950  ier = 1
      go to 990
 960  ier = -2
      go to 990
 970  ier = -1
      fp = 0.
 980  if(ncof.ne.rank) ier = -rank
 990  return
      end
      subroutine fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wu,wv,lu,lv)
c  ..scalar arguments..
      integer idim,nu,nv,mu,mv
c  ..array arguments..
      integer lu(mu),lv(mv)
      real tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv),
     * f(mu*mv*idim),wu(mu,4),wv(mv,4)
c  ..local scalars..
      integer i,i1,j,j1,k,l,l1,l2,l3,m,nuv,nu4,nv4
      real arg,sp,tb,te
c  ..local arrays..
      real h(4)
c  ..subroutine references..
c    fpbspl
c  ..
      nu4 = nu-4
      tb = tu(4)
      te = tu(nu4+1)
      l = 4
      l1 = l+1
      do 40 i=1,mu
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tu(l1) .or. l.eq.nu4) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tu,nu,3,arg,l,h)
        lu(i) = l-4
        do 30 j=1,4
          wu(i,j) = h(j)
  30    continue
  40  continue
      nv4 = nv-4
      tb = tv(4)
      te = tv(nv4+1)
      l = 4
      l1 = l+1
      do 80 i=1,mv
        arg = v(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.tv(l1) .or. l.eq.nv4) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(tv,nv,3,arg,l,h)
        lv(i) = l-4
        do 70 j=1,4
          wv(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      nuv = nu4*nv4
      do 140 k=1,idim
        l3 = (k-1)*nuv
        do 130 i=1,mu
          l = lu(i)*nv4+l3
          do 90 i1=1,4
            h(i1) = wu(i,i1)
  90      continue
          do 120 j=1,mv
            l1 = l+lv(j)
            sp = 0.
            do 110 i1=1,4
              l2 = l1
              do 100 j1=1,4
                l2 = l2+1
                sp = sp+c(l2)*h(i1)*wv(j,j1)
 100          continue
              l1 = l1+nv4
 110        continue
            m = m+1
            f(m) = sp
 120      continue
 130    continue
 140  continue
      return
      end
      subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest,
     * nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest,
     * nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h,
     * index,nummer,wrk,lwrk,ier)
c  ..
c  ..scalar arguments..
      real xb,xe,yb,ye,s,eta,tol,fp,fp0
      integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3,
     * nc,intest,nrest,nx0,ny0,lwrk,ier
c  ..array arguments..
      real x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest),
     * coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2),
     * by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
      integer index(nrest),nummer(m)
c  ..local scalars..
      real acc,arg,cos,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv,
     * piv,p1,p2,p3,sigma,sin,sq,store,wi,x0,x1,y0,y1,zi,eps,
     * rn,one,con1,con9,con4,half,ten
      integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii,
     * in,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l,
     * la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy,nreg,
     * nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
c  ..local arrays..
      real hx(6),hy(6)
c  ..function references..
      real abs,fprati,sqrt
      integer min0
c  ..subroutine references..
c    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
c  ..
c  set constants
      one = 0.1e+01
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
      half = 0.5e0
      ten = 0.1e+02
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 1: determination of the number of knots and their position.     c
c ****************************************************************     c
c given a set of knots we compute the least-squares spline sinf(x,y),  c
c and the corresponding weighted sum of squared residuals fp=f(p=inf). c
c if iopt=-1  sinf(x,y) is the requested approximation.                c
c if iopt=0 or iopt=1 we check whether we can accept the knots:        c
c   if fp <=s we will continue with the current set of knots.          c
c   if fp > s we will increase the number of knots and compute the     c
c      corresponding least-squares spline until finally  fp<=s.        c
c the initial choice of knots depends on the value of s and iopt.      c
c   if iopt=0 we first compute the least-squares polynomial of degree  c
c     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
c     fp0=f(0) denotes the corresponding weighted sum of squared       c
c     residuals                                                        c
c   if iopt=1 we start with the knots found at the last call of the    c
c     routine, except for the case that s>=fp0; then we can compute    c
c     the least-squares polynomial directly.                           c
c eventually the independent variables x and y (and the corresponding  c
c parameters) will be switched if this can reduce the bandwidth of the c
c system to be solved.                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ichang denotes whether(1) or not(-1) the directions have been inter-
c  changed.
      ichang = -1
      x0 = xb
      x1 = xe
      y0 = yb
      y1 = ye
      kx = kxx
      ky = kyy
      kx1 = kx+1
      ky1 = ky+1
      nxe = nxest
      nye = nyest
      eps = sqrt(eta)
      if(iopt.lt.0) go to 20
c  calculation of acc, the absolute tolerance for the root of f(p)=s.
      acc = tol*s
      if(iopt.eq.0) go to 10
      if(fp0.gt.s) go to 20
c  initialization for the least-squares polynomial.
  10  nminx = 2*kx1
      nminy = 2*ky1
      nx = nminx
      ny = nminy
      ier = -2
      go to 30
  20  nx = nx0
      ny = ny0
c  main loop for the different sets of knots. m is a save upper bound
c  for the number of trials.
  30  do 420 iter=1,m
c  find the position of the additional knots which are needed for the
c  b-spline representation of s(x,y).
        l = nx
        do 40 i=1,kx1
          tx(i) = x0
          tx(l) = x1
          l = l-1
  40    continue
        l = ny
        do 50 i=1,ky1
          ty(i) = y0
          ty(l) = y1
          l = l-1
  50    continue
c  find nrint, the total number of knot intervals and nreg, the number
c  of panels in which the approximation domain is subdivided by the
c  intersection of knots.
        nxx = nx-2*kx1+1
        nyy = ny-2*ky1+1
        nrint = nxx+nyy
        nreg = nxx*nyy
c  find the bandwidth of the observation matrix a.
c  if necessary, interchange the variables x and y, in order to obtain
c  a minimal bandwidth.
        iband1 = kx*(ny-ky1)+ky
        l = ky*(nx-kx1)+kx
        if(iband1.le.l) go to 130
        iband1 = l
        ichang = -ichang
        do 60 i=1,m
          store = x(i)
          x(i) = y(i)
          y(i) = store
  60    continue
        store = x0
        x0 = y0
        y0 = store
        store = x1
        x1 = y1
        y1 = store
        n = min0(nx,ny)
        do 70 i=1,n
          store = tx(i)
          tx(i) = ty(i)
          ty(i) = store
  70    continue
        n1 = n+1
        if(nx-ny) 80,120,100
  80    do 90 i=n1,ny
          tx(i) = ty(i)
  90    continue
        go to 120
 100    do 110 i=n1,nx
          ty(i) = tx(i)
 110    continue
 120    l = nx
        nx = ny
        ny = l
        l = nxe
        nxe = nye
        nye = l
        l = nxx
        nxx = nyy
        nyy = l
        l = kx
        kx = ky
        ky = l
        kx1 = kx+1
        ky1 = ky+1
 130    iband = iband1+1
c  arrange the data points according to the panel they belong to.
        call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,index,nreg)
c  find ncof, the number of b-spline coefficients.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
c  initialize the observation matrix a.
        do 140 i=1,ncof
          f(i) = 0.
          do 140 j=1,iband
            a(i,j) = 0.
 140    continue
c  initialize the sum of squared residuals.
        fp = 0.
c  fetch the data points in the new order. main loop for the
c  different panels.
        do 250 num=1,nreg
c  fix certain constants for the current panel; jrot records the column
c  number of the first non-zero element in a row of the observation
c  matrix according to a data point of the panel.
          num1 = num-1
          lx = num1/nyy
          l1 = lx+kx1
          ly = num1-lx*nyy
          l2 = ly+ky1
          jrot = lx*nk1y+ly
c  test whether there are still data points in the panel.
          in = index(num)
 150      if(in.eq.0) go to 250
c  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
c  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
c  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
c  store the value of these b-splines in spx and spy respectively.
          do 160 i=1,kx1
            spx(in,i) = hx(i)
 160      continue
          do 170 i=1,ky1
            spy(in,i) = hy(i)
 170      continue
c  initialize the new row of observation matrix.
          do 180 i=1,iband
            h(i) = 0.
 180      continue
c  calculate the non-zero elements of the new row by making the cross
c  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do 200 i=1,kx1
            hxi = hx(i)
            j1 = i1
            do 190 j=1,ky1
              j1 = j1+1
              h(j1) = hxi*hy(j)*wi
 190        continue
            i1 = i1+nk1y
 200      continue
c  rotate the row into triangle by givens transformations .
          irot = jrot
          do 220 i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 220
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),cos,sin)
c  apply that transformation to the right hand side.
            call fprota(cos,sin,zi,f(irot))
            if(i.eq.iband) go to 230
c  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do 210 j=i3,iband
              i2 = i2+1
              call fprota(cos,sin,h(j),a(irot,i2))
 210        continue
 220      continue
c  add the contribution of the row to the sum of squares of residual
c  right hand sides.
 230      fp = fp+zi**2
c  find the number of the next data point in the panel.
 240      in = nummer(in)
          go to 150
 250    continue
c  find dmax, the maximum value for the diagonal elements in the reduced
c  triangle.
        dmax = 0.
        do 260 i=1,ncof
          if(a(i,1).le.dmax) go to 260
          dmax = a(i,1)
 260    continue
c  check whether the observation matrix is rank deficient.
        sigma = eps*dmax
        do 270 i=1,ncof
          if(a(i,1).le.sigma) go to 280
 270    continue
c  backward substitution in case of full rank.
        call fpback(a,f,ncof,iband,c,nc)
        rank = ncof
        do 275 i=1,ncof
          q(i,1) = a(i,1)/dmax
 275    continue
        go to 300
c  in case of rank deficiency, find the minimum norm solution.
c  check whether there is sufficient working space
 280    lwest = ncof*iband+ncof+iband
        if(lwrk.lt.lwest) go to 780
        do 290 i=1,ncof
          ff(i) = f(i)
          do 290 j=1,iband
            q(i,j) = a(i,j)
 290    continue
        lf =1
        lh = lf+ncof
        la = lh+iband
        call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la),
     *    wrk(lf),wrk(lh))
        do 295 i=1,ncof
          q(i,1) = q(i,1)/dmax
 295    continue
c  add to the sum of squared residuals, the contribution of reducing
c  the rank.
        fp = fp+sq
 300    if(ier.eq.(-2)) fp0 = fp
c  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 820
        fpms = fp-s
        if(abs(fpms).le.acc) if(fp) 815,815,820
c  test whether we can accept the choice of knots.
        if(fpms.lt.0.) go to 430
c  test whether we cannot further increase the number of knots.
        if(ncof.gt.m) go to 790
        ier = 0
c  search where to add a new knot.
c  find for each interval the sum of squared residuals fpint for the
c  data points having the coordinate belonging to that knot interval.
c  calculate also coord which is the same sum, weighted by the position
c  of the data points considered.
 310    do 320 i=1,nrint
          fpint(i) = 0.
          coord(i) = 0.
 320    continue
        do 360 num=1,nreg
          num1 = num-1
          lx = num1/nyy
          l1 = lx+1
          ly = num1-lx*nyy
          l2 = ly+1+nxx
          jrot = lx*nk1y+ly
          in = index(num)
 330      if(in.eq.0) go to 360
          store = 0.
          i1 = jrot
          do 350 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 340 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 340        continue
            i1 = i1+nk1y
 350      continue
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*x(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*y(in)
          in = nummer(in)
          go to 330
 360    continue
c  find the interval for which fpint is maximal on the condition that
c  there still can be added a knot.
 370    l = 0
        fpmax = 0.
        l1 = 1
        l2 = nrint
        if(nx.eq.nxe) l1 = nxx+1
        if(ny.eq.nye) l2 = nxx
        if(l1.gt.l2) go to 810
        do 380 i=l1,l2
          if(fpmax.ge.fpint(i)) go to 380
          l = i
          fpmax = fpint(i)
 380    continue
c  test whether we cannot further increase the number of knots.
        if(l.eq.0) go to 785
c  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
c  test in what direction the new knot is going to be added.
        if(l.gt.nxx) go to 400
c  addition in the x-direction.
        jxy = l+kx1
        fpint(l) = 0.
        fac1 = tx(jxy)-arg
        fac2 = arg-tx(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = nx
        do 390 i=jxy,nx
          tx(j+1) = tx(j)
          j = j-1
 390    continue
        tx(jxy) = arg
        nx = nx+1
        go to 420
c  addition in the y-direction.
 400    jxy = l+ky1-nxx
        fpint(l) = 0.
        fac1 = ty(jxy)-arg
        fac2 = arg-ty(jxy-1)
        if(fac1.gt.(ten*fac2) .or. fac2.gt.(ten*fac1)) go to 370
        j = ny
        do 410 i=jxy,ny
          ty(j+1) = ty(j)
          j = j-1
 410    continue
        ty(jxy) = arg
        ny = ny+1
c  restart the computations with the new set of knots.
 420  continue
c  test whether the least-squares polynomial is a solution of our
c  approximation problem.
 430  if(ier.eq.(-2)) go to 830
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c part 2: determination of the smoothing spline sp(x,y)                c
c *****************************************************                c
c we have determined the number of knots and their position. we now    c
c compute the b-spline coefficients of the smoothing spline sp(x,y).   c
c the observation matrix a is extended by the rows of a matrix,        c
c expressing that sp(x,y) must be a polynomial of degree kx in x and   c
c ky in y. the corresponding weights of these additional rows are set  c
c to 1./p.  iteratively we than have to determine the value of p       c
c such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
c we already know that the least-squares polynomial corresponds to     c
c p=0  and that the least-squares spline corresponds to p=infinity.    c
c the iteration process which is proposed here makes use of rational   c
c interpolation. since f(p) is a convex and strictly decreasing        c
c function of p, it can be approximated by a rational function r(p)=   c
c (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
c of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
c new value of p such that r(p)=s. convergence is guaranteed by taking c
c f1 > 0 and f3 < 0.                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kx2 = kx1+1
c  test whether there are interior knots in the x-direction.
      if(nk1x.eq.kx1) go to 440
c  evaluate the discotinuity jumps of the kx-th order derivative of
c  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
      call fpdisc(tx,nx,kx2,bx,nmax)
 440  ky2 = ky1 + 1
c  test whether there are interior knots in the y-direction.
      if(nk1y.eq.ky1) go to 450
c  evaluate the discontinuity jumps of the ky-th order derivative of
c  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
      call fpdisc(ty,ny,ky2,by,nmax)
c  initial value for p.
 450  p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = 0.
      do 460 i=1,ncof
        p = p+a(i,1)
 460  continue
      rn = ncof
      p = rn/p
c  find the bandwidth of the extended observation matrix.
      iband3 = kx1*nk1y
      iband4 = iband3 +1
      ich1 = 0
      ich3 = 0
c  iteration process to find the root of f(p)=s.
      do 770 iter=1,maxit
        pinv = one/p
c  store the triangularized observation matrix into q.
        do 480 i=1,ncof
          ff(i) = f(i)
          do 470 j=1,iband
            q(i,j) = a(i,j)
 470      continue
          ibb = iband+1
          do 480 j=ibb,iband4
            q(i,j) = 0.
 480    continue
        if(nk1y.eq.ky1) go to 560
c  extend the observation matrix with the rows of a matrix, expressing
c  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        do 550 i=ky2,nk1y
          ii = i-ky1
          do 550 j=1,nk1x
c  initialize the new row.
            do 490 l=1,iband
              h(l) = 0.
 490        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            do 500 l=1,ky2
              h(l) = by(ii,l)*pinv
 500        continue
            zi = 0.
            jrot = (j-1)*nk1y+ii
c  rotate the new row into triangle by givens transformations without
c  square roots.
            do 540 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv.eq.0.) if(i2) 550,550,520
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 550
c  apply that givens transformation to the left hand side.
              do 510 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 510          continue
 520          do 530 l=1,i2
                h(l) = h(l+1)
 530          continue
              h(i2+1) = 0.
 540        continue
 550    continue
 560    if(nk1x.eq.kx1) go to 640
c  extend the observation matrix with the rows of a matrix expressing
c  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        do 630 i=kx2,nk1x
          ii = i-kx1
          do 630 j=1,nk1y
c  initialize the new row
            do 570 l=1,iband4
              h(l) = 0.
 570        continue
c  fill in the non-zero elements of the row. jrot records the column
c  number of the first non-zero element in the row.
            j1 = 1
            do 580 l=1,kx2
              h(j1) = bx(ii,l)*pinv
              j1 = j1+nk1y
 580        continue
            zi = 0.
            jrot = (i-kx2)*nk1y+j
c  rotate the new row into triangle by givens transformations .
            do 620 irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv.eq.0.) if(i2) 630,630,600
c  calculate the parameters of the givens transformation.
              call fpgivs(piv,q(irot,1),cos,sin)
c  apply that givens transformation to the right hand side.
              call fprota(cos,sin,zi,ff(irot))
              if(i2.eq.0) go to 630
c  apply that givens transformation to the left hand side.
              do 590 l=1,i2
                l1 = l+1
                call fprota(cos,sin,h(l1),q(irot,l1))
 590          continue
 600          do 610 l=1,i2
                h(l) = h(l+1)
 610          continue
              h(i2+1) = 0.
 620        continue
 630    continue
c  find dmax, the maximum value for the diagonal elements in the
c  reduced triangle.
 640    dmax = 0.
        do 650 i=1,ncof
          if(q(i,1).le.dmax) go to 650
          dmax = q(i,1)
 650    continue
c  check whether the matrix is rank deficient.
        sigma = eps*dmax
        do 660 i=1,ncof
          if(q(i,1).le.sigma) go to 670
 660    continue
c  backward substitution in case of full rank.
        call fpback(q,ff,ncof,iband4,c,nc)
        rank = ncof
        go to 675
c  in case of rank deficiency, find the minimum norm solution.
 670    lwest = ncof*iband4+ncof+iband4
        if(lwrk.lt.lwest) go to 780
        lf = 1
        lh = lf+ncof
        la = lh+iband4
        call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la),
     *   wrk(lf),wrk(lh))
 675    do 680 i=1,ncof
          q(i,1) = q(i,1)/dmax
 680    continue
c  compute f(p).
        fp = 0.
        do 720 num = 1,nreg
          num1 = num-1
          lx = num1/nyy
          ly = num1-lx*nyy
          jrot = lx*nk1y+ly
          in = index(num)
 690      if(in.eq.0) go to 720
          store = 0.
          i1 = jrot
          do 710 i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do 700 j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
 700        continue
            i1 = i1+nk1y
 710      continue
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
          go to 690
 720    continue
c  test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).le.acc) go to 820
c  test whether the maximum allowable number of iterations has been
c  reached.
        if(iter.eq.maxit) go to 795
c  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 740
        if((f2-f3).gt.acc) go to 730
c  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 770
 730    if(f2.lt.0.) ich3 = 1
 740    if(ich1.ne.0) go to 760
        if((f1-f2).gt.acc) go to 750
c  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 770
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 770
 750    if(f2.gt.0.) ich1 = 1
c  test whether the iteration process proceeds as theoretically
c  expected.
 760    if(f2.ge.f1 .or. f2.le.f3) go to 800
c  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 770  continue
c  error codes and messages.
 780  ier = lwest
      go to 830
 785  ier = 5
      go to 830
 790  ier = 4
      go to 830
 795  ier = 3
      go to 830
 800  ier = 2
      go to 830
 810  ier = 1
      go to 830
 815  ier = -1
      fp = 0.
 820  if(ncof.ne.rank) ier = -rank
c  test whether x and y are in the original order.
 830  if(ichang.lt.0) go to 930
c  if not, interchange x and y once more.
      l1 = 1
      do 840 i=1,nk1x
        l2 = i
        do 840 j=1,nk1y
          f(l2) = c(l1)
          l1 = l1+1
          l2 = l2+nk1x
 840  continue
      do 850 i=1,ncof
        c(i) = f(i)
 850  continue
      do 860 i=1,m
        store = x(i)
        x(i) = y(i)
        y(i) = store
 860  continue
      n = min0(nx,ny)
      do 870 i=1,n
        store = tx(i)
        tx(i) = ty(i)
        ty(i) = store
 870  continue
      n1 = n+1
      if(nx-ny) 880,920,900
 880  do 890 i=n1,ny
        tx(i) = ty(i)
 890  continue
      go to 920
 900  do 910 i=n1,nx
        ty(i) = tx(i)
 910  continue
 920  l = nx
      nx = ny
      ny = l
 930  if(iopt.lt.0) go to 940
      nx0 = nx
      ny0 = ny
 940  return
      end

      subroutine fpsysy(a,n,g)
c subroutine fpsysy solves a linear n x n symmetric system
c    (a) * (b) = (g)
c on input, vector g contains the right hand side ; on output it will
c contain the solution (b).
c  ..
c  ..scalar arguments..
      integer n
c  ..array arguments..
      real a(6,6),g(6)
c  ..local scalars..
      real fac
      integer i,i1,j,k
c  ..
      g(1) = g(1)/a(1,1)
      if(n.eq.1) return
c  decomposition of the symmetric matrix (a) = (l) * (d) *(l)'
c  with (l) a unit lower triangular matrix and (d) a diagonal
c  matrix
      do 10 k=2,n
         a(k,1) = a(k,1)/a(1,1)
  10  continue
      do 40 i=2,n
         i1 = i-1
         do 30 k=i,n
            fac = a(k,i)
            do 20 j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
  20        continue
            a(k,i) = fac
            if(k.gt.i) a(k,i) = fac/a(i,i)
  30     continue
  40  continue
c  solve the system (l)*(d)*(l)'*(b) = (g).
c  first step : solve (l)*(d)*(c) = (g).
      do 60 i=2,n
         i1 = i-1
         fac = g(i)
         do 50 j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
  50     continue
         g(i) = fac/a(i,i)
  60  continue
c  second step : solve (l)'*(b) = (c)
      i = n
      do 80 j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do 70 k=i1,n
            fac = fac-g(k)*a(k,i)
  70     continue
         g(i) = fac
  80  continue
      return
      end
      subroutine fptrnp(m,mm,idim,n,nr,sp,p,b,z,a,q,right)
c  subroutine fptrnp reduces the (m+n-7) x (n-4) matrix a to upper
c  triangular form and applies the same givens transformations to
c  the (m) x (mm) x (idim) matrix z to obtain the (n-4) x (mm) x
c  (idim) matrix q
c  ..
c  ..scalar arguments..
      real p
      integer m,mm,idim,n
c  ..array arguments..
      real sp(m,4),b(n,5),z(m*mm*idim),a(n,5),q((n-4)*mm*idim),
     * right(mm*idim)
      integer nr(m)
c  ..local scalars..
      real cos,pinv,piv,sin,one
      integer i,iband,irot,it,ii,i2,i3,j,jj,l,mid,nmd,m2,m3,
     * nrold,n4,number,n1
c  ..local arrays..
      real h(7)
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      one = 1
      if(p.gt.0.) pinv = one/p
      n4 = n-4
      mid = mm*idim
      m2 = m*mm
      m3 = n4*mm
c  reduce the matrix (a) to upper triangular form (r) using givens
c  rotations. apply the same transformations to the rows of matrix z
c  to obtain the mm x (n-4) matrix g.
c  store matrix (r) into (a) and g into q.
c  initialization.
      nmd = n4*mid
      do 50 i=1,nmd
        q(i) = 0.
  50  continue
      do 100 i=1,n4
        do 100 j=1,5
          a(i,j) = 0.
 100  continue
      nrold = 0
c  iband denotes the bandwidth of the matrices (a) and (r).
      iband = 4
      do 750 it=1,m
        number = nr(it)
 150    if(nrold.eq.number) go to 300
        if(p.le.0.) go to 700
        iband = 5
c  fetch a new row of matrix (b).
        n1 = nrold+1
        do 200 j=1,5
          h(j) = b(n1,j)*pinv
 200    continue
c  find the appropriate column of q.
        do 250 j=1,mid
          right(j) = 0.
 250    continue
        irot = nrold
        go to 450
c  fetch a new row of matrix (sp).
 300    h(iband) = 0.
        do 350 j=1,4
          h(j) = sp(it,j)
 350    continue
c  find the appropriate column of q.
        j = 0
        do 400 ii=1,idim
          l = (ii-1)*m2+(it-1)*mm
          do 400 jj=1,mm
            j = j+1
            l = l+1
            right(j) = z(l)
 400    continue
        irot = number
c  rotate the new row of matrix (a) into triangle.
 450    do 600 i=1,iband
          irot = irot+1
          piv = h(i)
          if(piv.eq.0.) go to 600
c  calculate the parameters of the givens transformation.
          call fpgivs(piv,a(irot,1),cos,sin)
c  apply that transformation to the rows of matrix q.
          j = 0
          do 500 ii=1,idim
            l = (ii-1)*m3+irot
            do 500 jj=1,mm
              j = j+1
              call fprota(cos,sin,right(j),q(l))
              l = l+n4
 500      continue
c  apply that transformation to the columns of (a).
          if(i.eq.iband) go to 650
          i2 = 1
          i3 = i+1
          do 550 j=i3,iband
            i2 = i2+1
            call fprota(cos,sin,h(j),a(irot,i2))
 550      continue
 600    continue
 650    if(nrold.eq.number) go to 750
 700    nrold = nrold+1
        go to 150
 750  continue
      return
      end
      subroutine fptrpe(m,mm,idim,n,nr,sp,p,b,z,a,aa,q,right)
c  subroutine fptrpe reduces the (m+n-7) x (n-7) cyclic bandmatrix a
c  to upper triangular form and applies the same givens transformations
c  to the (m) x (mm) x (idim) matrix z to obtain the (n-7) x (mm) x
c  (idim) matrix q.
c  ..
c  ..scalar arguments..
      real p
      integer m,mm,idim,n
c  ..array arguments..
      real sp(m,4),b(n,5),z(m*mm*idim),a(n,5),aa(n,4),q((n-7)*mm*idim),
     * right(mm*idim)
      integer nr(m)
c  ..local scalars..
      real co,pinv,piv,si,one
      integer i,iband,irot,it,ii,i2,i3,j,jj,l,mid,nmd,m2,m3,
     * nrold,n4,number,n1,n7,n11,m1
c  ..local arrays..
      real h(5),h1(5),h2(4)
c  ..subroutine references..
c    fpgivs,fprota
c  ..
      one = 1
      if(p.gt.0.) pinv = one/p
      n4 = n-4
      n7 = n-7
      n11 = n-11
      mid = mm*idim
      m2 = m*mm
      m3 = n7*mm
      m1 = m-1
c  we determine the matrix (a) and then we reduce her to
c  upper triangular form (r) using givens rotations.
c  we apply the same transformations to the rows of matrix
c  z to obtain the (mm) x (n-7) matrix g.
c  we store matrix (r) into a and aa, g into q.
c  the n7 x n7 upper triangular matrix (r) has the form
c             | a1 '     |
c       (r) = |    ' a2  |
c             |  0 '     |
c  with (a2) a n7 x 4 matrix and (a1) a n11 x n11 upper
c  triangular matrix of bandwidth 5.
c  initialization.
      nmd = n7*mid
      do 50 i=1,nmd
        q(i) = 0.
  50  continue
      do 100 i=1,n4
        a(i,5) = 0.
        do 100 j=1,4
          a(i,j) = 0.
          aa(i,j) = 0.
 100  continue
      jper = 0
      nrold = 0
      do 760 it=1,m1
        number = nr(it)
 120    if(nrold.eq.number) go to 180
        if(p.le.0.) go to 740
c  fetch a new row of matrix (b).
        n1 = nrold+1
        do 140 j=1,5
          h(j) = b(n1,j)*pinv
 140    continue
c  find the appropiate row of q.
        do 160 j=1,mid
          right(j) = 0.
 160    continue
        go to 240
c  fetch a new row of matrix (sp)
 180    h(5) = 0.
        do 200 j=1,4
          h(j) = sp(it,j)
 200    continue
c  find the appropiate row of q.
        j = 0
        do 220 ii=1,idim
          l = (ii-1)*m2+(it-1)*mm
          do 220 jj=1,mm
            j = j+1
            l = l+1
            right(j) = z(l)
 220    continue
c  test whether there are non-zero values in the new row of (a)
c  corresponding to the b-splines n(j,*),j=n7+1,...,n4.
 240     if(nrold.lt.n11) go to 640
         if(jper.ne.0) go to 320
c  initialize the matrix (aa).
         jk = n11+1
         do 300 i=1,4
            ik = jk
            do 260 j=1,5
               if(ik.le.0) go to 280
               aa(ik,i) = a(ik,j)
               ik = ik-1
 260        continue
 280        jk = jk+1
 300     continue
         jper = 1
c  if one of the non-zero elements of the new row corresponds to one of
c  the b-splines n(j;*),j=n7+1,...,n4,we take account of the periodicity
c  conditions for setting up this row of (a).
 320     do 340 i=1,4
            h1(i) = 0.
            h2(i) = 0.
 340     continue
         h1(5) = 0.
         j = nrold-n11
         do 420 i=1,5
            j = j+1
            l0 = j
 360        l1 = l0-4
            if(l1.le.0) go to 400
            if(l1.le.n11) go to 380
            l0 = l1-n11
            go to 360
 380        h1(l1) = h(i)
            go to 420
 400        h2(l0) = h2(l0) + h(i)
 420     continue
c  rotate the new row of (a) into triangle.
         if(n11.le.0) go to 560
c  rotations with the rows 1,2,...,n11 of (a).
         do 540 irot=1,n11
            piv = h1(1)
            i2 = min0(n11-irot,4)
            if(piv.eq.0.) go to 500
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the columns of matrix q.
            j = 0
            do 440 ii=1,idim
               l = (ii-1)*m3+irot
               do 440 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 440        continue
c  apply that transformation to the rows of (a) with respect to aa.
            do 460 i=1,4
               call fprota(co,si,h2(i),aa(irot,i))
 460        continue
c  apply that transformation to the rows of (a) with respect to a.
            if(i2.eq.0) go to 560
            do 480 i=1,i2
               i1 = i+1
               call fprota(co,si,h1(i1),a(irot,i1))
 480        continue
 500        do 520 i=1,i2
               h1(i) = h1(i+1)
 520        continue
            h1(i2+1) = 0.
 540     continue
c  rotations with the rows n11+1,...,n7 of a.
 560     do 620 irot=1,4
            ij = n11+irot
            if(ij.le.0) go to 620
            piv = h2(irot)
            if(piv.eq.0.) go to 620
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,aa(ij,irot),co,si)
c  apply that transformation to the columns of matrix q.
            j = 0
            do 580 ii=1,idim
               l = (ii-1)*m3+ij
               do 580 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 580        continue
            if(irot.eq.4) go to 620
c  apply that transformation to the rows of (a) with respect to aa.
            j1 = irot+1
            do 600 i=j1,4
               call fprota(co,si,h2(i),aa(ij,i))
 600        continue
 620     continue
         go to 720
c  rotation into triangle of the new row of (a), in case the elements
c  corresponding to the b-splines n(j;*),j=n7+1,...,n4 are all zero.
 640     irot =nrold
         do 700 i=1,5
            irot = irot+1
            piv = h(i)
            if(piv.eq.0.) go to 700
c  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),co,si)
c  apply that transformation to the columns of matrix g.
            j = 0
            do 660 ii=1,idim
               l = (ii-1)*m3+irot
               do 660 jj=1,mm
                 j = j+1
                 call fprota(co,si,right(j),q(l))
                 l = l+n7
 660        continue
c  apply that transformation to the rows of (a).
            if(i.eq.5) go to 700
            i2 = 1
            i3 = i+1
            do 680 j=i3,5
               i2 = i2+1
               call fprota(co,si,h(j),a(irot,i2))
 680        continue
 700     continue
 720     if(nrold.eq.number) go to 760
 740     nrold = nrold+1
         go to 120
 760  continue
      return
      end
      subroutine insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
c  subroutine insert inserts a new knot x into a spline function s(x)
c  of degree k and calculates the b-spline representation of s(x) with
c  respect to the new set of knots. in addition, if iopt.ne.0, s(x)
c  will be considered as a periodic spline with period per=t(n-k)-t(k+1)
c  satisfying the boundary constraints
c       t(i+n-2*k-1) = t(i)+per  ,i=1,2,...,2*k+1
c       c(i+n-2*k-1) = c(i)      ,i=1,2,...,k
c  in that case, the knots and b-spline coefficients returned will also
c  satisfy these boundary constraints, i.e.
c       tt(i+nn-2*k-1) = tt(i)+per  ,i=1,2,...,2*k+1
c       cc(i+nn-2*k-1) = cc(i)      ,i=1,2,...,k
c
c  calling sequence:
c     call insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
c
c  input parameters:
c    iopt : integer flag, specifying whether (iopt.ne.0) or not (iopt=0)
c           the given spline must be considered as being periodic.
c    t    : array,length nest, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length nest, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : real, which gives the location of the knot to be inserted.
c    nest : integer specifying the dimension of the arrays t,c,tt and cc
c           nest > n.
c
c  output parameters:
c    tt   : array,length nest, which contains the position of the knots
c           after insertion.
c    nn   : integer, giving the total number of knots after insertion
c    cc   : array,length nest, which contains the b-spline coefficients
c           of s(x) with respect to the new set of knots.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    nest > n
c    t(k+1) <= x <= t(n-k)
c    in case of a periodic spline (iopt.ne.0) there must be
c       either at least k interior knots t(j) satisfying t(k+1)<t(j)<=x
c       or at least k interior knots t(j) satisfying x<=t(j)<t(n-k)
c
c  other subroutines required: fpinst.
c
c  further comments:
c   subroutine insert may be called as follows
c        call insert(iopt,t,n,c,k,x,t,n,c,nest,ier)
c   in which case the new representation will simply replace the old one
c
c  references :
c    boehm w : inserting new knots into b-spline curves. computer aided
c              design 12 (1980) 199-201.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer iopt,n,k,nn,nest,ier
      real x
c  ..array arguments..
      real t(nest),c(nest),tt(nest),cc(nest)
c  ..local scalars..
      integer kk,k1,l,nk,nk1
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nest.le.n) go to 40
      k1 = k+1
      nk = n-k
      if(x.lt.t(k1) .or. x.gt.t(nk)) go to 40
c  search for knot interval t(l) <= x < t(l+1).
      nk1 = nk-1
      l = k1
  10  if(x.lt.t(l+1) .or. l.eq.nk1) go to 20
      l = l+1
      go to 10
  20  if(t(l).ge.t(l+1)) go to 40
      if(iopt.eq.0) go to 30
      kk = 2*k
      if(l.le.kk .and. l.ge.(n-kk)) go to 40
  30  ier = 0
c  insert the new knot.
      call fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)
  40  return
      end
      subroutine parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,
     * nc,c,fp,wrk,lwrk,iwrk,ier)
c  given the ordered set of m points x(i) in the idim-dimensional space
c  and given also a corresponding set of strictly increasing values u(i)
c  and the set of positive numbers w(i),i=1,2,...,m, subroutine parcur
c  determines a smooth approximating spline curve s(u), i.e.
c      x1 = s1(u)
c      x2 = s2(u)       ub <= u <= ue
c      .........
c      xidim = sidim(u)
c  with sj(u),j=1,2,...,idim spline functions of degree k with common
c  knots t(j),j=1,2,...,n.
c  if ipar=1 the values ub,ue and u(i),i=1,2,...,m must be supplied by
c  the user. if ipar=0 these values are chosen automatically by parcur
c  as  v(1) = 0
c      v(i) = v(i-1) + dist(x(i),x(i-1)) ,i=2,3,...,m
c      u(i) = v(i)/v(m) ,i=1,2,...,m
c      ub = u(1) = 0, ue = u(m) = 1.
c  if iopt=-1 parcur calculates the weighted least-squares spline curve
c  according to a given set of knots.
c  if iopt>=0 the number of knots of the splines sj(u) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(u) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(u) is given in the b-spline representation and can be
c  evaluated by means of subroutine curev.
c
c  calling sequence:
c     call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,
c    * fp,wrk,lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline curve (iopt=-1) or a smoothing spline
c           curve (iopt=0 or 1) must be determined.if iopt=0 the routine
c           will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
c           i=1,2,...,k+1. if iopt=1 the routine will continue with the
c           knots found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   ipar  : integer flag. on entry ipar must specify whether (ipar=1)
c           the user will supply the parameter values u(i),ub and ue
c           or whether (ipar=0) these values are to be calculated by
c           parcur. unchanged on exit.
c   idim  : integer. on entry idim must specify the dimension of the
c           curve. 0 < idim < 11.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > k. unchanged on exit.
c   u     : real array of dimension at least (m). in case ipar=1,before
c           entry, u(i) must be set to the i-th value of the parameter
c           variable u for i=1,2,...,m. these values must then be
c           supplied in strictly ascending order and will be unchanged
c           on exit. in case ipar=0, on exit,array u will contain the
c           values u(i) as determined by parcur.
c   mx    : integer. on entry mx must specify the actual dimension of
c           the array x as declared in the calling (sub)program. mx must
c           not be too small (see x). unchanged on exit.
c   x     : real array of dimension at least idim*m.
c           before entry, x(idim*(i-1)+j) must contain the j-th coord-
c           inate of the i-th data point for i=1,2,...,m and j=1,2,...,
c           idim. unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. unchanged on exit.
c           see also further comments.
c   ub,ue : real values. on entry (in case ipar=1) ub and ue must
c           contain the lower and upper bound for the parameter u.
c           ub <=u(1), ue>= u(m). if ipar = 0 these values will
c           automatically be set to 0 and 1 by parcur.
c   k     : integer. on entry k must specify the degree of the splines.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the splines returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+k+1, the number of knots
c           needed for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the smoothing spline curve returned
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline curve,i.e. the position of the interior knots t(k+2),
c           t(k+3),..,t(n-k-1) as well as the position of the additional
c           t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   nc    : integer. on entry nc must specify the actual dimension of
c           the array c as declared in the calling (sub)program. nc
c           must not be too small (see c). unchanged on exit.
c   c     : real array of dimension at least (nest*idim).
c           on succesful exit, this array will contain the coefficients
c           in the b-spline representation of the spline curve s(u),i.e.
c           the b-spline coefficients of the spline sj(u) will be given
c           in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
c   fp    : real. unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline curve returned.
c   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the curve returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the curve returned is an interpolating
c             spline curve (fp=0).
c    ier=-2 : normal return. the curve returned is the weighted least-
c             squares polynomial curve of degree k.in this extreme case
c             fp gives the upper bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares spline
c             curve according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline curve
c             with fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing curve
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
c             0<=ipar<=1, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
c             nc>=nest*idim
c             if ipar=0: sum j=1,idim (x(idim*i+j)-x(idim*(i-1)+j))**2>0
c                        i=1,2,...,m-1.
c             if ipar=1: ub<=u(1)<u(2)<...<u(m)<=ue
c             if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
c                         ub<t(k+2)<t(k+3)<...<t(n-k-1)<ue
c                            (ub=0 and ue=1 in case ipar=0)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points uu(j) such that
c                         t(j) < uu(j) < t(j+k+1), j=1,2,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+k+1
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the curve will be too smooth and signal will be
c   lost ; if s is too small the curve will pick up too much noise. in
c   the extreme cases the program will return an interpolating curve if
c   s=0 and the least-squares polynomial curve of degree k if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in x(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial curve and the upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximating curve shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if parcur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   curve underlying the data. but, if the computation mode iopt=1 is
c   used, the knots returned may also depend on the s-values at previous
c   calls (if these were smaller). therefore, if after a number of
c   trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   parcur once more with the selected value for s but now with iopt=0.
c   indeed, parcur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c   the form of the approximating curve can strongly be affected by
c   the choice of the parameter values u(i). if there is no physical
c   reason for choosing a particular parameter u, often good results
c   will be obtained with the choice of parcur (in case ipar=0), i.e.
c        v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
c   where
c        q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
c   other possibilities for q(i) are
c        q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
c        q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= max j=1,idim abs(xj(i)-xj(i-1))
c        q(i)= 1
c
c  other subroutines required:
c    fpback,fpbspl,fpchec,fppara,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real ub,ue,s,fp
      integer iopt,ipar,idim,m,mx,k,nest,n,nc,lwrk,ier
c  ..array arguments..
      real u(m),x(mx),w(m),t(nest),c(nc),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real tol,dist
      integer i,ia,ib,ifp,ig,iq,iz,i1,i2,j,k1,k2,lwest,maxit,nmin,ncc
c ..function references
      real sqrt
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(ipar.lt.0 .or. ipar.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      k2 = k1+1
      nmin = 2*k1
      if(m.lt.k1 .or. nest.lt.nmin) go to 90
      ncc = nest*idim
      if(mx.lt.m*idim .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if(lwrk.lt.lwest) go to 90
      if(ipar.ne.0 .or. iopt.gt.0) go to 40
      i1 = 0
      i2 = idim
      u(1) = 0.
      do 20 i=2,m
         dist = 0.
         do 10 j=1,idim
            i1 = i1+1
            i2 = i2+1
            dist = dist+(x(i2)-x(i1))**2
  10     continue
         u(i) = u(i-1)+sqrt(dist)
  20  continue
      if(u(m).le.0.) go to 90
      do 30 i=2,m
         u(i) = u(i)/u(m)
  30  continue
      ub = 0.
      ue = 1.
      u(m) = ue
  40  if(ub.gt.u(1) .or. ue.lt.u(m) .or. w(1).le.0.) go to 90
      do 50 i=2,m
         if(u(i-1).ge.u(i) .or. w(i).le.0.) go to 90
  50  continue
      if(iopt.ge.0) go to 70
      if(n.lt.nmin .or. n.gt.nest) go to 90
      j = n
      do 60 i=1,k1
         t(i) = ub
         t(j) = ue
         j = j-1
  60  continue
      call fpchec(u,m,t,n,k,ier)
      if(ier) 90,80,90
  70  if(s.lt.0.) go to 90
      if(s.eq.0. .and. nest.lt.(m+k1)) go to 90
      ier = 0
c we partition the working space and determine the spline curve.
  80  ifp = 1
      iz = ifp+nest
      ia = iz+ncc
      ib = ia+nest*k1
      ig = ib+nest*k2
      iq = ig+nest*k2
      call fppara(iopt,idim,m,u,mx,x,w,ub,ue,k,s,nest,tol,maxit,k1,k2,
     * n,t,ncc,c,fp,wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),
     * iwrk,ier)
  90  return
      end
      subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,
     * wrk,lwrk,iwrk,kwrk,ier)
c  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my the partial derivative ( order nux,nuy) of a bivariate spline
c  s(x,y) of degrees kx and ky, given in the b-spline representation.
c
c  calling sequence:
c     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   nux   : integer values, specifying the order of the partial
c   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
c   x     : real array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : real array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of the
c           specified partial derivative of s(x,y) at the point
c           (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
c   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,m,m0,m1,
     * nc,nkx1,nky1,nxx,nyy
      real ak,fac
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux.lt.0 .or. nux.ge.kx) go to 400
      if(nuy.lt.0 .or. nuy.ge.ky) go to 400
      lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
      if(lwrk.lt.lwest) go to 400
      if(kwrk.lt.(mx+my)) go to 400
      if(mx-1) 400,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 400
  20  continue
  30  if(my-1) 400,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 400
  50  continue
  60  ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
c  the partial derivative of order (nux,nuy) of a bivariate spline of
c  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
c  we calculate the b-spline coefficients of this spline
      do 70 i=1,nc
        wrk(i) = c(i)
  70  continue
      if(nux.eq.0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac.le.0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy.eq.0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac.le.0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
c  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+mx*(kx1-nux)
      call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky,
     * x,mx,y,my,z,wrk(iwx),wrk(iwy),iwrk(1),iwrk(mx+1))
 400  return
      end

      subroutine parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,
     * nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c  given the set of ordered points f(i,j) in the idim-dimensional space,
c  corresponding to grid values (u(i),v(j)) ,i=1,...,mu ; j=1,...,mv,
c  parsur determines a smooth approximating spline surface s(u,v) , i.e.
c    f1 = s1(u,v)
c      ...                u(1) <= u <= u(mu) ; v(1) <= v <= v(mv)
c    fidim = sidim(u,v)
c  with sl(u,v), l=1,2,...,idim bicubic spline functions with common
c  knots tu(i),i=1,...,nu in the u-variable and tv(j),j=1,...,nv in the
c  v-variable.
c  in addition, these splines will be periodic in the variable u if
c  ipar(1) = 1 and periodic in the variable v if ipar(2) = 1.
c  if iopt=-1, parsur determines the least-squares bicubic spline
c  surface according to a given set of knots.
c  if iopt>=0, the number of knots of s(u,v) and their position
c  is chosen automatically by the routine. the smoothness of s(u,v) is
c  achieved by minimalizing the discontinuity jumps of the derivatives
c  of the splines at the knots. the amount of smoothness of s(u,v) is
c  determined by the condition that
c  fp=sumi=1,mu(sumj=1,mv(dist(f(i,j)-s(u(i),v(j)))**2))<=s,
c  with s a given non-negative constant.
c  the fit s(u,v) is given in its b-spline representation and can be
c  evaluated by means of routine surev.
c
c calling sequence:
c     call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu,
c    *  nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. unchanged on exit.
c          on entry iopt must specify whether a least-squares surface
c          (iopt=-1) or a smoothing surface (iopt=0 or 1)must be
c          determined.
c          if iopt=0 the routine will start with the initial set of
c          knots needed for determining the least-squares polynomial
c          surface.
c          if iopt=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately
c          preceded by another call with iopt = 1 or iopt = 0.
c  ipar  : integer array of dimension 2. unchanged on exit.
c          on entry ipar(1) must specify whether (ipar(1)=1) or not
c          (ipar(1)=0) the splines must be periodic in the variable u.
c          on entry ipar(2) must specify whether (ipar(2)=1) or not
c          (ipar(2)=0) the splines must be periodic in the variable v.
c  idim  : integer. on entry idim must specify the dimension of the
c          surface. 1 <= idim <= 3. unchanged on exit.
c  mu    : integer. on entry mu must specify the number of grid points
c          along the u-axis. unchanged on exit.
c          mu >= mumin where mumin=4-2*ipar(1)
c  u     : real array of dimension at least (mu). before entry, u(i)
c          must be set to the u-co-ordinate of the i-th grid point
c          along the u-axis, for i=1,2,...,mu. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  mv    : integer. on entry mv must specify the number of grid points
c          along the v-axis. unchanged on exit.
c          mv >= mvmin where mvmin=4-2*ipar(2)
c  v     : real array of dimension at least (mv). before entry, v(j)
c          must be set to the v-co-ordinate of the j-th grid point
c          along the v-axis, for j=1,2,...,mv. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  f     : real array of dimension at least (mu*mv*idim).
c          before entry, f(mu*mv*(l-1)+mv*(i-1)+j) must be set to the
c          l-th co-ordinate of the data point corresponding to the
c          the grid point (u(i),v(j)) for l=1,...,idim ,i=1,...,mu
c          and j=1,...,mv. unchanged on exit.
c          if ipar(1)=1 it is expected that f(mu*mv*(l-1)+mv*(mu-1)+j)
c          = f(mu*mv*(l-1)+j), l=1,...,idim ; j=1,...,mv
c          if ipar(2)=1 it is expected that f(mu*mv*(l-1)+mv*(i-1)+mv)
c          = f(mu*mv*(l-1)+mv*(i-1)+1), l=1,...,idim ; i=1,...,mu
c  s     : real. on entry (if iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = mu/2, nvest=mv/2, will
c          be sufficient. always large enough are nuest=mu+4+2*ipar(1),
c          nvest = mv+4+2*ipar(2), the number of knots needed for
c          interpolation (s=0). see also further comments.
c  nu    : integer.
c          unless ier=10 (in case iopt>=0), nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          surface returned. if the computation mode iopt=1 is used,
c          the value of nu should be left unchanged between subsequent
c          calls. in case iopt=-1, the value of nu should be specified
c          on entry.
c  tu    : real array of dimension at least (nuest).
c          on succesful exit, this array will contain the knots of the
c          splines with respect to the u-variable, i.e. the position of
c          the interior knots tu(5),...,tu(nu-4) as well as the position
c          of the additional knots tu(1),...,tu(4) and tu(nu-3),...,
c          tu(nu) needed for the b-spline representation.
c          if the computation mode iopt=1 is used,the values of tu(1)
c          ...,tu(nu) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tu(5),
c          ...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          surface returned. if the computation mode iopt=1 is used,
c          the value of nv should be left unchanged between subsequent
c          calls. in case iopt=-1, the value of nv should be specified
c          on entry.
c  tv    : real array of dimension at least (nvest).
c          on succesful exit, this array will contain the knots of the
c          splines with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt=1 is used,the values of tv(1)
c          ...,tv(nv) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tv(5),
c          ...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nuest-4)*(nvest-4)*idim.
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v)
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline surface returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of
c          wrk(1),...,wrk(4) should be left unchanged between subsequent
c          calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+
c           4*(mu+mv)+q*idim where q is the larger of mv and nuest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of
c          iwrk(1),.,iwrk(3) should be left unchanged between subsequent
c          calls.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 3+mu+mv+nuest+nvest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the surface returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline surface returned is an
c            interpolating surface (fp=0).
c   ier=-2 : normal return. the surface returned is the least-squares
c            polynomial surface. in this extreme case fp gives the
c            upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares surface
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing surface with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing surface
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 0<=ipar(1)<=1, 0<=ipar(2)<=1, 1 <=idim<=3
c            mu >= 4-2*ipar(1),mv >= 4-2*ipar(2), nuest >=8, nvest >= 8,
c            kwrk>=3+mu+mv+nuest+nvest,
c            lwrk >= 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))
c             +4*(mu+mv)+max(nuest,mv)*idim
c            u(i-1)<u(i),i=2,..,mu, v(i-1)<v(i),i=2,...,mv
c            if iopt=-1: 8<=nu<=min(nuest,mu+4+2*ipar(1))
c                        u(1)<tu(5)<tu(6)<...<tu(nu-4)<u(mu)
c                        8<=nv<=min(nvest,mv+4+2*ipar(2))
c                        v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(mv)
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates uu(p) and vv(q) such
c                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
c                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
c                     (see fpchec or fpchep)
c            if iopt>=0: s>=0
c                       if s=0: nuest>=mu+4+2*ipar(1)
c                               nvest>=mv+4+2*ipar(2)
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the surface will be too smooth and signal will be
c   lost ; if s is too small the surface will pick up too much noise. in
c   the extreme cases the program will return an interpolating surface
c   if s=0 and the constrained least-squares polynomial surface if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, parsur will return a surface s(u,v) which
c   exactly reproduces the surface underlying the data he can evaluate
c   the sum(dist(f(i,j)-s(u(i),v(j)))**2) to find a good estimate for s.
c   for example, if he knows that the statistical errors on his f(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mu*mv*(0.1)**2.
c   if nothing is known about the statistical error in f(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the le-sq polynomial surface and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt = 1 the program will continue with the knots found at
c   the last call of the routine. this will save a lot of computation
c   time if parsur is called repeatedly for different values of s.
c   the number of knots of the surface returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   surface underlying the data. if the computation mode iopt = 1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   parsur once more with the chosen value for s but now with iopt=0.
c   indeed, parsur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in parsur the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nuest=8 (the lowest allowable value for
c   nuest), the user can indicate that he wants an approximation with
c   splines which are simple cubic polynomials in the variable u.
c
c  other subroutines required:
c    fppasu,fpchec,fpchep,fpknot,fprati,fpgrpa,fptrnp,fpback,
c    fpbacp,fpbspl,fptrpe,fpdisc,fpgivs,fprota
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real s,fp
      integer iopt,idim,mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
c  ..array arguments..
      real u(mu),v(mv),f(mu*mv*idim),tu(nuest),tv(nvest),
     * c((nuest-4)*(nvest-4)*idim),wrk(lwrk)
      integer ipar(2),iwrk(kwrk)
c  ..local scalars..
      real tol,ub,ue,vb,ve,peru,perv
      integer i,j,jwrk,kndu,kndv,knru,knrv,kwest,l1,l2,l3,l4,
     * lfpu,lfpv,lwest,lww,maxit,nc,mf,mumin,mvmin
c  ..function references..
      integer max0
c  ..subroutine references..
c    fppasu,fpchec,fpchep
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 200
      if(ipar(1).lt.0 .or. ipar(1).gt.1) go to 200
      if(ipar(2).lt.0 .or. ipar(2).gt.1) go to 200
      if(idim.le.0 .or. idim.gt.3) go to 200
      mumin = 4-2*ipar(1)
      if(mu.lt.mumin .or. nuest.lt.8) go to 200
      mvmin = 4-2*ipar(2)
      if(mv.lt.mvmin .or. nvest.lt.8) go to 200
      mf = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 4+nuest*(mv*idim+11+4*ipar(1))+nvest*(11+4*ipar(2))+
     * 4*(mu+mv)+max0(nuest,mv)*idim
      kwest = 3+mu+mv+nuest+nvest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 200
      do 10 i=2,mu
        if(u(i-1).ge.u(i)) go to 200
  10  continue
      do 20 i=2,mv
        if(v(i-1).ge.v(i)) go to 200
  20  continue
      if(iopt.ge.0) go to 100
      if(nu.lt.8 .or. nu.gt.nuest) go to 200
      ub = u(1)
      ue = u(mu)
      if (ipar(1).ne.0) go to 40
      j = nu
      do 30 i=1,4
        tu(i) = ub
        tu(j) = ue
        j = j-1
  30  continue
      call fpchec(u,mu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
      go to 60
  40  l1 = 4
      l2 = l1
      l3 = nu-3
      l4 = l3
      peru = ue-ub
      tu(l2) = ub
      tu(l3) = ue
      do 50 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tu(l2) = tu(l4)-peru
        tu(l3) = tu(l1)+peru
  50  continue
      call fpchep(u,mu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
  60  if(nv.lt.8 .or. nv.gt.nvest) go to 200
      vb = v(1)
      ve = v(mv)
      if (ipar(2).ne.0) go to 80
      j = nv
      do 70 i=1,4
        tv(i) = vb
        tv(j) = ve
        j = j-1
  70  continue
      call fpchec(v,mv,tv,nv,3,ier)
      if(ier.ne.0) go to 200
      go to 150
  80  l1 = 4
      l2 = l1
      l3 = nv-3
      l4 = l3
      perv = ve-vb
      tv(l2) = vb
      tv(l3) = ve
      do 90 j=1,3
        l1 = l1+1
        l2 = l2-1
        l3 = l3+1
        l4 = l4-1
        tv(l2) = tv(l4)-perv
        tv(l3) = tv(l1)+perv
  90  continue
      call fpchep(v,mv,tv,nv,3,ier)
      if(ier) 200,150,200
 100  if(s.lt.0.) go to 200
      if(s.eq.0. .and. (nuest.lt.(mu+4+2*ipar(1)) .or.
     * nvest.lt.(mv+4+2*ipar(2))) )go to 200
      ier = 0
c  we partition the working space and determine the spline approximation
 150  lfpu = 5
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-4-nuest-nvest
      knru = 4
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fppasu(iopt,ipar,idim,u,mu,v,mv,f,mf,s,nuest,nvest,
     * tol,maxit,nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
     * wrk(lfpu),wrk(lfpv),iwrk(1),iwrk(2),iwrk(3),iwrk(knru),
     * iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),jwrk,ier)
 200  return
      end

      subroutine percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,
     * wrk,lwrk,iwrk,ier)
c  given the set of data points (x(i),y(i)) and the set of positive
c  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth
c  periodic spline approximation of degree k with period per=x(m)-x(1).
c  if iopt=-1 percur calculates the weighted least-squares periodic
c  spline according to a given set of knots.
c  if iopt>=0 the number of knots of the spline s(x) and the position
c  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
c  ness of s(x) is then achieved by minimalizing the discontinuity
c  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
c  n-k-1. the amount of smoothness is determined by the condition that
c  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
c  negative constant, called the smoothing factor.
c  the fit s(x) is given in the b-spline representation (b-spline coef-
c  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
c  subroutine splev.
c
c  calling sequence:
c     call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,
c    * lwrk,iwrk,ier)
c
c  parameters:
c   iopt  : integer flag. on entry iopt must specify whether a weighted
c           least-squares spline (iopt=-1) or a smoothing spline (iopt=
c           0 or 1) must be determined. if iopt=0 the routine will start
c           with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1),
c           i=1,2,...,2*k+2. if iopt=1 the routine will continue with
c           the knots found at the last call of the routine.
c           attention: a call with iopt=1 must always be immediately
c           preceded by another call with iopt=1 or iopt=0.
c           unchanged on exit.
c   m     : integer. on entry m must specify the number of data points.
c           m > 1. unchanged on exit.
c   x     : real array of dimension at least (m). before entry, x(i)
c           must be set to the i-th value of the independent variable x,
c           for i=1,2,...,m. these values must be supplied in strictly
c           ascending order. x(m) only indicates the length of the
c           period of the spline, i.e per=x(m)-x(1).
c           unchanged on exit.
c   y     : real array of dimension at least (m). before entry, y(i)
c           must be set to the i-th value of the dependent variable y,
c           for i=1,2,...,m-1. the element y(m) is not used.
c           unchanged on exit.
c   w     : real array of dimension at least (m). before entry, w(i)
c           must be set to the i-th value in the set of weights. the
c           w(i) must be strictly positive. w(m) is not used.
c           see also further comments. unchanged on exit.
c   k     : integer. on entry k must specify the degree of the spline.
c           1<=k<=5. it is recommended to use cubic splines (k=3).
c           the user is strongly dissuaded from choosing k even,together
c           with a small s-value. unchanged on exit.
c   s     : real.on entry (in case iopt>=0) s must specify the smoothing
c           factor. s >=0. unchanged on exit.
c           for advice on the choice of s see further comments.
c   nest  : integer. on entry nest must contain an over-estimate of the
c           total number of knots of the spline returned, to indicate
c           the storage space available to the routine. nest >=2*k+2.
c           in most practical situation nest=m/2 will be sufficient.
c           always large enough is nest=m+2*k,the number of knots needed
c           for interpolation (s=0). unchanged on exit.
c   n     : integer.
c           unless ier = 10 (in case iopt >=0), n will contain the
c           total number of knots of the spline approximation returned.
c           if the computation mode iopt=1 is used this value of n
c           should be left unchanged between subsequent calls.
c           in case iopt=-1, the value of n must be specified on entry.
c   t     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the knots of the
c           spline,i.e. the position of the interior knots t(k+2),t(k+3)
c           ...,t(n-k-1) as well as the position of the additional knots
c           t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for
c           the b-spline representation.
c           if the computation mode iopt=1 is used, the values of t(1),
c           t(2),...,t(n) should be left unchanged between subsequent
c           calls. if the computation mode iopt=-1 is used, the values
c           t(k+2),...,t(n-k-1) must be supplied by the user, before
c           entry. see also the restrictions (ier=10).
c   c     : real array of dimension at least (nest).
c           on succesful exit, this array will contain the coefficients
c           c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
c   fp    : real. unless ier = 10, fp contains the weighted sum of
c           squared residuals of the spline approximation returned.
c   wrk   : real array of dimension at least (m*(k+1)+nest*(8+5*k)).
c           used as working space. if the computation mode iopt=1 is
c           used, the values wrk(1),...,wrk(n) should be left unchanged
c           between subsequent calls.
c   lwrk  : integer. on entry,lwrk must specify the actual dimension of
c           the array wrk as declared in the calling (sub)program. lwrk
c           must not be too small (see wrk). unchanged on exit.
c   iwrk  : integer array of dimension at least (nest).
c           used as working space. if the computation mode iopt=1 is
c           used,the values iwrk(1),...,iwrk(n) should be left unchanged
c           between subsequent calls.
c   ier   : integer. unless the routine detects an error, ier contains a
c           non-positive value on exit, i.e.
c    ier=0  : normal return. the spline returned has a residual sum of
c             squares fp such that abs(fp-s)/s <= tol with tol a relat-
c             ive tolerance set to 0.001 by the program.
c    ier=-1 : normal return. the spline returned is an interpolating
c             periodic spline (fp=0).
c    ier=-2 : normal return. the spline returned is the weighted least-
c             squares constant. in this extreme case fp gives the upper
c             bound fp0 for the smoothing factor s.
c    ier=1  : error. the required storage space exceeds the available
c             storage space, as specified by the parameter nest.
c             probably causes : nest too small. if nest is already
c             large (say nest > m/2), it may also indicate that s is
c             too small
c             the approximation returned is the least-squares periodic
c             spline according to the knots t(1),t(2),...,t(n). (n=nest)
c             the parameter fp gives the corresponding weighted sum of
c             squared residuals (fp>s).
c    ier=2  : error. a theoretically impossible result was found during
c             the iteration proces for finding a smoothing spline with
c             fp = s. probably causes : s too small.
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=3  : error. the maximal number of iterations maxit (set to 20
c             by the program) allowed for finding a smoothing spline
c             with fp=s has been reached. probably causes : s too small
c             there is an approximation returned but the corresponding
c             weighted sum of squared residuals does not satisfy the
c             condition abs(fp-s)/s < tol.
c    ier=10 : error. on entry, the input data are controlled on validity
c             the following restrictions must be satisfied.
c             -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1
c             x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k)
c             if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
c                         x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m)
c                       the schoenberg-whitney conditions, i.e. there
c                       must be a subset of data points xx(j) with
c                       xx(j) = x(i) or x(i)+(x(m)-x(1)) such that
c                         t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1
c             if iopt>=0: s>=0
c                         if s=0 : nest >= m+2*k
c             if one of these conditions is found to be violated,control
c             is immediately repassed to the calling program. in that
c             case there is no approximation returned.
c
c  further comments:
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating periodic
c   spline if s=0 and the weighted least-squares constant if s is very
c   large. between these extremes, a properly chosen s will result in
c   a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in y(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   constant and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if percur is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. but, if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   percur once more with the selected value for s but now with iopt=0.
c   indeed, percur may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c
c  other subroutines required:
c    fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota
c
c  references:
c   dierckx p. : algorithms for smoothing data with periodic and
c                parametric splines, computer graphics and image
c                processing 20 (1982) 171-184.
c   dierckx p. : algorithms for smoothing data with periodic and param-
c                etric splines, report tw55, dept. computer science,
c                k.u.leuven, 1981.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real s,fp
      integer iopt,m,k,nest,n,lwrk,ier
c  ..array arguments..
      real x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      integer iwrk(nest)
c  ..local scalars..
      real per,tol
      integer i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,lwest,
     * maxit,m1,nmin
c  ..subroutine references..
c    perper,pcheck
c  ..
c  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.2 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(8+5*k)
      if(lwrk.lt.lwest) go to 50
      m1 = m-1
      do 10 i=1,m1
         if(x(i).ge.x(i+1) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.le.nmin .or. n.gt.nest) go to 50
      per = x(m)-x(1)
      j1 = k1
      t(j1) = x(1)
      i1 = n-k
      t(i1) = x(m)
      j2 = j1
      i2 = i1
      do 20 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  20  continue
      call fpchep(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+2*k)) go to 50
      ier = 0
c we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia1 = iz+nest
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,
     * wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1),wrk(ig2),
     * wrk(iq),iwrk,ier)
  50  return
      end
      subroutine pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,
     * nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c  subroutine pogrid fits a function f(x,y) to a set of data points
c  z(i,j) given at the nodes (x,y)=(u(i)*cos(v(j)),u(i)*sin(v(j))),
c  i=1,...,mu ; j=1,...,mv , of a radius-angle grid over a disc
c    x ** 2  +  y ** 2  <=  r ** 2 .
c
c  this approximation problem is reduced to the determination of a
c  bicubic spline s(u,v) smoothing the data (u(i),v(j),z(i,j)) on the
c  rectangle 0<=u<=r, v(1)<=v<=v(1)+2*pi
c  in order to have continuous partial derivatives
c              i+j
c             d   f(0,0)
c    g(i,j) = ----------
c                i   j
c              dx  dy
c
c  s(u,v)=f(x,y) must satisfy the following conditions
c
c    (1) s(0,v) = g(0,0)   v(1)<=v<= v(1)+2*pi
c
c        d s(0,v)
c    (2) -------- = cos(v)*g(1,0)+sin(v)*g(0,1)  v(1)<=v<= v(1)+2*pi
c        d u
c
c  moreover, s(u,v) must be periodic in the variable v, i.e.
c
c         j            j
c        d s(u,vb)   d s(u,ve)
c    (3) ---------- = ---------   0 <=u<= r, j=0,1,2 , vb=v(1),
c           j            j                             ve=vb+2*pi
c        d v          d v
c
c  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
c  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
c  smoothness of s(u,v) is achieved by minimalizing the discontinuity
c  jumps of the derivatives of the spline at the knots. the amount of
c  smoothness of s(u,v) is determined by the condition that
c  fp=sumi=1,mu(sumj=1,mv((z(i,j)-s(u(i),v(j)))**2))+(z0-g(0,0))**2<=s,
c  with s a given non-negative constant.
c  the fit s(u,v) is given in its b-spline representation and can be
c  evaluated by means of routine bispev. f(x,y) = s(u,v) can also be
c  evaluated by means of function program evapol.
c
c calling sequence:
c     call pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,nu,tu,
c    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer array of dimension 3, specifying different options.
c          unchanged on exit.
c  iopt(1):on entry iopt(1) must specify whether a least-squares spline
c          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be
c          determined.
c          if iopt(1)=0 the routine will start with an initial set of
c          knots tu(i)=0,tu(i+4)=r,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi,i=1,.
c          ...,8.
c          if iopt(1)=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt(1)=1 must always be immediately
c          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
c  iopt(2):on entry iopt(2) must specify the requested order of conti-
c          nuity for f(x,y) at the origin.
c          if iopt(2)=0 only condition (1) must be fulfilled and
c          if iopt(2)=1 conditions (1)+(2) must be fulfilled.
c  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
c          (iopt(3)=0) the approximation f(x,y) must vanish at the
c          boundary of the approximation domain.
c  ider  : integer array of dimension 2, specifying different options.
c          unchanged on exit.
c  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not
c          (ider(1)=-1) there is a data value z0 at the origin.
c          if ider(1)=1, z0 will be considered to be the right function
c          value, and it will be fitted exactly (g(0,0)=z0=c(1)).
c          if ider(1)=0, z0 will be considered to be a data value just
c          like the other data values z(i,j).
c  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not
c          (ider(2)=0) f(x,y) must have vanishing partial derivatives
c          g(1,0) and g(0,1) at the origin. (in case iopt(2)=1)
c  mu    : integer. on entry mu must specify the number of grid points
c          along the u-axis. unchanged on exit.
c          mu >= mumin where mumin=4-iopt(3)-ider(2) if ider(1)<0
c                                 =3-iopt(3)-ider(2) if ider(1)>=0
c  u     : real array of dimension at least (mu). before entry, u(i)
c          must be set to the u-co-ordinate of the i-th grid point
c          along the u-axis, for i=1,2,...,mu. these values must be
c          positive and supplied in strictly ascending order.
c          unchanged on exit.
c  mv    : integer. on entry mv must specify the number of grid points
c          along the v-axis. mv > 3 . unchanged on exit.
c  v     : real array of dimension at least (mv). before entry, v(j)
c          must be set to the v-co-ordinate of the j-th grid point
c          along the v-axis, for j=1,2,...,mv. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c          -pi <= v(1) < pi , v(mv) < v(1)+2*pi.
c  z     : real array of dimension at least (mu*mv).
c          before entry, z(mv*(i-1)+j) must be set to the data value at
c          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv.
c          unchanged on exit.
c  z0    : real value. on entry (if ider(1) >=0 ) z0 must specify the
c          data value at the origin. unchanged on exit.
c  r     : real value. on entry r must specify the radius of the disk.
c          r>=u(mu) (>u(mu) if iopt(3)=1). unchanged on exit.
c  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = mu/2, nvest=mv/2, will
c          be sufficient. always large enough are nuest=mu+5+iopt(2)+
c          iopt(3), nvest = mv+7, the number of knots needed for
c          interpolation (s=0). see also further comments.
c  nu    : integer.
c          unless ier=10 (in case iopt(1)>=0), nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1 is
c          used, the value of nu should be left unchanged between sub-
c          sequent calls. in case iopt(1)=-1, the value of nu should be
c          specified on entry.
c  tu    : real array of dimension at least (nuest).
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the u-variable, i.e. the position of
c          the interior knots tu(5),...,tu(nu-4) as well as the position
c          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...=
c          tu(nu)=r needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used,the values of tu(1)
c          ...,tu(nu) should be left unchanged between subsequent calls.
c          if the computation mode iopt(1)=-1 is used, the values tu(5),
c          ...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt(1)>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1 is
c          used, the value of nv should be left unchanged between sub-
c          sequent calls. in case iopt(1) = -1, the value of nv should
c          be specified on entry.
c  tv    : real array of dimension at least (nvest).
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used,the values of tv(1)
c          ...,tv(nv) should be left unchanged between subsequent calls.
c          if the computation mode iopt(1)=-1 is used, the values tv(5),
c          ...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nuest-4)*(nvest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v)
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline approximation returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt(1)=1 is used the values of
c          wrk(1),...,wrk(8) should be left unchanged between subsequent
c          calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 8+nuest*(mv+nvest+3)+nvest*21+4*mu+6*mv+q
c           where q is the larger of (mv+nvest) and nuest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt(1)=1 is used the values of
c          iwrk(1),.,iwrk(4) should be left unchanged between subsequent
c          calls.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 4+mu+mv+nuest+nvest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the least-squares
c            constrained polynomial. in this extreme case fp gives the
c            upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares spline
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1,
c            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0.
c            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8,
c            kwrk>=4+mu+mv+nuest+nvest,
c            lwrk >= 8+nuest*(mv+nvest+3)+nvest*21+4*mu+6*mv+
c             max(nuest,mv+nvest)
c            0< u(i-1)<u(i)<=r,i=2,..,mu, (< r if iopt(3)=1)
c            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv
c            if iopt(1)=-1: 8<=nu<=min(nuest,mu+5+iopt(2)+iopt(3))
c                           0<tu(5)<tu(6)<...<tu(nu-4)<r
c                           8<=nv<=min(nvest,mv+7)
c                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates uu(p) and vv(q) such
c                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
c                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value
c                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
c                     (vv(q) is either a value v(j) or v(j)+2*pi)
c            if iopt(1)>=0: s>=0
c                       if s=0: nuest>=mu+5+iopt(2)+iopt(3), nvest>=mv+7
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   pogrid does not allow individual weighting of the data-values.
c   so, if these were determined to widely different accuracies, then
c   perhaps the general data set routine polar should rather be used
c   in spite of efficiency.
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, pogrid will return a spline s(u,v) which
c   exactly reproduces the function underlying the data he can evaluate
c   the sum((z(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s
c   for example, if he knows that the statistical errors on his z(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mu*mv*(0.1)**2.
c   if nothing is known about the statistical error in z(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the least-squares polynomial and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt(1)=0.
c   if iopt(1) = 1 the program will continue with the knots found at
c   the last call of the routine. this will save a lot of computation
c   time if pogrid is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt(1) = 1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt(1)=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   pogrid once more with the chosen value for s but now with iopt(1)=0.
c   indeed, pogrid may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in pogrid the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nuest=8 (the lowest allowable value for
c   nuest), the user can indicate that he wants an approximation which
c   is a simple cubic polynomial in the variable u.
c
c  other subroutines required:
c    fppogr,fpchec,fpchep,fpknot,fpopdi,fprati,fpgrdi,fpsysy,fpback,
c    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota
c
c  references:
c   dierckx p. : fast algorithms for smoothing data over a disc or a
c                sphere using tensor product splines, in "algorithms
c                for approximation", ed. j.c.mason and m.g.cox,
c                clarendon press oxford, 1987, pp. 51-65
c   dierckx p. : fast algorithms for smoothing data over a disc or a
c                sphere using tensor product splines, report tw73, dept.
c                computer science,k.u.leuven, 1985.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : july 1985
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real z0,r,s,fp
      integer mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
c  ..array arguments..
      integer iopt(3),ider(2),iwrk(kwrk)
      real u(mu),v(mv),z(mu*mv),c((nuest-4)*(nvest-4)),tu(nuest),
     * tv(nvest),wrk(lwrk)
c  ..local scalars..
      real per,pi,tol,uu,ve,zmax,zmin,one,half,rn,zb
      integer i,i1,i2,j,jwrk,j1,j2,kndu,kndv,knru,knrv,kwest,l,
     * ldz,lfpu,lfpv,lwest,lww,m,maxit,mumin,muu,nc
c  ..function references..
      real atan2
      integer max0
c  ..subroutine references..
c    fpchec,fpchep,fppogr
c  ..
c  set constants
      one = 1
      half = 0.5e0
      pi = atan2(0.,-one)
      per = pi+pi
      ve = v(1)+per
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations, a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt(1).lt.(-1) .or. iopt(1).gt.1) go to 200
      if(iopt(2).lt.0 .or. iopt(2).gt.1) go to 200
      if(iopt(3).lt.0 .or. iopt(3).gt.1) go to 200
      if(ider(1).lt.(-1) .or. ider(1).gt.1) go to 200
      if(ider(2).lt.0 .or. ider(2).gt.1) go to 200
      if(ider(2).eq.1 .and. iopt(2).eq.0) go to 200
      mumin = 4-iopt(3)-ider(2)
      if(ider(1).ge.0) mumin = mumin-1
      if(mu.lt.mumin .or. mv.lt.4) go to 200
      if(nuest.lt.8 .or. nvest.lt.8) go to 200
      m = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 8+nuest*(mv+nvest+3)+21*nvest+4*mu+6*mv+
     * max0(nuest,mv+nvest)
      kwest = 4+mu+mv+nuest+nvest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 200
      if(u(1).le.0. .or. u(mu).gt.r) go to 200
      if(iopt(3).eq.0) go to 10
      if(u(mu).eq.r) go to 200
  10  if(mu.eq.1) go to 30
      do 20 i=2,mu
        if(u(i-1).ge.u(i)) go to 200
  20  continue
  30  if(v(1).lt. (-pi) .or. v(1).ge.pi ) go to 200
      if(v(mv).ge.v(1)+per) go to 200
      do 40 i=2,mv
        if(v(i-1).ge.v(i)) go to 200
  40  continue
      if(iopt(1).gt.0) go to 140
c  if not given, we compute an estimate for z0.
      if(ider(1).lt.0) go to 50
      zb = z0
      go to 70
  50  zb = 0.
      do 60 i=1,mv
         zb = zb+z(i)
  60  continue
      rn = mv
      zb = zb/rn
c  we determine the range of z-values.
  70  zmin = zb
      zmax = zb
      do 80 i=1,m
         if(z(i).lt.zmin) zmin = z(i)
         if(z(i).gt.zmax) zmax = z(i)
  80  continue
      wrk(5) = zb
      wrk(6) = 0.
      wrk(7) = 0.
      wrk(8) = zmax -zmin
      iwrk(4) = mu
      if(iopt(1).eq.0) go to 140
      if(nu.lt.8 .or. nu.gt.nuest) go to 200
      if(nv.lt.11 .or. nv.gt.nvest) go to 200
      j = nu
      do 90 i=1,4
        tu(i) = 0.
        tu(j) = r
        j = j-1
  90  continue
      l = 9
      wrk(l) = 0.
      if(iopt(2).eq.0) go to 100
      l = l+1
      uu = u(1)
      if(uu.gt.tu(5)) uu = tu(5)
      wrk(l) = uu*half
 100  do 110 i=1,mu
        l = l+1
        wrk(l) = u(i)
 110  continue
      if(iopt(3).eq.0) go to 120
      l = l+1
      wrk(l) = r
 120  muu = l-8
      call fpchec(wrk(9),muu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
      j1 = 4
      tv(j1) = v(1)
      i1 = nv-3
      tv(i1) = ve
      j2 = j1
      i2 = i1
      do 130 i=1,3
        i1 = i1+1
        i2 = i2-1
        j1 = j1+1
        j2 = j2-1
        tv(j2) = tv(i2)-per
        tv(i1) = tv(j1)+per
 130  continue
      l = 9
      do 135 i=1,mv
        wrk(l) = v(i)
        l = l+1
 135  continue
      wrk(l) = ve
      call fpchep(wrk(9),mv+1,tv,nv,3,ier)
      if(ier) 200,150,200
 140  if(s.lt.0.) go to 200
      if(s.eq.0. .and. (nuest.lt.(mu+5+iopt(2)+iopt(3)) .or.
     * nvest.lt.(mv+7)) ) go to 200
c  we partition the working space and determine the spline approximation
 150  ldz = 5
      lfpu = 9
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-8-nuest-nvest
      knru = 5
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fppogr(iopt,ider,u,mu,v,mv,z,m,zb,r,s,nuest,nvest,tol,maxit,
     * nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),wrk(lfpu),
     * wrk(lfpv),wrk(ldz),wrk(8),iwrk(1),iwrk(2),iwrk(3),iwrk(4),
     * iwrk(knru),iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),jwrk,ier)
 200  return
      end

      subroutine polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
     *  nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  subroutine polar fits a smooth function f(x,y) to a set of data
c  points (x(i),y(i),z(i)) scattered arbitrarily over an approximation
c  domain  x**2+y**2 <= rad(atan(y/x))**2. through the transformation
c    x = u*rad(v)*cos(v) , y = u*rad(v)*sin(v)
c  the approximation problem is reduced to the determination of a bi-
c  cubic spline s(u,v) fitting a corresponding set of data points
c  (u(i),v(i),z(i)) on the rectangle 0<=u<=1,-pi<=v<=pi.
c  in order to have continuous partial derivatives
c              i+j
c             d   f(0,0)
c    g(i,j) = ----------
c                i   j
c              dx  dy
c
c  s(u,v)=f(x,y) must satisfy the following conditions
c
c    (1) s(0,v) = g(0,0)   -pi <=v<= pi.
c
c        d s(0,v)
c    (2) -------- = rad(v)*(cos(v)*g(1,0)+sin(v)*g(0,1))
c        d u
c                                                    -pi <=v<= pi
c         2
c        d s(0,v)         2       2             2
c    (3) -------- = rad(v)*(cos(v)*g(2,0)+sin(v)*g(0,2)+sin(2*v)*g(1,1))
c           2
c        d u                                         -pi <=v<= pi
c
c  moreover, s(u,v) must be periodic in the variable v, i.e.
c
c         j            j
c        d s(u,-pi)   d s(u,pi)
c    (4) ---------- = ---------   0 <=u<= 1, j=0,1,2
c           j           j
c        d v         d v
c
c  if iopt(1) < 0 circle calculates a weighted least-squares spline
c  according to a given set of knots in u- and v- direction.
c  if iopt(1) >=0, the number of knots in each direction and their pos-
c  ition tu(j),j=1,2,...,nu ; tv(j),j=1,2,...,nv are chosen automatical-
c  ly by the routine. the smoothness of s(u,v) is then achieved by mini-
c  malizing the discontinuity jumps of the derivatives of the spline
c  at the knots. the amount of smoothness of s(u,v) is determined  by
c  the condition that fp = sum((w(i)*(z(i)-s(u(i),v(i))))**2) be <= s,
c  with s a given non-negative constant.
c  the bicubic spline is given in its standard b-spline representation
c  and the corresponding function f(x,y) can be evaluated by means of
c  function program evapol.
c
c calling sequence:
c     call polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
c    *  nv,tv,u,v,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer array of dimension 3, specifying different options.
c          unchanged on exit.
c  iopt(1):on entry iopt(1) must specify whether a weighted
c          least-squares polar spline (iopt(1)=-1) or a smoothing
c          polar spline (iopt(1)=0 or 1) must be determined.
c          if iopt(1)=0 the routine will start with an initial set of
c          knots tu(i)=0,tu(i+4)=1,i=1,...,4;tv(i)=(2*i-9)*pi,i=1,...,8.
c          if iopt(1)=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt(1)=1 must always be immediately
c          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
c  iopt(2):on entry iopt(2) must specify the requested order of conti-
c          nuity for f(x,y) at the origin.
c          if iopt(2)=0 only condition (1) must be fulfilled,
c          if iopt(2)=1 conditions (1)+(2) must be fulfilled and
c          if iopt(2)=2 conditions (1)+(2)+(3) must be fulfilled.
c  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
c          (iopt(3)=0) the approximation f(x,y) must vanish at the
c          boundary of the approximation domain.
c  m     : integer. on entry m must specify the number of data points.
c          m >= 4-iopt(2)-iopt(3) unchanged on exit.
c  x     : real array of dimension at least (m).
c  y     : real array of dimension at least (m).
c  z     : real array of dimension at least (m).
c          before entry, x(i),y(i),z(i) must be set to the co-ordinates
c          of the i-th data point, for i=1,...,m. the order of the data
c          points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  rad   : real function subprogram defining the boundary of the approx-
c          imation domain, i.e   x = rad(v)*cos(v) , y = rad(v)*sin(v),
c          -pi <= v <= pi.
c          must be declared external in the calling (sub)program.
c  s     : real. on entry (in case iopt(1) >=0) s must specify the
c          smoothing factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions resp.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = nvest = 8+sqrt(m/2) will
c          be sufficient. see also further comments.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nu    : integer.
c          unless ier=10 (in case iopt(1) >=0),nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1
c          is used, the value of nu should be left unchanged between
c          subsequent calls.
c          in case iopt(1)=-1,the value of nu must be specified on entry
c  tu    : real array of dimension at least nuest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the u-variable, i.e. the position
c          of the interior knots tu(5),...,tu(nu-4) as well as the
c          position of the additional knots tu(1)=...=tu(4)=0 and
c          tu(nu-3)=...=tu(nu)=1 needed for the b-spline representation
c          if the computation mode iopt(1)=1 is used,the values of
c          tu(1),...,tu(nu) should be left unchanged between subsequent
c          calls. if the computation mode iopt(1)=-1 is used,the values
c          tu(5),...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt(1)>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1
c          is used, the value of nv should be left unchanged between
c          subsequent calls. in case iopt(1)=-1, the value of nv should
c          be specified on entry.
c  tv    : real array of dimension at least nvest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used, the values of
c          tv(1),...,tv(nv) should be left unchanged between subsequent
c          calls. if the computation mode iopt(1)=-1 is used,the values
c          tv(5),...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  u     : real array of dimension at least (m).
c  v     : real array of dimension at least (m).
c          on succesful exit, u(i),v(i) contains the co-ordinates of
c          the i-th data point with respect to the transformed rectan-
c          gular approximation domain, for i=1,2,...,m.
c          if the computation mode iopt(1)=1 is used the values of
c          u(i),v(i) should be left unchanged between subsequent calls.
c  c     : real array of dimension at least (nuest-4)*(nvest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v).
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt(1)=1 is used the value of
c          wrk1(1) should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
c          values d(i)/max(d(i)),i=1,...,ncof=1+iopt(2)*(iopt(2)+3)/2+
c          (nv-7)*(nu-5-iopt(2)-iopt(3)) with d(i) the i-th diagonal el-
c          ement of the triangular matrix for calculating the b-spline
c          coefficients.it includes those elements whose square is < eps
c          which are treated as 0 in the case of rank deficiency(ier=-2)
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            k = nuest-7, l = nvest-7, p = 1+iopt(2)*(iopt(2)+3)/2,
c            q = k+2-iopt(2)-iopt(3) then
c          lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper bound  for lwrk2 = (p+l*q+1)*(4*l+p)
c          +p+l*q where p,l,q are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(nuest-7)*(nvest-7).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares constrained polynomial . in this extreme case
c            fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)-iopt(3))+ier
c            is large the results may be inaccurate.
c            they could also seriously depend on the value of eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small or badly chosen eps.
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=4  : error. no more knots can be added because the dimension
c            of the spline 1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)
c            -iopt(3)) already exceeds the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            polar spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt(1)<=1 , 0<=iopt(2)<=2 , 0<=iopt(3)<=1 ,
c            m>=4-iopt(2)-iopt(3) , nuest>=8 ,nvest >=8, 0<eps<1,
c            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
c            lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
c            kwrk >= m+(nuest-7)*(nvest-7)
c            if iopt(1)=-1:9<=nu<=nuest,9+iopt(2)*(iopt(2)+1)<=nv<=nvest
c                          0<tu(5)<tu(6)<...<tu(nu-4)<1
c                          -pi<tv(5)<tv(6)<...<tv(nv-4)<pi
c            if iopt(1)>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nu,
c            nv,tu,tv,wrk1,u,v and having adjusted the value of lwrk2
c            and the dimension of the array wrk2 accordingly, the user
c            can continue at the point the program was left, by calling
c            polar with iopt(1)=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the constrained weighted least-squares polynomial if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in z(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to choose s very small is strongly discouraged. this considerably
c   increases computation time and memory requirements. it may also
c   cause rank-deficiency (ier<-2) and endager numerical stability.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt(1)=0.
c   if iopt(1)=1 the program will continue with the set of knots found
c   at the last call of the routine. this will save a lot of computation
c   time if polar is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt(1)=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt(1)=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   polar once more with the selected value for s but now with iopt(1)=0
c   indeed, polar may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in polar the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c
c  other subroutines required:
c    fpback,fpbspl,fppola,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
c    fprppo
c
c  references:
c   dierckx p.: an algorithm for fitting data over a circle using tensor
c               product splines,j.comp.appl.maths 15 (1986) 161-173.
c   dierckx p.: an algorithm for fitting data on a circle using tensor
c               product splines, report tw68, dept. computer science,
c               k.u.leuven, 1984.
c   dierckx p.: curve and surface fitting with splines, monographs on
c               numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : june 1984
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real s,eps,fp
      integer m,nuest,nvest,nu,nv,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real x(m),y(m),z(m),w(m),tu(nuest),tv(nvest),u(m),v(m),
     * c((nuest-4)*(nvest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iopt(3),iwrk(kwrk)
c  ..user specified function
      real rad
c  ..local scalars..
      real tol,pi,dist,r,one
      integer i,ib1,ib3,ki,kn,kwest,la,lbu,lcc,lcs,lro,j
     * lbv,lco,lf,lff,lfp,lh,lq,lsu,lsv,lwest,maxit,ncest,ncc,nuu,
     * nvv,nreg,nrint,nu4,nv4,iopt1,iopt2,iopt3,ipar,nvmin
c  ..function references..
      real atan2,sqrt
      external rad
c  ..subroutine references..
c    fppola
c  ..
c  set up constants
      one = 1
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 60
      iopt1 = iopt(1)
      if(iopt1.lt.(-1) .or. iopt1.gt.1) go to 60
      iopt2 = iopt(2)
      if(iopt2.lt.0 .or. iopt2.gt.2) go to 60
      iopt3 = iopt(3)
      if(iopt3.lt.0 .or. iopt3.gt.1) go to 60
      if(m.lt.(4-iopt2-iopt3)) go to 60
      if(nuest.lt.8 .or. nvest.lt.8) go to 60
      nu4 = nuest-4
      nv4 = nvest-4
      ncest = nu4*nv4
      nuu = nuest-7
      nvv = nvest-7
      ipar = 1+iopt2*(iopt2+3)/2
      ncc = ipar+nvv*(nuest-5-iopt2-iopt3)
      nrint = nuu+nvv
      nreg = nuu*nvv
      ib1 = 4*nvv
      ib3 = ib1+ipar
      lwest = ncc*(1+ib1+ib3)+2*nrint+ncest+m*8+ib3+5*nuest+12*nvest
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 60
      if(iopt1.gt.0) go to 40
      do 10 i=1,m
        if(w(i).le.0.) go to 60
        dist = x(i)**2+y(i)**2
        u(i) = 0.
        v(i) = 0.
        if(dist.le.0.) go to 10
        v(i) = atan2(y(i),x(i))
        r = rad(v(i))
        if(r.le.0.) go to 60
        u(i) = sqrt(dist)/r
        if(u(i).gt.one) go to 60
  10  continue
      if(iopt1.eq.0) go to 40
      nuu = nu-8
      if(nuu.lt.1 .or. nu.gt.nuest) go to 60
      tu(4) = 0.
      do 20 i=1,nuu
         j = i+4
         if(tu(j).le.tu(j-1) .or. tu(j).ge.one) go to 60
  20  continue
      nvv = nv-8
      nvmin = 9+iopt2*(iopt2+1)
      if(nv.lt.nvmin .or. nv.gt.nvest) go to 60
      pi = atan2(0.,-one)
      tv(4) = -pi
      do 30 i=1,nvv
         j = i+4
         if(tv(j).le.tv(j-1) .or. tv(j).ge.pi) go to 60
  30  continue
      go to 50
  40  if(s.lt.0.) go to 60
  50  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbu = lh+ib3
      lbv = lbu+5*nuest
      lro = lbv+5*nvest
      lcc = lro+nvest
      lcs = lcc+nvest
      lsu = lcs+nvest*5
      lsv = lsu+m*4
      call fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s,nuest,nvest,eps,tol,
     * maxit,ib1,ib3,ncest,ncc,nrint,nreg,nu,tu,nv,tv,c,fp,wrk1(1),
     * wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),
     * wrk1(lcs),wrk1(la),wrk1(lq),wrk1(lbu),wrk1(lbv),wrk1(lsu),
     * wrk1(lsv),wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  60  return
      end

      subroutine profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
c  if iopt=0 subroutine profil calculates the b-spline coefficients of
c  the univariate spline f(y) = s(u,y) with s(x,y) a bivariate spline of
c  degrees kx and ky, given in the b-spline representation.
c  if iopt = 1 it calculates the b-spline coefficients of the univariate
c  spline g(x) = s(x,u)
c
c  calling sequence:
c     call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
c
c  input parameters:
c   iopt  : integer flag, specifying whether the profile f(y) (iopt=0)
c           or the profile g(x) (iopt=1) must be determined.
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   u     : real value, specifying the requested profile.
c           tx(kx+1)<=u<=tx(nx-kx), if iopt=0.
c           ty(ky+1)<=u<=ty(ny-ky), if iopt=1.
c   nu    : on entry nu must specify the dimension of the array cu.
c           nu >= ny if iopt=0, nu >= nx if iopt=1.
c
c  output parameters:
c   cu    : real array of dimension (nu).
c           on succesful exit this array contains the b-spline
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   if iopt=0 : tx(kx+1) <= u <= tx(nx-kx), nu >=ny.
c   if iopt=1 : ty(ky+1) <= u <= ty(ny-ky), nu >=nx.
c
c  other subroutines required:
c    fpbspl
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer iopt,nx,ny,kx,ky,nu,ier
      real u
c  ..array arguments..
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),cu(nu)
c  ..local scalars..
      integer i,j,kx1,ky1,l,l1,m,m0,nkx1,nky1
      real sum
c  ..local array
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      ier = 10
      if(iopt.ne.0) go to 200
      if(nu.lt.ny) go to 300
      if(u.lt.tx(kx1) .or. u.gt.tx(nkx1+1)) go to 300
c  the b-splinecoefficients of f(y) = s(u,y).
      ier = 0
      l = kx1
      l1 = l+1
 110  if(u.lt.tx(l1) .or. l.eq.nkx1) go to 120
      l = l1
      l1 = l+1
      go to 110
 120  call fpbspl(tx,nx,kx,u,l,h)
      m0 = (l-kx1)*nky1+1
      do 140 i=1,nky1
        m = m0
        sum = 0.
        do 130 j=1,kx1
          sum = sum+h(j)*c(m)
          m = m+nky1
 130    continue
        cu(i) = sum
        m0 = m0+1
 140  continue
      go to 300
 200  if(nu.lt.nx) go to 300
      if(u.lt.ty(ky1) .or. u.gt.ty(nky1+1)) go to 300
c  the b-splinecoefficients of g(x) = s(x,u).
      ier = 0
      l = ky1
      l1 = l+1
 210  if(u.lt.ty(l1) .or. l.eq.nky1) go to 220
      l = l1
      l1 = l+1
      go to 210
 220  call fpbspl(ty,ny,ky,u,l,h)
      m0 = l-ky
      do 240 i=1,nkx1
        m = m0
        sum = 0.
        do 230 j=1,ky1
          sum = sum+h(j)*c(m)
          m = m+1
 230    continue
        cu(i) = sum
        m0 = m0+nky1
 240  continue
 300  return
      end

      subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,
     * nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
c i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 regrid calculates the least-squares spline according
c to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
c constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a least-
c          squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
c          must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0 and
c                     s.ne.0.
c          unchanged on exit.
c  mx    : integer. on entry mx must specify the number of grid points
c          along the x-axis. mx > kx . unchanged on exit.
c  x     : real array of dimension at least (mx). before entry, x(i)
c          must be set to the x-co-ordinate of the i-th grid point
c          along the x-axis, for i=1,2,...,mx. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  my    : integer. on entry my must specify the number of grid points
c          along the y-axis. my > ky . unchanged on exit.
c  y     : real array of dimension at least (my). before entry, y(j)
c          must be set to the y-co-ordinate of the j-th grid point
c          along the y-axis, for j=1,2,...,my. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c  z     : real array of dimension at least (mx*my).
c          before entry, z(my*(i-1)+j) must be set to the data value at
c          the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
c          unchanged on exit.
c  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
c          unchanged on exit.
c  kx,ky : integer values. on entry kx and ky must specify the degrees
c          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
c          (kx=ky=3) splines. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nxest : integer. unchanged on exit.
c  nyest : integer. unchanged on exit.
c          on entry, nxest and nyest must specify an upper bound for the
c          number of knots required in the x- and y-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
c          in most practical situation nxest = mx/2, nyest=my/2, will
c          be sufficient. always large enough are nxest=mx+kx+1, nyest=
c          my+ky+1, the number of knots needed for interpolation (s=0).
c          see also further comments.
c  nx    : integer.
c          unless ier=10 (in case iopt >=0), nx will contain the total
c          number of knots with respect to the x-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of nx should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of nx should be specified on entry
c  tx    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the x-variable, i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
c          position of the additional knots tx(1)=...=tx(kx+1)=xb and
c          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of tx(1),
c          ...,tx(nx) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tx(kx+2),
c          ...tx(nx-kx-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  ny    : integer.
c          unless ier=10 (in case iopt >=0), ny will contain the total
c          number of knots with respect to the y-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of ny should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of ny should be specified on entry
c  ty    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the y-variable, i.e. the position of
c          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
c          position of the additional knots ty(1)=...=ty(ky+1)=yb and
c          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of ty(1),
c          ...,ty(ny) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values ty(ky+2),
c          ...ty(ny-ky-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(x,y)
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline approximation returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of wrk(1),
c          ...,wrk(4) should be left unchanged between subsequent calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c            my*(ky+1) +u
c           where u is the larger of my and nxest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt=1 is used the values of iwrk(1),
c          ...,iwrk(3) should be left unchanged between subsequent calls
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 3+mx+my+nxest+nyest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the least-squares
c            polynomial of degrees kx and ky. in this extreme case fp
c            gives the upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares spline
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
c            nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
c            lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c             my*(ky+1) +max(my,nxest),
c            xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
c            if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=min(nyest,my+ky+1)
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates xx(p) and yy(q) such
c                    that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
c                           ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
c            if iopt>=0: s>=0
c                        if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   regrid does not allow individual weighting of the data-values.
c   so, if these were determined to widely different accuracies, then
c   perhaps the general data set routine surfit should rather be used
c   in spite of efficiency.
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the least-squares polynomial (degrees kx,ky) if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, regrid will return a spline s(x,y) which
c   exactly reproduces the function underlying the data he can evaluate
c   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
c   for example, if he knows that the statistical errors on his z(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mx*my*(0.1)**2.
c   if nothing is known about the statistical error in z(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the least-squares polynomial and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if regrid is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   regrid once more with the selected value for s but now with iopt=0.
c   indeed, regrid may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in regrid the number of knots
c   in one direction (say nx) has reached the value of its upper bound
c   (nxest), then from that moment on all subsequent knots are added
c   in the other (y) direction. this may indicate that the value of
c   nxest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nxest=2*kx+2 (the lowest allowable value for
c   nxest), the user can indicate that he wants an approximation which
c   is a simple polynomial of degree kx in the variable x.
c
c  other subroutines required:
c    fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
c    fpknot
c
c  references:
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, siam j.numer.anal.
c                19 (1982) 1286-1304.
c   dierckx p. : a fast algorithm for smoothing data on a rectangular
c                grid while using spline functions, report tw53, dept.
c                computer science,k.u.leuven, 1980.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real xb,xe,yb,ye,s,fp
      integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
c  ..array arguments..
      real x(mx),y(my),z(mx*my),tx(nxest),ty(nyest),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk)
      integer iwrk(kwrk)
c  ..local scalars..
      real tol
      integer i,j,jwrk,kndx,kndy,knrx,knry,kwest,kx1,kx2,ky1,ky2,
     * lfpx,lfpy,lwest,lww,maxit,nc,nminx,nminy,mz
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpregr,fpchec
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(kx.le.0 .or. kx.gt.5) go to 70
      kx1 = kx+1
      kx2 = kx1+1
      if(ky.le.0 .or. ky.gt.5) go to 70
      ky1 = ky+1
      ky2 = ky1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
      nminx = 2*kx1
      if(mx.lt.kx1 .or. nxest.lt.nminx) go to 70
      nminy = 2*ky1
      if(my.lt.ky1 .or. nyest.lt.nminy) go to 70
      mz = mx*my
      nc = (nxest-kx1)*(nyest-ky1)
      lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+
     * my*ky1+max0(nxest,my)
      kwest = 3+mx+my+nxest+nyest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 70
      if(xb.gt.x(1) .or. xe.lt.x(mx)) go to 70
      do 10 i=2,mx
        if(x(i-1).ge.x(i)) go to 70
  10  continue
      if(yb.gt.y(1) .or. ye.lt.y(my)) go to 70
      do 20 i=2,my
        if(y(i-1).ge.y(i)) go to 70
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 70
      j = nx
      do 30 i=1,kx1
        tx(i) = xb
        tx(j) = xe
        j = j-1
  30  continue
      call fpchec(x,mx,tx,nx,kx,ier)
      if(ier.ne.0) go to 70
      if(ny.lt.nminy .or. ny.gt.nyest) go to 70
      j = ny
      do 40 i=1,ky1
        ty(i) = yb
        ty(j) = ye
        j = j-1
  40  continue
      call fpchec(y,my,ty,ny,ky,ier)
      if(ier) 70,60,70
  50  if(s.lt.0.) go to 70
      if(s.eq.0. .and. (nxest.lt.(mx+kx1) .or. nyest.lt.(my+ky1)) )
     * go to 70
      ier = 0
c  we partition the working space and determine the spline approximation
  60  lfpx = 5
      lfpy = lfpx+nxest
      lww = lfpy+nyest
      jwrk = lwrk-4-nxest-nyest
      knrx = 4
      knry = knrx+mx
      kndx = knry+my
      kndy = kndx+nxest
      call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),
     * wrk(lfpx),wrk(lfpy),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx),
     * iwrk(knry),iwrk(kndx),iwrk(kndy),wrk(lww),jwrk,ier)
  70  return
      end

      subroutine spalde(t,n,c,k1,x,d,ier)
c  subroutine spalde evaluates at a point x all the derivatives
c              (j-1)
c      d(j) = s     (x) , j=1,2,...,k1
c  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline
c  representation.
c
c  calling sequence:
c     call spalde(t,n,c,k1,x,d,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k1   : integer, giving the order of s(x) (order=degree+1)
c    x    : real, which contains the point where the derivatives must
c           be evaluated.
c
c  output parameters:
c    d    : array,length k1, containing the derivative values of s(x).
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    t(k1) <= x <= t(n-k1+1)
c
c  further comments:
c    if x coincides with a knot, right derivatives are computed
c    ( left derivatives if x = t(n-k1+1) ).
c
c  other subroutines required: fpader.
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k1,ier
      real x
c  ..array arguments..
      real t(n),c(n),d(k1)
c  ..local scalars..
      integer l,nk1
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      nk1 = n-k1
      if(x.lt.t(k1) .or. x.gt.t(nk1+1)) go to 300
c  search for knot interval t(l) <= x < t(l+1)
      l = k1
 100  if(x.lt.t(l+1) .or. l.eq.nk1) go to 200
      l = l+1
      go to 100
 200  if(t(l).ge.t(l+1)) go to 300
      ier = 0
c  calculate the derivatives.
      call fpader(t,n,c,k1,x,l,d)
 300  return
      end
      subroutine spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,
     * nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c  given the function values r(i,j) on the latitude-longitude grid
c  (u(i),v(j)), i=1,...,mu ; j=1,...,mv , spgrid determines a smooth
c  bicubic spline approximation on the rectangular domain 0<=u<=pi,
c  vb<=v<=ve (vb = v(1), ve=vb+2*pi).
c  this approximation s(u,v) will satisfy the properties
c
c    (1) s(0,v) = s(0,0) = dr(1)
c
c        d s(0,v)           d s(0,0)           d s(0,pi/2)
c    (2) -------- = cos(v)* -------- + sin(v)* -----------
c        d u                d u                d u
c
c                 = cos(v)*dr(2)+sin(v)*dr(3)
c                                                     vb <= v <= ve
c    (3) s(pi,v) = s(pi,0) = dr(4)
c
c        d s(pi,v)           d s(pi,0)           d s(pi,pi/2)
c    (4) -------- = cos(v)*  --------- + sin(v)* ------------
c        d u                 d u                 d u
c
c                 = cos(v)*dr(5)+sin(v)*dr(6)
c
c  and will be periodic in the variable v, i.e.
c
c         j           j
c        d s(u,vb)   d s(u,ve)
c    (5) --------- = ---------   0 <=u<= pi , j=0,1,2
c           j           j
c        d v         d v
c
c  the number of knots of s(u,v) and their position tu(i),i=1,2,...,nu;
c  tv(j),j=1,2,...,nv, is chosen automatically by the routine. the
c  smoothness of s(u,v) is achieved by minimalizing the discontinuity
c  jumps of the derivatives of the spline at the knots. the amount of
c  smoothness of s(u,v) is determined by the condition that
c  fp=sumi=1,mu(sumj=1,mv((r(i,j)-s(u(i),v(j)))**2))+(r0-s(0,v))**2
c  + (r1-s(pi,v))**2 <= s, with s a given non-negative constant.
c  the fit s(u,v) is given in its b-spline representation and can be
c  evaluated by means of routine bispev
c
c calling sequence:
c     call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,
c    *  ,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer array of dimension 3, specifying different options.
c          unchanged on exit.
c  iopt(1):on entry iopt(1) must specify whether a least-squares spline
c          (iopt(1)=-1) or a smoothing spline (iopt(1)=0 or 1) must be
c          determined.
c          if iopt(1)=0 the routine will start with an initial set of
c          knots tu(i)=0,tu(i+4)=pi,i=1,...,4;tv(i)=v(1)+(i-4)*2*pi,
c          i=1,...,8.
c          if iopt(1)=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt(1)=1 must always be immediately
c          preceded by another call with iopt(1) = 1 or iopt(1) = 0.
c  iopt(2):on entry iopt(2) must specify the requested order of conti-
c          nuity at the pole u=0.
c          if iopt(2)=0 only condition (1) must be fulfilled and
c          if iopt(2)=1 conditions (1)+(2) must be fulfilled.
c  iopt(3):on entry iopt(3) must specify the requested order of conti-
c          nuity at the pole u=pi.
c          if iopt(3)=0 only condition (3) must be fulfilled and
c          if iopt(3)=1 conditions (3)+(4) must be fulfilled.
c  ider  : integer array of dimension 4, specifying different options.
c          unchanged on exit.
c  ider(1):on entry ider(1) must specify whether (ider(1)=0 or 1) or not
c          (ider(1)=-1) there is a data value r0 at the pole u=0.
c          if ider(1)=1, r0 will be considered to be the right function
c          value, and it will be fitted exactly (s(0,v)=r0).
c          if ider(1)=0, r0 will be considered to be a data value just
c          like the other data values r(i,j).
c  ider(2):on entry ider(2) must specify whether (ider(2)=1) or not
c          (ider(2)=0) the approximation has vanishing derivatives
c          dr(2) and dr(3) at the pole u=0  (in case iopt(2)=1)
c  ider(3):on entry ider(3) must specify whether (ider(3)=0 or 1) or not
c          (ider(3)=-1) there is a data value r1 at the pole u=pi.
c          if ider(3)=1, r1 will be considered to be the right function
c          value, and it will be fitted exactly (s(pi,v)=r1).
c          if ider(3)=0, r1 will be considered to be a data value just
c          like the other data values r(i,j).
c  ider(4):on entry ider(4) must specify whether (ider(4)=1) or not
c          (ider(4)=0) the approximation has vanishing derivatives
c          dr(5) and dr(6) at the pole u=pi (in case iopt(3)=1)
c  mu    : integer. on entry mu must specify the number of grid points
c          along the u-axis. unchanged on exit.
c          mu >= 1, mu >=mumin=4-i0-i1-ider(2)-ider(4) with
c            i0=min(1,ider(1)+1), i1=min(1,ider(3)+1)
c  u     : real array of dimension at least (mu). before entry, u(i)
c          must be set to the u-co-ordinate of the i-th grid point
c          along the u-axis, for i=1,2,...,mu. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c          0 < u(i) < pi.
c  mv    : integer. on entry mv must specify the number of grid points
c          along the v-axis. mv > 3 . unchanged on exit.
c  v     : real array of dimension at least (mv). before entry, v(j)
c          must be set to the v-co-ordinate of the j-th grid point
c          along the v-axis, for j=1,2,...,mv. these values must be
c          supplied in strictly ascending order. unchanged on exit.
c          -pi <= v(1) < pi , v(mv) < v(1)+2*pi.
c  r     : real array of dimension at least (mu*mv).
c          before entry, r(mv*(i-1)+j) must be set to the data value at
c          the grid point (u(i),v(j)) for i=1,...,mu and j=1,...,mv.
c          unchanged on exit.
c  r0    : real value. on entry (if ider(1) >=0 ) r0 must specify the
c          data value at the pole u=0. unchanged on exit.
c  r1    : real value. on entry (if ider(1) >=0 ) r1 must specify the
c          data value at the pole u=pi. unchanged on exit.
c  s     : real. on entry (if iopt(1)>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nuest : integer. unchanged on exit.
c  nvest : integer. unchanged on exit.
c          on entry, nuest and nvest must specify an upper bound for the
c          number of knots required in the u- and v-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nuest >= 8, nvest >= 8.
c          in most practical situation nuest = mu/2, nvest=mv/2, will
c          be sufficient. always large enough are nuest=mu+6+iopt(2)+
c          iopt(3), nvest = mv+7, the number of knots needed for
c          interpolation (s=0). see also further comments.
c  nu    : integer.
c          unless ier=10 (in case iopt(1)>=0), nu will contain the total
c          number of knots with respect to the u-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1 is
c          used, the value of nu should be left unchanged between sub-
c          sequent calls. in case iopt(1)=-1, the value of nu should be
c          specified on entry.
c  tu    : real array of dimension at least (nuest).
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the u-variable, i.e. the position of
c          the interior knots tu(5),...,tu(nu-4) as well as the position
c          of the additional knots tu(1)=...=tu(4)=0 and tu(nu-3)=...=
c          tu(nu)=pi needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used,the values of tu(1)
c          ...,tu(nu) should be left unchanged between subsequent calls.
c          if the computation mode iopt(1)=-1 is used, the values tu(5),
c          ...tu(nu-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  nv    : integer.
c          unless ier=10 (in case iopt(1)>=0), nv will contain the total
c          number of knots with respect to the v-variable, of the spline
c          approximation returned. if the computation mode iopt(1)=1 is
c          used, the value of nv should be left unchanged between sub-
c          sequent calls. in case iopt(1) = -1, the value of nv should
c          be specified on entry.
c  tv    : real array of dimension at least (nvest).
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the v-variable, i.e. the position of
c          the interior knots tv(5),...,tv(nv-4) as well as the position
c          of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
c          tv(nv) needed for the b-spline representation.
c          if the computation mode iopt(1)=1 is used,the values of tv(1)
c          ...,tv(nv) should be left unchanged between subsequent calls.
c          if the computation mode iopt(1)=-1 is used, the values tv(5),
c          ...tv(nv-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nuest-4)*(nvest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(u,v)
c  fp    : real. unless ier=10, fp contains the sum of squared
c          residuals of the spline approximation returned.
c  wrk   : real array of dimension (lwrk). used as workspace.
c          if the computation mode iopt(1)=1 is used the values of
c          wrk(1),..,wrk(12) should be left unchanged between subsequent
c          calls.
c  lwrk  : integer. on entry lwrk must specify the actual dimension of
c          the array wrk as declared in the calling (sub)program.
c          lwrk must not be too small.
c           lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+q
c           where q is the larger of (mv+nvest) and nuest.
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c          if the computation mode iopt(1)=1 is used the values of
c          iwrk(1),.,iwrk(5) should be left unchanged between subsequent
c          calls.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= 5+mu+mv+nuest+nvest.
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the least-squares
c            constrained polynomial. in this extreme case fp gives the
c            upper bound for the smoothing factor s.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nuest and
c            nvest.
c            probably causes : nuest or nvest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the least-squares spline
c            according to the current set of knots. the parameter fp
c            gives the corresponding sum of squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small.
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            sum of squared residuals does not satisfy the condition
c            abs(fp-s)/s < tol.
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt(1)<=1, 0<=iopt(2)<=1, 0<=iopt(3)<=1,
c            -1<=ider(1)<=1, 0<=ider(2)<=1, ider(2)=0 if iopt(2)=0.
c            -1<=ider(3)<=1, 0<=ider(4)<=1, ider(4)=0 if iopt(3)=0.
c            mu >= mumin (see above), mv >= 4, nuest >=8, nvest >= 8,
c            kwrk>=5+mu+mv+nuest+nvest,
c            lwrk >= 12+nuest*(mv+nvest+3)+nvest*24+4*mu+8*mv+
c             max(nuest,mv+nvest)
c            0< u(i-1)<u(i)< pi,i=2,..,mu,
c            -pi<=v(1)< pi, v(1)<v(i-1)<v(i)<v(1)+2*pi, i=3,...,mv
c            if iopt(1)=-1: 8<=nu<=min(nuest,mu+6+iopt(2)+iopt(3))
c                           0<tu(5)<tu(6)<...<tu(nu-4)< pi
c                           8<=nv<=min(nvest,mv+7)
c                           v(1)<tv(5)<tv(6)<...<tv(nv-4)<v(1)+2*pi
c                    the schoenberg-whitney conditions, i.e. there must
c                    be subset of grid co-ordinates uu(p) and vv(q) such
c                    that   tu(p) < uu(p) < tu(p+4) ,p=1,...,nu-4
c                     (iopt(2)=1 and iopt(3)=1 also count for a uu-value
c                           tv(q) < vv(q) < tv(q+4) ,q=1,...,nv-4
c                     (vv(q) is either a value v(j) or v(j)+2*pi)
c            if iopt(1)>=0: s>=0
c                       if s=0: nuest>=mu+6+iopt(2)+iopt(3), nvest>=mv+7
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c
c further comments:
c   spgrid does not allow individual weighting of the data-values.
c   so, if these were determined to widely different accuracies, then
c   perhaps the general data set routine sphere should rather be used
c   in spite of efficiency.
c   by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the constrained least-squares polynomial(degrees 3,0)if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the accuracy of the data values.
c   if the user has an idea of the statistical errors on the data, he
c   can also find a proper estimate for s. for, by assuming that, if he
c   specifies the right s, spgrid will return a spline s(u,v) which
c   exactly reproduces the function underlying the data he can evaluate
c   the sum((r(i,j)-s(u(i),v(j)))**2) to find a good estimate for this s
c   for example, if he knows that the statistical errors on his r(i,j)-
c   values is not greater than 0.1, he may expect that a good s should
c   have a value not larger than mu*mv*(0.1)**2.
c   if nothing is known about the statistical error in r(i,j), s must
c   be determined by trial and error, taking account of the comments
c   above. the best is then to start with a very large value of s (to
c   determine the least-squares polynomial and the corresponding upper
c   bound fp0 for s) and then to progressively decrease the value of s
c   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
c   and more carefully as the approximation shows more detail) to
c   obtain closer fits.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt(1)=0.
c   if iopt(1) = 1 the program will continue with the knots found at
c   the last call of the routine. this will save a lot of computation
c   time if spgrid is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt(1) = 1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt(1)=1,the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   spgrid once more with the chosen value for s but now with iopt(1)=0.
c   indeed, spgrid may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nuest and
c   nvest. indeed, if at a certain stage in spgrid the number of knots
c   in one direction (say nu) has reached the value of its upper bound
c   (nuest), then from that moment on all subsequent knots are added
c   in the other (v) direction. this may indicate that the value of
c   nuest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nuest=8 (the lowest allowable value for
c   nuest), the user can indicate that he wants an approximation which
c   is a simple cubic polynomial in the variable u.
c
c  other subroutines required:
c    fpspgr,fpchec,fpchep,fpknot,fpopsp,fprati,fpgrsp,fpsysy,fpback,
c    fpbacp,fpbspl,fpcyt1,fpcyt2,fpdisc,fpgivs,fprota
c
c  references:
c   dierckx p. : fast algorithms for smoothing data over a disc or a
c                sphere using tensor product splines, in "algorithms
c                for approximation", ed. j.c.mason and m.g.cox,
c                clarendon press oxford, 1987, pp. 51-65
c   dierckx p. : fast algorithms for smoothing data over a disc or a
c                sphere using tensor product splines, report tw73, dept.
c                computer science,k.u.leuven, 1985.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : july 1985
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real r0,r1,s,fp
      integer mu,mv,nuest,nvest,nu,nv,lwrk,kwrk,ier
c  ..array arguments..
      integer iopt(3),ider(4),iwrk(kwrk)
      real u(mu),v(mv),r(mu*mv),c((nuest-4)*(nvest-4)),tu(nuest),
     * tv(nvest),wrk(lwrk)
c  ..local scalars..
      real per,pi,tol,uu,ve,rmax,rmin,one,half,rn,rb,re
      integer i,i1,i2,j,jwrk,j1,j2,kndu,kndv,knru,knrv,kwest,l,
     * ldr,lfpu,lfpv,lwest,lww,m,maxit,mumin,muu,nc
c  ..function references..
      real atan2
      integer max0
c  ..subroutine references..
c    fpchec,fpchep,fpspgr
c  ..
c  set constants
      one = 1
      half = 0.5e0
      pi = atan2(0.,-one)
      per = pi+pi
      ve = v(1)+per
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations, a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt(1).lt.(-1) .or. iopt(1).gt.1) go to 200
      if(iopt(2).lt.0 .or. iopt(2).gt.1) go to 200
      if(iopt(3).lt.0 .or. iopt(3).gt.1) go to 200
      if(ider(1).lt.(-1) .or. ider(1).gt.1) go to 200
      if(ider(2).lt.0 .or. ider(2).gt.1) go to 200
      if(ider(2).eq.1 .and. iopt(2).eq.0) go to 200
      if(ider(3).lt.(-1) .or. ider(3).gt.1) go to 200
      if(ider(4).lt.0 .or. ider(4).gt.1) go to 200
      if(ider(4).eq.1 .and. iopt(3).eq.0) go to 200
      mumin = 4
      if(ider(1).ge.0) mumin = mumin-1
      if(iopt(2).eq.1 .and. ider(2).eq.1) mumin = mumin-1
      if(ider(3).ge.0) mumin = mumin-1
      if(iopt(3).eq.1 .and. ider(4).eq.1) mumin = mumin-1
      if(mumin.eq.0) mumin = 1
      if(mu.lt.mumin .or. mv.lt.4) go to 200
      if(nuest.lt.8 .or. nvest.lt.8) go to 200
      m = mu*mv
      nc = (nuest-4)*(nvest-4)
      lwest = 12+nuest*(mv+nvest+3)+24*nvest+4*mu+8*mv+
     * max0(nuest,mv+nvest)
      kwest = 5+mu+mv+nuest+nvest
      if(lwrk.lt.lwest .or. kwrk.lt.kwest) go to 200
      if(u(1).le.0. .or. u(mu).ge.pi) go to 200
      if(mu.eq.1) go to 30
      do 20 i=2,mu
        if(u(i-1).ge.u(i)) go to 200
  20  continue
  30  if(v(1).lt. (-pi) .or. v(1).ge.pi ) go to 200
      if(v(mv).ge.v(1)+per) go to 200
      do 40 i=2,mv
        if(v(i-1).ge.v(i)) go to 200
  40  continue
      if(iopt(1).gt.0) go to 140
c  if not given, we compute an estimate for r0.
      rn = mv
      if(ider(1).lt.0) go to 45
      rb = r0
      go to 55
  45  rb = 0.
      do 50 i=1,mv
         rb = rb+r(i)
  50  continue
      rb = rb/rn
c  if not given, we compute an estimate for r1.
  55  if(ider(3).lt.0) go to 60
      re = r1
      go to 70
  60  re = 0.
      j = m
      do 65 i=1,mv
         re = re+r(j)
         j = j-1
  65  continue
      re = re/rn
c  we determine the range of r-values.
  70  rmin = rb
      rmax = re
      do 80 i=1,m
         if(r(i).lt.rmin) rmin = r(i)
         if(r(i).gt.rmax) rmax = r(i)
  80  continue
      wrk(5) = rb
      wrk(6) = 0.
      wrk(7) = 0.
      wrk(8) = re
      wrk(9) = 0.
      wrk(10) = 0.
      wrk(11) = rmax -rmin
      wrk(12) = wrk(11)
      iwrk(4) = mu
      iwrk(5) = mu
      if(iopt(1).eq.0) go to 140
      if(nu.lt.8 .or. nu.gt.nuest) go to 200
      if(nv.lt.11 .or. nv.gt.nvest) go to 200
      j = nu
      do 90 i=1,4
        tu(i) = 0.
        tu(j) = pi
        j = j-1
  90  continue
      l = 13
      wrk(l) = 0.
      if(iopt(2).eq.0) go to 100
      l = l+1
      uu = u(1)
      if(uu.gt.tu(5)) uu = tu(5)
      wrk(l) = uu*half
 100  do 110 i=1,mu
        l = l+1
        wrk(l) = u(i)
 110  continue
      if(iopt(3).eq.0) go to 120
      l = l+1
      uu = u(mu)
      if(uu.lt.tu(nu-4)) uu = tu(nu-4)
      wrk(l) = uu+(pi-uu)*half
 120  l = l+1
      wrk(l) = pi
      muu = l-12
      call fpchec(wrk(13),muu,tu,nu,3,ier)
      if(ier.ne.0) go to 200
      j1 = 4
      tv(j1) = v(1)
      i1 = nv-3
      tv(i1) = ve
      j2 = j1
      i2 = i1
      do 130 i=1,3
        i1 = i1+1
        i2 = i2-1
        j1 = j1+1
        j2 = j2-1
        tv(j2) = tv(i2)-per
        tv(i1) = tv(j1)+per
 130  continue
      l = 13
      do 135 i=1,mv
        wrk(l) = v(i)
        l = l+1
 135  continue
      wrk(l) = ve
      call fpchep(wrk(13),mv+1,tv,nv,3,ier)
      if(ier) 200,150,200
 140  if(s.lt.0.) go to 200
      if(s.eq.0. .and. (nuest.lt.(mu+6+iopt(2)+iopt(3)) .or.
     * nvest.lt.(mv+7)) ) go to 200
c  we partition the working space and determine the spline approximation
 150  ldr = 5
      lfpu = 13
      lfpv = lfpu+nuest
      lww = lfpv+nvest
      jwrk = lwrk-12-nuest-nvest
      knru = 6
      knrv = knru+mu
      kndu = knrv+mv
      kndv = kndu+nuest
      call fpspgr(iopt,ider,u,mu,v,mv,r,m,rb,re,s,nuest,nvest,tol,maxit,
     * nc,nu,tu,nv,tv,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),wrk(lfpu),
     * wrk(lfpv),wrk(ldr),wrk(11),iwrk(1),iwrk(2),iwrk(3),iwrk(4),
     * iwrk(5),iwrk(knru),iwrk(knrv),iwrk(kndu),iwrk(kndv),wrk(lww),
     * jwrk,ier)
 200  return
      end
      subroutine sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,
     *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c  subroutine sphere determines a smooth bicubic spherical spline
c  approximation s(teta,phi), 0 <= teta <= pi ; 0 <= phi <= 2*pi
c  to a given set of data points (teta(i),phi(i),r(i)),i=1,2,...,m.
c  such a spline has the following specific properties
c
c    (1) s(0,phi)  = constant   0 <=phi<= 2*pi.
c
c    (2) s(pi,phi) = constant   0 <=phi<= 2*pi
c
c         j             j
c        d s(teta,0)   d s(teta,2*pi)
c    (3) ----------- = ------------   0 <=teta<=pi, j=0,1,2
c             j             j
c        d phi         d phi
c
c        d s(0,phi)    d s(0,0)             d s(0,pi/2)
c    (4) ----------  = -------- *cos(phi) + ----------- *sin(phi)
c        d teta        d teta               d teta
c
c        d s(pi,phi)   d s(pi,0)            d s(pi,pi/2)
c    (5) ----------- = ---------*cos(phi) + ------------*sin(phi)
c        d teta        d teta               d teta
c
c  if iopt =-1 sphere calculates a weighted least-squares spherical
c  spline according to a given set of knots in teta- and phi- direction.
c  if iopt >=0, the number of knots in each direction and their position
c  tt(j),j=1,2,...,nt ; tp(j),j=1,2,...,np are chosen automatically by
c  the routine. the smoothness of s(teta,phi) is then achieved by mini-
c  malizing the discontinuity jumps of the derivatives of the spline
c  at the knots. the amount of smoothness of s(teta,phi) is determined
c  by the condition that fp = sum((w(i)*(r(i)-s(teta(i),phi(i))))**2)
c  be <= s, with s a given non-negative constant.
c  the spherical spline is given in the standard b-spline representation
c  of bicubic splines and can be evaluated by means of subroutine bispev
c
c calling sequence:
c     call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,
c    *  nt,tt,np,tp,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a weighted
c          least-squares spherical spline (iopt=-1) or a smoothing
c          spherical spline (iopt=0 or 1) must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tt(i)=0,tt(i+4)=pi,i=1,...,4;tp(i)=0,tp(i+4)=2*pi,i=1,...,4.
c          if iopt=1 the routine will continue with the set of knots
c          found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c  m     : integer. on entry m must specify the number of data points.
c          m >= 2. unchanged on exit.
c  teta  : real array of dimension at least (m).
c  phi   : real array of dimension at least (m).
c  r     : real array of dimension at least (m).
c          before entry,teta(i),phi(i),r(i) must be set to the spherical
c          co-ordinates of the i-th data point, for i=1,...,m.the order
c          of the data points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  ntest : integer. unchanged on exit.
c  npest : integer. unchanged on exit.
c          on entry, ntest and npest must specify an upper bound for the
c          number of knots required in the teta- and phi-directions.
c          these numbers will also determine the storage space needed by
c          the routine. ntest >= 8, npest >= 8.
c          in most practical situation ntest = npest = 8+sqrt(m/2) will
c          be sufficient. see also further comments.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nt    : integer.
c          unless ier=10 (in case iopt >=0), nt will contain the total
c          number of knots with respect to the teta-variable, of the
c          spline approximation returned. if the computation mode iopt=1
c          is used, the value of nt should be left unchanged between
c          subsequent calls.
c          in case iopt=-1, the value of nt should be specified on entry
c  tt    : real array of dimension at least ntest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the teta-variable, i.e. the position
c          of the interior knots tt(5),...,tt(nt-4) as well as the
c          position of the additional knots tt(1)=...=tt(4)=0 and
c          tt(nt-3)=...=tt(nt)=pi needed for the b-spline representation
c          if the computation mode iopt=1 is used, the values of tt(1),
c          ...,tt(nt) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tt(5),
c          ...tt(nt-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  np    : integer.
c          unless ier=10 (in case iopt >=0), np will contain the total
c          number of knots with respect to the phi-variable, of the
c          spline approximation returned. if the computation mode iopt=1
c          is used, the value of np should be left unchanged between
c          subsequent calls.
c          in case iopt=-1, the value of np (>=9) should be specified
c          on entry.
c  tp    : real array of dimension at least npest.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the phi-variable, i.e. the position of
c          the interior knots tp(5),...,tp(np-4) as well as the position
c          of the additional knots tp(1),...,tp(4) and tp(np-3),...,
c          tp(np) needed for the b-spline representation.
c          if the computation mode iopt=1 is used, the values of tp(1),
c          ...,tp(np) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tp(5),
c          ...tp(np-4) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (ntest-4)*(npest-4).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(teta,phi).
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt=1 is used the value of wrk1(1)
c          should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
c          values d(i)/max(d(i)),i=1,...,ncof=6+(np-7)*(nt-8)
c          with d(i) the i-th diagonal element of the reduced triangular
c          matrix for calculating the b-spline coefficients. it includes
c          those elements whose square is less than eps,which are treat-
c          ed as 0 in the case of presumed rank deficiency (ier<-2).
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            u = ntest-7, v = npest-7, then
c          lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper bound  for lwrk2 = 48+21*v+7*u*v+
c          4*(u-1)*v**2 where u,v are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(ntest-7)*(npest-7).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is a spherical
c            interpolating spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares constrained polynomial . in this extreme case
c            fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            6+(nt-8)*(np-7)+ier, is large the results may be inaccurate
c            they could also seriously depend on the value of eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters ntest and
c            npest.
c            probably causes : ntest or npest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small or badly chosen eps.
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=4  : error. no more knots can be added because the dimension
c            of the spherical spline 6+(nt-8)*(np-7) already exceeds
c            the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            spherical spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1,  m>=2, ntest>=8 ,npest >=8, 0<eps<1,
c            0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
c            lwrk1 >= 185+52*v+10*u+14*u*v+8*(u-1)*v**2+8*m
c            kwrk >= m+(ntest-7)*(npest-7)
c            if iopt=-1: 8<=nt<=ntest , 9<=np<=npest
c                        0<tt(5)<tt(6)<...<tt(nt-4)<pi
c                        0<tp(5)<tp(6)<...<tp(np-4)<2*pi
c            if iopt>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nt,
c            np,tt,tp,wrk1, and having adjusted the value of lwrk2 and
c            the dimension of the array wrk2 accordingly, the user can
c            continue at the point the program was left, by calling
c            sphere with iopt=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the constrained weighted least-squares polynomial if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   r(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in r(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to choose s very small is strongly discouraged. this considerably
c   increases computation time and memory requirements. it may also
c   cause rank-deficiency (ier<-2) and endager numerical stability.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if sphere is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   sphere once more with the selected value for s but now with iopt=0.
c   indeed, sphere may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds ntest and
c   npest. indeed, if at a certain stage in sphere the number of knots
c   in one direction (say nt) has reached the value of its upper bound
c   (ntest), then from that moment on all subsequent knots are added
c   in the other (phi) direction. this may indicate that the value of
c   ntest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting ntest=8 (the lowest allowable value for
c   ntest), the user can indicate that he wants an approximation which
c   is a cubic polynomial in the variable teta.
c
c  other subroutines required:
c    fpback,fpbspl,fpsphe,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
c    fprpsp
c
c  references:
c   dierckx p. : algorithms for smoothing data on the sphere with tensor
c                product splines, computing 32 (1984) 319-342.
c   dierckx p. : algorithms for smoothing data on the sphere with tensor
c                product splines, report tw62, dept. computer science,
c                k.u.leuven, 1983.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : july 1983
c  latest update : march 1989
c
c  ..
c  ..scalar arguments..
      real s,eps,fp
      integer iopt,m,ntest,npest,nt,np,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real teta(m),phi(m),r(m),w(m),tt(ntest),tp(npest),
     * c((ntest-4)*(npest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
c  ..local scalars..
      real tol,pi,pi2,one
      integer i,ib1,ib3,ki,kn,kwest,la,lbt,lcc,lcs,lro,j
     * lbp,lco,lf,lff,lfp,lh,lq,lst,lsp,lwest,maxit,ncest,ncc,ntt,
     * npp,nreg,nrint,ncof,nt4,np4
c  ..function references..
      real atan
c  ..subroutine references..
c    fpsphe
c  ..
c  set constants
      one = 0.1e+01
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 80
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 80
      if(m.lt.2) go to 80
      if(ntest.lt.8 .or. npest.lt.8) go to 80
      nt4 = ntest-4
      np4 = npest-4
      ncest = nt4*np4
      ntt = ntest-7
      npp = npest-7
      ncc = 6+npp*(ntt-1)
      nrint = ntt+npp
      nreg = ntt*npp
      ncof = 6+3*npp
      ib1 = 4*npp
      ib3 = ib1+3
      if(ncof.gt.ib1) ib1 = ncof
      if(ncof.gt.ib3) ib3 = ncof
      lwest = 185+52*npp+10*ntt+14*ntt*npp+8*(m+(ntt-1)*npp**2)
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 80
      if(iopt.gt.0) go to 60
      pi = atan(one)*4
      pi2 = pi+pi
      do 20 i=1,m
        if(w(i).le.0.) go to 80
        if(teta(i).lt.0. .or. teta(i).gt.pi) go to 80
        if(phi(i) .lt.0. .or. phi(i).gt.pi2) go to 80
  20  continue
      if(iopt.eq.0) go to 60
      ntt = nt-8
      if(ntt.lt.0 .or. nt.gt.ntest) go to 80
      if(ntt.eq.0) go to 40
      tt(4) = 0.
      do 30 i=1,ntt
         j = i+4
         if(tt(j).le.tt(j-1) .or. tt(j).ge.pi) go to 80
  30  continue
  40  npp = np-8
      if(npp.lt.1 .or. np.gt.npest) go to 80
      tp(4) = 0.
      do 50 i=1,npp
         j = i+4
         if(tp(j).le.tp(j-1) .or. tp(j).ge.pi2) go to 80
  50  continue
      go to 70
  60  if(s.lt.0.) go to 80
  70  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbt = lh+ib3
      lbp = lbt+5*ntest
      lro = lbp+5*npest
      lcc = lro+npest
      lcs = lcc+npest
      lst = lcs+npest
      lsp = lst+m*4
      call fpsphe(iopt,m,teta,phi,r,w,s,ntest,npest,eps,tol,maxit,
     * ib1,ib3,ncest,ncc,nrint,nreg,nt,tt,np,tp,c,fp,wrk1(1),wrk1(lfp),
     * wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),wrk1(lcs),
     * wrk1(la),wrk1(lq),wrk1(lbt),wrk1(lbp),wrk1(lst),wrk1(lsp),
     * wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  80  return
      end

      subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
c  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
c  the derivative of order nu of a spline s(x) of degree k,given in
c  its b-spline representation.
c
c  calling sequence:
c     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    nu   : integer, specifying the order of the derivative. 0<=nu<=k
c    x    : array,length m, which contains the points where the deriv-
c           ative of s(x) must be evaluated.
c    m    : integer, giving the number of points where the derivative
c           of s(x) must be evaluated
c    wrk  : real array of dimension n. used as working space.
c
c  output parameters:
c    y    : array,length m, giving the value of the derivative of s(x)
c           at the different points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    0 <= nu <= k
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,nu,m,ier
c  ..array arguments..
      real t(n),c(n),x(m),y(m),wrk(n)
c  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      real ak,arg,fac,sp,tb,te
c  ..local arrays ..
      real h(6)
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
      if(m-1) 200,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 200
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
c  the derivative of order nu of a spline of degree k is a spline of
c  degree k-nu,the b-spline coefficients wrk(i) of which can be found
c  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
c  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
         arg = x(i)
  70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
         l = l+1
         j = j+1
         go to 70
  80     y(i) = wrk(j)
  90  continue
      go to 200
 100  l = k1
      l1 = l+1
      k2 = k1-nu
c  main loop for the different points.
      do 180 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
c  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
c  find the value of the derivative at x=arg.
        sp = 0.
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end
      subroutine splev(t,n,c,k,x,y,m,ier)
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      real t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      real arg,sp,tb,te
c  ..local array..
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
      real function splint(t,n,c,k,a,b,wrk)
c  function splint calculates the integral of a spline function s(x)
c  of degree k, which is given in its normalized b-spline representation
c
c  calling sequence:
c     aint = splint(t,n,c,k,a,b,wrk)
c
c  input parameters:
c    t    : array,length n,which contains the position of the knots
c           of s(x).
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, containing the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    a,b  : real values, containing the end points of the integration
c           interval. s(x) is considered to be identically zero outside
c           the interval (t(k+1),t(n-k)).
c
c  output parameter:
c    aint : real, containing the integral of s(x) between a and b.
c    wrk  : real array, length n.  used as working space
c           on output, wrk will contain the integrals of the normalized
c           b-splines defined on the set of knots.
c
c  other subroutines required: fpintb.
c
c  references :
c    gaffney p.w. : the calculation of indefinite integrals of b-splines
c                   j. inst. maths applics 17 (1976) 37-41.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      real a,b
      integer n,k
c  ..array arguments..
      real t(n),c(n),wrk(n)
c  ..local scalars..
      integer i,nk1
c  ..
      nk1 = n-k-1
c  calculate the integrals wrk(i) of the normalized b-splines
c  ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)
c  calculate the integral of s(x).
      splint = 0.
      do 10 i=1,nk1
        splint = splint+c(i)*wrk(i)
  10  continue
      return
      end
      subroutine sproot(t,n,c,zero,mest,m,ier)
c  subroutine sproot finds the zeros of a cubic spline s(x),which is
c  given in its normalized b-spline representation.
c
c  calling sequence:
c     call sproot(t,n,c,zero,mest,m,ier)
c
c  input parameters:
c    t    : real array,length n, containing the knots of s(x).
c    n    : integer, containing the number of knots.  n>=8
c    c    : real array,length n, containing the b-spline coefficients.
c    mest : integer, specifying the dimension of array zero.
c
c  output parameters:
c    zero : real array,lenth mest, containing the zeros of s(x).
c    m    : integer,giving the number of zeros.
c    ier  : error flag:
c      ier = 0: normal return.
c      ier = 1: the number of zeros exceeds mest.
c      ier =10: invalid input data (see restrictions).
c
c  other subroutines required: fpcuro
c
c  restrictions:
c    1) n>= 8.
c    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
c       t(1) <= t(2) <= t(3) <= t(4)
c       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c ..
c ..scalar arguments..
      integer n,mest,m,ier
c  ..array arguments..
      real t(n),c(n),zero(mest)
c  ..local scalars..
      integer i,j,j1,l,n4
      real ah,a0,a1,a2,a3,bh,b0,b1,c1,c2,c3,c4,c5,d4,d5,h1,h2,
     * three,two,t1,t2,t3,t4,t5,zz
      logical z0,z1,z2,z3,z4,nz0,nz1,nz2,nz3,nz4
c  ..local array..
      real y(3)
c  ..
c  set some constants
      two = 0.2e+01
      three = 0.3e+01
c  before starting computations a data check is made. if the input data
c  are invalid, control is immediately repassed to the calling program.
      n4 = n-4
      ier = 10
      if(n.lt.8) go to 800
      j = n
      do 10 i=1,3
        if(t(i).gt.t(i+1)) go to 800
        if(t(j).lt.t(j-1)) go to 800
        j = j-1
  10  continue
      do 20 i=4,n4
        if(t(i).ge.t(i+1)) go to 800
  20  continue
c  the problem considered reduces to finding the zeros of the cubic
c  polynomials pl(x) which define the cubic spline in each knot
c  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
c  the condition that it belongs to the knot interval.
c  the cubic polynomial pl(x) is determined by computing s(t(l)),
c  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
c  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
c  splines and their derivatives, the value of s(t(l)) and s'(t(l))
c  is already known from the foregoing knot interval.
      ier = 0
c  evaluate some constants for the first knot interval
      h1 = t(4)-t(3)
      h2 = t(5)-t(4)
      t1 = t(4)-t(2)
      t2 = t(5)-t(3)
      t3 = t(6)-t(4)
      t4 = t(5)-t(2)
      t5 = t(6)-t(3)
c  calculate a0 = s(t(4)) and ah = s'(t(4)).
      c1 = c(1)
      c2 = c(2)
      c3 = c(3)
      c4 = (c2-c1)/t4
      c5 = (c3-c2)/t5
      d4 = (h2*c1+t1*c2)/t4
      d5 = (t3*c2+h1*c3)/t5
      a0 = (h2*d4+h1*d5)/t2
      ah = three*(h2*c4+h1*c5)/t2
      z1 = .true.
      if(ah.lt.0.) z1 = .false.
      nz1 = .not.z1
      m = 0
c  main loop for the different knot intervals.
      do 300 l=4,n4
c  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2
        h2 = t(l+2)-t(l+1)
        t1 = t2
        t2 = t3
        t3 = t(l+3)-t(l+1)
        t4 = t5
        t5 = t(l+3)-t(l)
c  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2
        c2 = c3
        c3 = c(l)
        c4 = c5
        c5 = (c3-c2)/t5
        d4 = (h2*c1+t1*c2)/t4
        d5 = (h1*c3+t3*c2)/t5
        b0 = (h2*d4+h1*d5)/t2
        bh = three*(h2*c4+h1*c5)/t2
c  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
c  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
        a1 = ah*h1
        b1 = bh*h1
        a2 = three*(b0-a0)-b1-two*a1
        a3 = two*(a0-b0)+b1+a1
c  test whether or not pl(x) could have a zero in the range
c  t(l) <= x <= t(l+1).
        z3 = .true.
        if(b1.lt.0.) z3 = .false.
        nz3 = .not.z3
        if(a0*b0.le.0.) go to 100
        z0 = .true.
        if(a0.lt.0.) z0 = .false.
        nz0 = .not.z0
        z2 = .true.
        if(a2.lt.0.) z2 = .false.
        nz2 = .not.z2
        z4 = .true.
        if(3.0*a3+a2.lt.0.) z4 = .false.
        nz4 = .not.z4
        if(.not.((z0.and.(nz1.and.(z3.or.z2.and.nz4).or.nz2.and.
     * z3.and.z4).or.nz0.and.(z1.and.(nz3.or.nz2.and.z4).or.z2.and.
     * nz3.and.nz4))))go to 200
c  find the zeros of ql(y).
 100    call fpcuro(a3,a2,a1,a0,y,j)
        if(j.eq.0) go to 200
c  find which zeros of pl(x) are zeros of s(x).
        do 150 i=1,j
          if(y(i).lt.0. .or. y(i).gt.1.0) go to 150
c  test whether the number of zeros of s(x) exceeds mest.
          if(m.ge.mest) go to 700
          m = m+1
          zero(m) = t(l)+h1*y(i)
 150    continue
 200    a0 = b0
        ah = bh
        z1 = z3
        nz1 = nz3
 300  continue
c  the zeros of s(x) are arranged in increasing order.
      if(m.lt.2) go to 800
      do 400 i=2,m
        j = i
 350    j1 = j-1
        if(j1.eq.0) go to 400
        if(zero(j).ge.zero(j1)) go to 400
        zz = zero(j)
        zero(j) = zero(j1)
        zero(j1) = zz
        j = j1
        go to 350
 400  continue
      j = m
      m = 1
      do 500 i=2,j
        if(zero(i).eq.zero(m)) go to 500
        m = m+1
        zero(m) = zero(i)
 500  continue
      go to 800
 700  ier = 1
 800  return
      end
      subroutine surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,
     * iwrk,kwrk,ier)
c  subroutine surev evaluates on a grid (u(i),v(j)),i=1,...,mu; j=1,...
c  ,mv a bicubic spline surface of dimension idim, given in the
c  b-spline representation.
c
c  calling sequence:
c     call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   idim  : integer, specifying the dimension of the spline surface.
c   tu    : real array, length nu, which contains the position of the
c           knots in the u-direction.
c   nu    : integer, giving the total number of knots in the u-direction
c   tv    : real array, length nv, which contains the position of the
c           knots in the v-direction.
c   nv    : integer, giving the total number of knots in the v-direction
c   c     : real array, length (nu-4)*(nv-4)*idim, which contains the
c           b-spline coefficients.
c   u     : real array of dimension (mu).
c           before entry u(i) must be set to the u co-ordinate of the
c           i-th grid point along the u-axis.
c           tu(4)<=u(i-1)<=u(i)<=tu(nu-3), i=2,...,mu.
c   mu    : on entry mu must specify the number of grid points along
c           the u-axis. mu >=1.
c   v     : real array of dimension (mv).
c           before entry v(j) must be set to the v co-ordinate of the
c           j-th grid point along the v-axis.
c           tv(4)<=v(j-1)<=v(j)<=tv(nv-3), j=2,...,mv.
c   mv    : on entry mv must specify the number of grid points along
c           the v-axis. mv >=1.
c   mf    : on entry, mf must specify the dimension of the array f.
c           mf >= mu*mv*idim
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= 4*(mu+mv)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mu+mv.
c
c  output parameters:
c   f     : real array of dimension (mf).
c           on succesful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the
c           l-th co-ordinate of the bicubic spline surface at the
c           point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mu >=1, mv >=1, lwrk>=4*(mu+mv), kwrk>=mu+mv , mf>=mu*mv*idim
c   tu(4) <= u(i-1) <= u(i) <= tu(nu-3), i=2,...,mu
c   tv(4) <= v(j-1) <= v(j) <= tv(nv-3), j=2,...,mv
c
c  other subroutines required:
c    fpsuev,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer idim,nu,nv,mu,mv,mf,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv),f(mf),
     * wrk(lwrk)
c  ..local scalars..
      integer i,muv
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(mf.lt.mu*mv*idim) go to 100
      muv = mu+mv
      if(lwrk.lt.4*muv) go to 100
      if(kwrk.lt.muv) go to 100
      if(mu-1) 100,30,10
  10  do 20 i=2,mu
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mv-1) 100,60,40
  40  do 50 i=2,mv
        if(v(i).lt.v(i-1)) go to 100
  50  continue
  60  ier = 0
      call fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wrk(1),wrk(4*mu+1),
     * iwrk(1),iwrk(mu+1))
 100  return
      end
      subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c given the set of data points (x(i),y(i),z(i)) and the set of positive
c numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 surfit calculates the weighted least-squares spline
c according to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
c ative constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
c     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
c    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
c  iopt  : integer flag. on entry iopt must specify whether a weighted
c          least-squares spline (iopt=-1) or a smoothing spline (iopt=0
c          or 1) must be determined.
c          if iopt=0 the routine will start with an initial set of knots
c          tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
c          1,...,ky+1. if iopt=1 the routine will continue with the set
c          of knots found at the last call of the routine.
c          attention: a call with iopt=1 must always be immediately pre-
c                     ceded by another call with iopt=1 or iopt=0.
c          unchanged on exit.
c  m     : integer. on entry m must specify the number of data points.
c          m >= (kx+1)*(ky+1). unchanged on exit.
c  x     : real array of dimension at least (m).
c  y     : real array of dimension at least (m).
c  z     : real array of dimension at least (m).
c          before entry, x(i),y(i),z(i) must be set to the co-ordinates
c          of the i-th data point, for i=1,...,m. the order of the data
c          points is immaterial. unchanged on exit.
c  w     : real array of dimension at least (m). before entry, w(i) must
c          be set to the i-th value in the set of weights. the w(i) must
c          be strictly positive. unchanged on exit.
c  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
c  yb,ye   aries of the rectangular approximation domain.
c          xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
c  kx,ky : integer values. on entry kx and ky must specify the degrees
c          of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
c          (kx=ky=3) splines. unchanged on exit.
c  s     : real. on entry (in case iopt>=0) s must specify the smoothing
c          factor. s >=0. unchanged on exit.
c          for advice on the choice of s see further comments
c  nxest : integer. unchanged on exit.
c  nyest : integer. unchanged on exit.
c          on entry, nxest and nyest must specify an upper bound for the
c          number of knots required in the x- and y-directions respect.
c          these numbers will also determine the storage space needed by
c          the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
c          in most practical situation nxest = kx+1+sqrt(m/2), nyest =
c          ky+1+sqrt(m/2) will be sufficient. see also further comments.
c  nmax  : integer. on entry nmax must specify the actual dimension of
c          the arrays tx and ty. nmax >= nxest, nmax >=nyest.
c          unchanged on exit.
c  eps   : real.
c          on entry, eps must specify a threshold for determining the
c          effective rank of an over-determined linear system of equat-
c          ions. 0 < eps < 1.  if the number of decimal digits in the
c          computer representation of a real number is q, then 10**(-q)
c          is a suitable value for eps in most practical applications.
c          unchanged on exit.
c  nx    : integer.
c          unless ier=10 (in case iopt >=0), nx will contain the total
c          number of knots with respect to the x-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of nx should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of nx should be specified on entry
c  tx    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the x-variable, i.e. the position of
c          the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
c          position of the additional knots tx(1)=...=tx(kx+1)=xb and
c          tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of tx(1),
c          ...,tx(nx) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values tx(kx+2),
c          ...tx(nx-kx-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  ny    : integer.
c          unless ier=10 (in case iopt >=0), ny will contain the total
c          number of knots with respect to the y-variable, of the spline
c          approximation returned. if the computation mode iopt=1 is
c          used, the value of ny should be left unchanged between sub-
c          sequent calls.
c          in case iopt=-1, the value of ny should be specified on entry
c  ty    : real array of dimension nmax.
c          on succesful exit, this array will contain the knots of the
c          spline with respect to the y-variable, i.e. the position of
c          the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
c          position of the additional knots ty(1)=...=ty(ky+1)=yb and
c          ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
c          if the computation mode iopt=1 is used, the values of ty(1),
c          ...,ty(ny) should be left unchanged between subsequent calls.
c          if the computation mode iopt=-1 is used, the values ty(ky+2),
c          ...ty(ny-ky-1) must be supplied by the user, before entry.
c          see also the restrictions (ier=10).
c  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
c          on succesful exit, c contains the coefficients of the spline
c          approximation s(x,y)
c  fp    : real. unless ier=10, fp contains the weighted sum of
c          squared residuals of the spline approximation returned.
c  wrk1  : real array of dimension (lwrk1). used as workspace.
c          if the computation mode iopt=1 is used the value of wrk1(1)
c          should be left unchanged between subsequent calls.
c          on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
c          contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
c          with d(i) the i-th diagonal element of the reduced triangular
c          matrix for calculating the b-spline coefficients. it includes
c          those elements whose square is less than eps,which are treat-
c          ed as 0 in the case of presumed rank deficiency (ier<-2).
c  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
c          the array wrk1 as declared in the calling (sub)program.
c          lwrk1 must not be too small. let
c            u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
c            ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
c            if(bx.le.by) b1 = bx, b2 = b1+v-ky
c            if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
c          lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c  wrk2  : real array of dimension (lwrk2). used as workspace, but
c          only in the case a rank deficient system is encountered.
c  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
c          the array wrk2 as declared in the calling (sub)program.
c          lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
c          where u,v and b2 are as above. if there are enough data
c          points, scattered uniformly over the approximation domain
c          and if the smoothing factor s is not too small, there is a
c          good chance that this extra workspace is not needed. a lot
c          of memory might therefore be saved by setting lwrk2=1.
c          (see also ier > 10)
c  iwrk  : integer array of dimension (kwrk). used as workspace.
c  kwrk  : integer. on entry kwrk must specify the actual dimension of
c          the array iwrk as declared in the calling (sub)program.
c          kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
c  ier   : integer. unless the routine detects an error, ier contains a
c          non-positive value on exit, i.e.
c   ier=0  : normal return. the spline returned has a residual sum of
c            squares fp such that abs(fp-s)/s <= tol with tol a relat-
c            ive tolerance set to 0.001 by the program.
c   ier=-1 : normal return. the spline returned is an interpolating
c            spline (fp=0).
c   ier=-2 : normal return. the spline returned is the weighted least-
c            squares polynomial of degrees kx and ky. in this extreme
c            case fp gives the upper bound for the smoothing factor s.
c   ier<-2 : warning. the coefficients of the spline returned have been
c            computed as the minimal norm least-squares solution of a
c            (numerically) rank deficient system. (-ier) gives the rank.
c            especially if the rank deficiency which can be computed as
c            (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
c            curate. they could also seriously depend on the value of
c            eps.
c   ier=1  : error. the required storage space exceeds the available
c            storage space, as specified by the parameters nxest and
c            nyest.
c            probably causes : nxest or nyest too small. if these param-
c            eters are already large, it may also indicate that s is
c            too small
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=2  : error. a theoretically impossible result was found during
c            the iteration proces for finding a smoothing spline with
c            fp = s. probably causes : s too small or badly chosen eps.
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=3  : error. the maximal number of iterations maxit (set to 20
c            by the program) allowed for finding a smoothing spline
c            with fp=s has been reached. probably causes : s too small
c            there is an approximation returned but the corresponding
c            weighted sum of squared residuals does not satisfy the
c            condition abs(fp-s)/s < tol.
c   ier=4  : error. no more knots can be added because the number of
c            b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
c            the number of data points m.
c            probably causes : either s or m too small.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=5  : error. no more knots can be added because the additional
c            knot would (quasi) coincide with an old one.
c            probably causes : s too small or too large a weight to an
c            inaccurate data point.
c            the approximation returned is the weighted least-squares
c            spline according to the current set of knots.
c            the parameter fp gives the corresponding weighted sum of
c            squared residuals (fp>s).
c   ier=10 : error. on entry, the input data are controlled on validity
c            the following restrictions must be satisfied.
c            -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
c            nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
c            xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
c            lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
c            kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
c            if iopt=-1: 2*kx+2<=nx<=nxest
c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
c                        2*ky+2<=ny<=nyest
c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
c            if iopt>=0: s>=0
c            if one of these conditions is found to be violated,control
c            is immediately repassed to the calling program. in that
c            case there is no approximation returned.
c   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
c            space for computing the minimal least-squares solution of
c            a rank deficient system of linear equations. ier gives the
c            requested value for lwrk2. there is no approximation re-
c            turned but, having saved the information contained in nx,
c            ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
c            the dimension of the array wrk2 accordingly, the user can
c            continue at the point the program was left, by calling
c            surfit with iopt=1.
c
c further comments:
c  by means of the parameter s, the user can control the tradeoff
c   between closeness of fit and smoothness of fit of the approximation.
c   if s is too large, the spline will be too smooth and signal will be
c   lost ; if s is too small the spline will pick up too much noise. in
c   the extreme cases the program will return an interpolating spline if
c   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
c   very large. between these extremes, a properly chosen s will result
c   in a good compromise between closeness of fit and smoothness of fit.
c   to decide whether an approximation, corresponding to a certain s is
c   satisfactory the user is highly recommended to inspect the fits
c   graphically.
c   recommended values for s depend on the weights w(i). if these are
c   taken as 1/d(i) with d(i) an estimate of the standard deviation of
c   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
c   sqrt(2*m)). if nothing is known about the statistical error in z(i)
c   each w(i) can be set equal to one and s determined by trial and
c   error, taking account of the comments above. the best is then to
c   start with a very large value of s ( to determine the least-squares
c   polynomial and the corresponding upper bound fp0 for s) and then to
c   progressively decrease the value of s ( say by a factor 10 in the
c   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
c   approximation shows more detail) to obtain closer fits.
c   to choose s very small is strongly discouraged. this considerably
c   increases computation time and memory requirements. it may also
c   cause rank-deficiency (ier<-2) and endager numerical stability.
c   to economize the search for a good s-value the program provides with
c   different modes of computation. at the first call of the routine, or
c   whenever he wants to restart with the initial set of knots the user
c   must set iopt=0.
c   if iopt=1 the program will continue with the set of knots found at
c   the last call of the routine. this will save a lot of computation
c   time if surfit is called repeatedly for different values of s.
c   the number of knots of the spline returned and their location will
c   depend on the value of s and on the complexity of the shape of the
c   function underlying the data. if the computation mode iopt=1
c   is used, the knots returned may also depend on the s-values at
c   previous calls (if these were smaller). therefore, if after a number
c   of trials with different s-values and iopt=1, the user can finally
c   accept a fit as satisfactory, it may be worthwhile for him to call
c   surfit once more with the selected value for s but now with iopt=0.
c   indeed, surfit may then return an approximation of the same quality
c   of fit but with fewer knots and therefore better if data reduction
c   is also an important objective for the user.
c   the number of knots may also depend on the upper bounds nxest and
c   nyest. indeed, if at a certain stage in surfit the number of knots
c   in one direction (say nx) has reached the value of its upper bound
c   (nxest), then from that moment on all subsequent knots are added
c   in the other (y) direction. this may indicate that the value of
c   nxest is too small. on the other hand, it gives the user the option
c   of limiting the number of knots the routine locates in any direction
c   for example, by setting nxest=2*kx+2 (the lowest allowable value for
c   nxest), the user can indicate that he wants an approximation which
c   is a simple polynomial of degree kx in the variable x.
c
c  other subroutines required:
c    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
c
c  references:
c   dierckx p. : an algorithm for surface fitting with spline functions
c                ima j. numer. anal. 1 (1981) 267-283.
c   dierckx p. : an algorithm for surface fitting with spline functions
c                report tw50, dept. computer science,k.u.leuven, 1980.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author:
c    p.dierckx
c    dept. computer science, k.u. leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  creation date : may 1979
c  latest update : march 1987
c
c  ..
c  ..scalar arguments..
      real xb,xe,yb,ye,s,eps,fp
      integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
c  ..array arguments..
      real x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
c  ..local scalars..
      real tol
      integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,
     * lby,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek,
     * nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
c  ..function references..
      integer max0
c  ..subroutine references..
c    fpsurf
c  ..
c  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
c  before starting computations a data check is made. if the input data
c  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 70
      if(kx.le.0 .or. kx.gt.5) go to 70
      kx1 = kx+1
      if(ky.le.0 .or. ky.gt.5) go to 70
      ky1 = ky+1
      kmax = max0(kx,ky)
      km1 = kmax+1
      km2 = km1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 70
      if(m.lt.(kx1*ky1)) go to 70
      nminx = 2*kx1
      if(nxest.lt.nminx .or. nxest.gt.nmax) go to 70
      nminy = 2*ky1
      if(nyest.lt.nminy .or. nyest.gt.nmax) go to 70
      nest = max0(nxest,nyest)
      nxk = nxest-kx1
      nyk = nyest-ky1
      ncest = nxk*nyk
      nmx = nxest-nminx+1
      nmy = nyest-nminy+1
      nrint = nmx+nmy
      nreg = nmx*nmy
      ib1 = kx*nyk+ky1
      jb1 = ky*nxk+kx1
      ib3 = kx1*nyk+1
      if(ib1.le.jb1) go to 10
      ib1 = jb1
      ib3 = ky1*nxk+1
  10  lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 70
      if(xb.ge.xe .or. yb.ge.ye) go to 70
      do 20 i=1,m
        if(w(i).le.0.) go to 70
        if(x(i).lt.xb .or. x(i).gt.xe) go to 70
        if(y(i).lt.yb .or. y(i).gt.ye) go to 70
  20  continue
      if(iopt.ge.0) go to 50
      if(nx.lt.nminx .or. nx.gt.nxest) go to 70
      nxk = nx-kx1
      tx(kx1) = xb
      tx(nxk+1) = xe
      do 30 i=kx1,nxk
        if(tx(i+1).le.tx(i)) go to 70
  30  continue
      if(ny.lt.nminy .or. ny.gt.nyest) go to 70
      nyk = ny-ky1
      ty(ky1) = yb
      ty(nyk+1) = ye
      do 40 i=ky1,nyk
        if(ty(i+1).le.ty(i)) go to 70
  40  continue
      go to 60
  50  if(s.lt.0.) go to 70
  60  ier = 0
c  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncest*ib3
      lf = la+ncest*ib1
      lff = lf+ncest
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbx = lh+ib3
      nek = nest*km2
      lby = lbx+nek
      lsx = lby+nek
      lsy = lsx+m*km1
      call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     * eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx,
     * ny,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),
     * wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy),
     * wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  70  return
      end
