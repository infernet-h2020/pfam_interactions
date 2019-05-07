cccc Subroutines for TM-align
cccc TMsearch, cal_tmscore

      subroutine TMsearch(dx,dxs,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm,isearch)

      parameter (maxres=3000)
      dimension x1(maxres),y1(maxres),z1(maxres),n1(maxres)
      dimension x2(maxres),y2(maxres),z2(maxres),n2(maxres)
      dimension xa(maxres),ya(maxres),za(maxres)
      dimension k_ali(maxres),k_ali0(maxres)
      dimension L_ini(100),iq(maxres)
      double precision score,score_max
      dimension iL0(maxres),i_ali(maxres)

      common /stru/xt(maxres),yt(maxres),zt(maxres),xb(maxres),yb(maxres),zb(maxres) !for cal_tmscore routine
      common /align/ iA(maxres),iB(maxres)  !shared in next subroutine
      common /length/ nseq1,nseq2
      logical nstat

ccc   RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/

cccc  convert input data. Special case here of L1.eq.L2
      nstat=.false. ! used in cal_tmscore routine for decision 
                    ! for including residues for TM-score calculation
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

cccc d0 parameters 
      d0_min=0.5
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
       d0_search=dxs

cccc  iterative parameters ----->
      n_it=20                   !maximum number of iterations
      n_init_max=6              !maximum number of L_init
      n_init=0                  ! number of fragments L/2**(n)
      L_ini_min=4

      if(n_ali.lt.4)L_ini_min=n_ali

      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

cccc  find the maximum score starting from local structures superposition

      score_max=-1              !TM-score
      if (isearch.eq.1) then
        ijump=40
        if (min(nseq1,nseq2).le.80)ijump=30
      endif
      if (isearch.ge.2) ijump=1
      if (isearch.eq.1.or.isearch.eq.2) nstat=.true.

      do 333 i_init=1,n_init
        L_init=L_ini(i_init)    ! Length to be used for scanning L
        iL_max=n_ali-L_init+1   ! maximum number of scans (Nali-L+1), shifting by one residue
c      if (isearch.gt.2) write(*,*)iL_max,score_max,d0_search,d0

        if (isearch.eq.1) then
         k=0
            do i=1,iL_max,ijump     ! this is the for the real quick result
                k=k+1               ! when isearch.eq.1
                iL0(k)=i            ! storing the starting residue for scan
            enddo

            if(iL0(k).lt.iL_max)then
                k=k+1
                iL0(k)=iL_max   ! fixing the last run, quite arbitarary!!
            endif
            n_shift=k           ! number of shifts
        else
           n_shift=iL_max
        endif    

              do 300 i_shift=1,n_shift
                    if (isearch.eq.1) then
                     iL=iL0(i_shift)
                    else
                     iL=i_shift
                    endif

                    LL=0
                    ka=0
                      do i=1,L_init
                          k=iL+i-1          ![1,n_ali] common aligned
                          r_1(1,i)=xa(iA(k))
                          r_1(2,i)=ya(iA(k))
                          r_1(3,i)=za(iA(k))
                          r_2(1,i)=xb(iB(k))
                          r_2(2,i)=yb(iB(k))
                          r_2(3,i)=zb(iB(k))
                          LL=LL+1
                          ka=ka+1
                          k_ali(ka)=k
                      enddo

                    call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

                    if(i_init.eq.1)then  !global superposition
                     armsd=dsqrt(rms/LL)
                     Rcomm=armsd
                    endif

                    do j=1,nseqA
                      xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                      yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                      zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
                    enddo

                    d=d0_search-1                                 ! reduce the search rms first
                    call cal_tmscore(d,d0,n_ali,nstat,score,i_ali,n_cut) !init, get scores, n_cut+i_ali(i) for iteration
                    if(score_max.lt.score)then
                        score_max=score
                        ka0=ka
                        do i=1,ka0
                            k_ali0(i)=k_ali(i)  ! alignment stored
                        enddo
                    endif
cccc  iteration for extending
                d=d0_search+1
                do 301 it=1,n_it
                    LL=0
                    ka=0
                    do i=1,n_cut
                        m=i_ali(i)          ![1,n_ali]
                        r_1(1,i)=xa(iA(m))
                        r_1(2,i)=ya(iA(m))
                        r_1(3,i)=za(iA(m))
                        r_2(1,i)=xb(iB(m))
                        r_2(2,i)=yb(iB(m))
                        r_2(3,i)=zb(iB(m))
                        ka=ka+1
                        k_ali(ka)=m
                        LL=LL+1
                    enddo
                    call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

                    do j=1,nseqA
                     xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                     yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                     zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
                    enddo

                    call cal_tmscore(d,d0,n_ali,nstat,score,i_ali,n_cut) !get scores, n_cut+i_ali(i) for iteration

                    if(score_max.lt.score)then
                        score_max=score
                        ka0=ka
                        do i=1,ka
                            k_ali0(i)=k_ali(i)
                        enddo
                    endif

                    if(it.eq.n_it)goto 302

                    if(n_cut.eq.ka)then
                        neq=0
                        do i=1,n_cut
                            if(i_ali(i).eq.k_ali(i))neq=neq+1
                        enddo
                        if(n_cut.eq.neq)goto 302
                    endif
 301            continue             !for iteration
 302           continue
 300          continue                !for shift
 333  continue                  !for initial length, L_ali/M

cccc  return the final rotation with rotated (xtm1 or x1)
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo

      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      
      TM=score_max

      return
      end

cccc rotuine to calculate TM-score given the superposed coordinates
cccc IMP: superposed coordinates are passed in common
cccc d is for selecting the residues with certain cutoff and d0 is for
cccc TM-score calculation. nstat is type of calculation d or d8

      subroutine cal_tmscore (d,d0,n_ali,nstat,score,i_ali,n_cut)
      parameter (maxres=3000)
      double precision score,score_max
      dimension i_ali(maxres)

      common/stru/xt(maxres),yt(maxres),zt(maxres),xb(maxres),yb(maxres),zb(maxres)
      common /align/ iA(maxres),iB(maxres)
      common /length/ nseq1,nseq2
      common /dinfo/ d8
      logical nstat

      d_cp=d        ! for preventing superposition of 2 residues and less
      nseqB=n_ali   ! normalized for TMscore

  116 continue
        n_cut=0                   !number of residue-pairs dis<d, for iteration
        score_sum=0               !TMscore
        do k=1,n_ali
            i=iA(k)                ![1,nseqA] reoder number of structureA
            j=iB(k)                ![1,nseqB]
            dis=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+(zt(i)-zb(j))**2)
                if(dis.lt.d_cp)then
                    n_cut=n_cut+1
                    i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
                endif

                if(nstat) then
                    if(dis.le.d8)then
                        score_sum=score_sum+1.0/(1+(dis/d0)**2)
                    endif
                else
                        score_sum=score_sum+1.0/(1+(dis/d0)**2)
                endif
        enddo

        if(n_cut.gt.0.and.n_cut.lt.3) then
            d_cp=d_cp+0.5
            goto 116
        endif
        score=score_sum !TM-score

      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(n), x(3, n), y(3, n), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end
