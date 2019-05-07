cccc Subroutines for Fr-TM-align 
cccc get_initial3,dp,make_iter,get_score1,get_score

*************************************************************
** 	DP routine (score matrix, gap open, alignment)
**  returns invmap_i
** Fix: In the Iteration routine if the number of residues aligned is
**      less than 5 residues the method aborts checking for any 
**      improvement in the alignment.
*************************************************************
      subroutine dp (score,gap_open,invmap_i,dpout)
      
      parameter (maxres=3000)             ! no. of residues  
      dimension score(maxres,maxres),invmap_i(maxres)
      logical*1 dir
      dimension dir(0:maxres,0:maxres),val(0:maxres,0:maxres)
      real h,v

      common /length/ nseq1,nseq2

cccc initialize matrix
      do ii=1,nseq2
        invmap_i(ii)=-1
      enddo

       val(0,0)=0.0
       do i=1,nseq1
         dir(i,0)=.false.
         val(i,0)=0.0
       enddo
       do j=1,nseq2
         dir(0,j)=.false.
         val(0,j)=0.0
         invmap_i(j)=-1
       enddo

cccc fill the matrix and path

       do j=1,nseq2
        do i=1,nseq1
           d=val(i-1,j-1)+score(i,j)
           h=val(i-1,j)
           v=val(i,j-1)

            if(dir(i-1,j))h=h+gap_open
            if(dir(i,j-1))v=v+gap_open
 
                if((d.ge.h).and.(d.ge.v)) then
                    dir(i,j)=.true.
                    val(i,j)=d
                else
                    dir(i,j)=.false.
                        if(v.ge.h)then
                            val(i,j)=v
                        else
                            val(i,j)=h
                        end if
                endif
        enddo
       enddo
 
cccc   extract the alignment:
       i=nseq1
       j=nseq2

       amax=0.0
       icho=0
       minseq=nseq1
       if (minseq .lt. nseq2) then
         icho=1
       endif
      
       if (icho .eq. 1) then
        do jj=nseq2,1,-1
           if (val(nseq1,jj).gt.amax) then
             amax=val(nseq1,jj)
             j=jj
           endif
      enddo
       else
        do ii=nseq1,1,-1
          if (val(ii,nseq2).gt.amax) then
             amax=val(ii,nseq2)
             i=ii
          endif
        enddo
       endif

       do while((i.gt.0).and.(j.gt.0))
         if(dir(i,j))then
           invmap_i(j)=i
           i=i-1
           j=j-1
         else
           h=val(i-1,j)
           v=val(i,j-1)

           if(dir(i-1,j))h=h+gap_open
           if(dir(i,j-1))v=v+gap_open

           if(v.ge.h) then
             j=j-1
           else
             i=i-1
           endif
         endif
       enddo
 
       j=nseq2
       do while(invmap_i(j).lt.1)
          j=j-1
       enddo
       dpout=val(invmap_i(j),j)
 
      return
      end 


**************************************************************
**  a kind of mixing of two previous alignments.
**  Takes the alignment and retuns the alignment!!
**************************************************************

      subroutine get_initial3 (d0,invmapi,score,scorer)
      parameter (maxres=3000)             ! no. of residues  
      dimension invmapr(maxres),invmap_i(maxres),invmapi(maxres)
      dimension score(maxres,maxres),scorer(maxres,maxres)

      common /coord/ xa(3,maxres,0:1)
      common /secstr/ isec(maxres),jsec(maxres)
      common /length/ nseq1,nseq2

cccc  score matrix 
      do j=1,nseq2
         invmap_i(j)=invmapi(j)
      enddo

      call get_score1 (d0,invmap_i,score)          !get score(i,j) using RMSD martix

      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               scorer(i,j)=0.5+score(i,j)
            else
               scorer(i,j)=score(i,j)
            endif
         enddo
      enddo

      return
      end

******************************************************************
**  with invmap_i(i) calculate score(i,j) using RMSD rotation
******************************************************************
      subroutine get_score1 (dx,invmapi,score)
      parameter (maxres=3000)             ! no. of residues  
      dimension invmapi(maxres)
      dimension score(maxres,maxres)

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/

cccc  calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,nseq2
         i=invmapi(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0) ! for rotation matrix
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo

cccc  calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      
      d0_min=0.5
      d01=dx+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01**2

      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+(dd/d02))
         enddo
      enddo

      return
      end

******************************************************************
**  with invmap_i(i) calculate TM-score and martix score(i,j) for rotation
**  isearch takes 3 values for 3 different kind of TM-score calculation
**  isearch=1, searching with long jumps of 40 residues in TMsearch routine
******************************************************************

      subroutine get_score (d0,d0_s,invmap_i,score,TM)
      parameter (maxres=3000)             ! no. of residues  
      dimension invmap(maxres),invmap_i(maxres)
      dimension score(maxres,maxres)
      dimension xtm1(maxres),ytm1(maxres),ztm1(maxres)
      dimension xtm2(maxres),ytm2(maxres),ztm2(maxres)

      common /coord/ xa(3,maxres,0:1)
      common /n1n2/ n1(maxres),n2(maxres)
      common /length/ nseq1,nseq2

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/

cccc  calculate RMSD between aligned structures and rotate the structures -->
      anseq=min(nseq1,nseq2)    ! for normalization of TM-score
      d02=d0**2

      n_al=0
      do j=1,nseq2
         i=invmap_i(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1

            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)

ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
      d0_input=d0
      isearch=1     ! VERY IMP. determines way TM-score is calculate and reported !

      call TMsearch (d0_input,d0_s,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm,isearch) !simplified search engine 
        
       TM=(TM/anseq)

cccc output score matrix score(i,j)
cccc protein 1 or r_1 is need to rotated since in previous the rotation is for aligned subset!!
        do i=1,n_al
            r_2(1,i)=xtm1(i)
            r_2(2,i)=ytm1(i)
            r_2(3,i)=ztm1(i)
        enddo

        call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2

          do i=1,nseq1
            xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
            yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
            zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
            do j=1,nseq2
                dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
                score(i,j)=1/(1+dd/d02)
            enddo   
          enddo

       return
       end

******************************************************************
** Make iteration for a given alignment with two gap opening's
**
******************************************************************

       subroutine make_iter (d0,d0_s,score,niter,invmapi,aTM,invmapr)
       parameter (maxres=3000)
       dimension gap(10),score(maxres,maxres),invmapi(maxres)
       dimension invmap_i(maxres),invmap(maxres),invmapr(maxres)

       common /length/ nseq1,nseq2

       ngap=2
       gap(1)=-0.6
       gap(2)=0.0

cccc copy the alignment
       do j=1,nseq2
        invmap(j)=invmapi(j)
       enddo

        do 111 i_gap=1,ngap
        gap_open=gap(i_gap)

            do 222 id=1,niter
                call dp(score,gap_open,invmap_i,dpout)  ! returns invmap_i
                ilc_count=0
                do ik=1,nseq2
                    if (invmap_i(ik).gt.0)ilc_count=ilc_count+1
                enddo
                if (ilc_count.le.5) goto 33
                call get_score (d0,d0_s,invmap_i,score,TM)   ! refills the score matrix using invmap_i
                    if (TM.gt.aTM) then
                        aTM=TM
                        do j=1,nseq2
                            invmap(j)=invmap_i(j)
                        enddo
                    endif
                    if (id.gt.1) then
                        diff=abs(TM-TM_old)
                        if (diff.lt.0.000001)goto 33
                    endif
              TM_old=TM   
  222       continue
   33     continue
  111   continue

        do j=1,nseq2
          invmapr(j)=invmap(j)
        enddo

       return
       end
