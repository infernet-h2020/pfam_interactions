ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  Subroutine to calculate GL-score and secondary structure identity
cccc  for a given fragment length
cccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_nGL (ist1,istp1,ist2,istp2,nlen,GL,ssa,ssan)  
      parameter (maxres=3000)           ! no. of residues  
      dimension invmap(maxres)
      dimension xo1(maxres),yo1(maxres),zo1(maxres)
      dimension xo2(maxres),yo2(maxres),zo2(maxres)
      dimension dis2(maxres)
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/

       da=0.5
       d02=(da)**2
       GL=0.0
       ssa=0.0

       idiff1=(istp1-ist1)+1
       idiff2=(istp2-ist2)+1
       istp2_mod=istp2
    
       if (idiff1.lt.idiff2) then   ! if first fragment is less than 
                                    ! second one make second one equal to first fragment
        istp2_mod=ist2+idiff1-1
       endif
        
       i=ist1
       n_al=0
       iss=0
       issn=0
       do j=ist2,istp2_mod
         n_al=n_al+1
         r_1(1,n_al)=xa(1,i,0)
         r_1(2,n_al)=xa(2,i,0)
         r_1(3,n_al)=xa(3,i,0)
         r_2(1,n_al)=xa(1,j,1)
         r_2(2,n_al)=xa(2,j,1)
         r_2(3,n_al)=xa(3,j,1)
         xo1(n_al)=xa(1,i,0)
         yo1(n_al)=xa(2,i,0)
         zo1(n_al)=xa(3,i,0)
         xo2(n_al)=xa(1,j,1)
         yo2(n_al)=xa(2,j,1)
         zo2(n_al)=xa(3,j,1)
         if (isec(i).eq.jsec(j).and.isec(i).gt.1)iss=iss+1
         if (isec(i).eq.jsec(j))issn=issn+1
         i=i+1
       enddo

       ssa=float(iss)/float(n_al)
       ssan=float(issn)/float(n_al)

       call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier)  !u rotate r_1 to r_2
       GL1=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
          yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
          zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL1=GL1+(1/(1+dis2(i)/(d02)))
       enddo
       GL1=GL1/float(n_al)

       dxs=4.5
       d002=(dxs**2)

   20  continue
       j=0
       do i=1,n_al
         if (dis2(i).le.d002) then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
        endif
       enddo
       if (j.le.3) then
        d002=d002+0.5
       goto 20
      endif
       L=j

       call u3b (w,r_1,r_2,L,1,rms,u,t,ier)  !u rotate r_1 to r_2

       GL2=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
          yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
          zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL2=GL2+(1/(1+dis2(i)/(d02)))
       enddo
       GL2=GL2/float(n_al)

   21  continue
       j=0
       do i=1,n_al
         if (dis2(i).le.d002) then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
        endif
       enddo
       if (j.le.3) then
        d002=d002+0.5
       goto 21
      endif
       L=j

       call u3b (w,r_1,r_2,L,1,rms,u,t,ier)  !u rotate r_1 to r_2

       GL3=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
          yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
          zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL3=GL3+(1/(1+dis2(i)/(d02)))
       enddo
       GL3=GL3/float(n_al)
       
        GL=GL1
        if (GL2.gt.GL) GL=GL2
        if (GL3.gt.GL) GL=GL3
        
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
cccc Subroutine for aligning the fragments based on GL-score matrix or
cccc SSA+GL-score matrix.
cccc
cccc A variation Waterman-Eggert algorithm is applied to generate
cccc suboptimal alignments after getting the best optimal alignments.
cccc The step of generating suboptimal alignment:
cccc 1. Recompute the F matrix with not allowing match for those position,
cccc    which are aligned in the previous alignment.
cccc 2. Backtrace the F matrix generating the alignment
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine fragdp (nfr1,nfr2,smtx,gap,ifrali,ntot)
       parameter (maxfr=1200)
       dimension smtx(maxfr,maxfr),ifrali(100,maxfr)
       dimension idir(0:maxfr,0:maxfr),lcmap(maxfr)
       dimension val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr)

       ! initialize the F matrix
         val(0,0)=0.0
         do i=1,nfr1
          idir(i,0)=0
          val(i,0)=0.0
         enddo

         do j=1,nfr2
          idir(0,j)=0
          val(0,j)=0.0
          lcmap(j)=-1
         enddo
        
       ! compute F matrix
            do j=1,nfr2
                do i=1,nfr1
                 d=val(i-1,j-1)+smtx(i,j)
                 h=val(i-1,j)+gap
                 v=val(i,j-1)+gap 
                 lf(i,j)=0
             
                 if (d.ge.h.and.d.ge.v) then
                     val(i,j)=d
                     idir(i,j)=1
                 else
                     if (h.ge.v) then
                         val(i,j)=h
                         idir(i,j)=2
                     else
                         val(i,j)=v
                         idir(i,j)=3
                     endif
                 endif
                enddo
            enddo
        
        amax=0.0
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,1,amax)
        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo
        
        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,2,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,3,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,4,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Routine to filter the redundant fragment alignments
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine filter (ifrali,ntot,nfr2,nret,jali) 
        parameter (maxfr=1200)
        dimension ifrali (100,maxfr),jali(100,maxfr),istate(100)

        do ii=1,ntot
            istate(ii)=1 ! all are unique alignment
        enddo
        
         do ii=1,ntot-1
            do jj=ii+1,ntot
             if (istate(jj).eq.1) then
               nlc=0
               do kk=1,nfr2
                if (ifrali(ii,kk).eq.ifrali(jj,kk)) nlc=nlc+1
               enddo
               if (nlc.eq.nfr2) istate(jj)=0
             endif
            enddo
         enddo

         nret=0
         do ii=1,ntot
            if (istate(ii).eq.1) then
                nret=nret+1
                do jj=1,nfr2
                    jali(nret,jj)=ifrali(ii,jj)
                enddo
            endif
         enddo
        
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Subroutine for calculating the F matrix and returning the traceback matrix
cccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
       parameter (maxfr=1200)
       dimension smtx(maxfr,maxfr)
       dimension idir(0:maxfr,0:maxfr)
       dimension val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr)

       ! initialize F matrix
         val(0,0)=0.0
         do i=1,nfr1
            idir(i,0)=0
            val(i,0)=0.0
         enddo

         do j=1,nfr2
            idir(0,j)=0
            val(0,j)=0.0
         enddo

      ! fill F matrix using a modified recursion relation

        do j=1,nfr2
             do i=1,nfr1
              if (lf(i,j).eq.0) then
                 d=val(i-1,j-1)+smtx(i,j)
              else
                 d=0.0
              endif
             
              h=val(i-1,j)+gap
              v=val(i,j-1)+gap 
             
               if (d.ge.h.and.d.ge.v) then
                 val(i,j)=d
                 idir(i,j)=1
               else
                 if (h.ge.v) then
                     val(i,j)=h
                     idir(i,j)=2
                 else
                     val(i,j)=v
                     idir(i,j)=3
                 endif
               endif
             enddo
        enddo

       return
       end

cccc Subroutine to traceback the filled F matrix 
cccc 
       subroutine traceback (nfr1,nfr2,val,idir,lf,lcmap,iround,amax)
       parameter (maxfr=1200)
       dimension idir(0:maxfr,0:maxfr),val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr),lcmap(maxfr)
       
        ! Traceback GLOBAL optimal/suboptimal alignment with no penalty for end gaps
        ! find the highest score on the last row and column and start
        ! alignment from that position

        i=nfr1
        j=nfr2
        amax2=0.0
        if (iround.eq.1) then
            amax=40000.00
        endif

        do ii=1,nfr1
            if (val(ii,nfr2).gt.amax2.and.val(ii,nfr2).lt.amax) then
                amax2=val(ii,nfr2)
                i=ii
            endif
        enddo

        do jj=1,nfr2
            if (val(nfr1,jj).gt.amax2.and.val(nfr1,jj).lt.amax) then
                amax2=val(nfr1,jj)
                j=jj
                i=nfr1
            endif
        enddo

        if (iround.eq.1) then
            amax=amax2
        endif
    
        do while ( (i.gt.0) .and. (j.gt.0) )
            if ( idir(i,j) .eq. 1 ) then
                lcmap(j)=i
                lf(i,j)=1
                i=i-1
                j=j-1
            elseif (idir(i,j) .eq. 2) then
                i=i-1
            elseif (idir(i,j) .eq. 3) then
                j=j-1
            endif
        enddo

       return
       end
