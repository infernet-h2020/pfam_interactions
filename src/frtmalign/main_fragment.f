************************************************************************
**
** The protein is divided into set of non-overlapping fragments of a
** given length and we would comapre one fragment against all others
** using cheap version of TM-score (GL-score). This score is finally 
** used for doing a repeated match DP to get many fragment alignment.
** The aligned sets of fragment is joined to make the initial seed 
** alignment. The initial alignment is expanded with a routine of
** extension using DP and TM-score. Now, we are using a new DP routine.
**
************************************************************************

       subroutine fragscan (dx,dxs,mapr,irt,TMmax)
       parameter (maxres=3000)
       parameter (maxfr=1200)   ! max. no. of fragments
       parameter (mfrali=100)   ! max. no. of aligned fragments

       dimension invmap(maxres),mapA(maxres),mapr(maxres)
       dimension dpst(maxfr,maxfr),sspt(maxfr,maxfr),iali(100,maxfr)
       dimension smtx(maxfr,maxfr),ifrali(100,maxfr),gap(4),sspta(maxfr,maxfr)
       dimension ilen(4),glm(2),gls(maxfr,maxfr,2),ssl(maxfr,maxfr,2)
       
       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)
        
       fr=float(min(nseq1,nseq2))/float(max(nseq1,nseq2))
       minseq=min(nseq1,nseq2)
       
       if (minseq.lt.100) then
        ipp=8
       else
        ipp=12
       endif

       lenfr=ipp
       nfr1=nseq1/lenfr ! no. of fragments for seq1
       nfr2=nseq2/lenfr ! no. of fragments for seq2
       glmax=0.0
       ntot=0
       slmax=0.0
       slmaxn=0.0
    
       do ii=1,nfr2
        istart2=(ii-1)*lenfr+1
        iend2=istart2+lenfr-1
        if (ii.eq.nfr2) iend2=nseq2

            do jj=1,nfr1
             istart1=(jj-1)*lenfr+1
             iend1=istart1+lenfr-1
             if (jj.eq.nfr1) iend1=nseq1

             call get_nGL(istart1,iend1,istart2,iend2,lenfr,GL,ssa,ssna)
              
              if (GL.gt.glmax) then
                glmax=GL
              endif
               if (ii.ne.nfr2.and.jj.ne.nfr1) then
                if (ssa.gt.slmax) slmax=ssa
                if (ssna.gt.slmaxn) slmaxn=ssna
               endif

              dpst(jj,ii)=GL
              sspt(jj,ii)=ssa
              sspta(jj,ii)=ssna
              ntot=ntot+1
            enddo
       enddo

       gap(1)=-0.6
       gap(2)=-0.1
       gap(3)=0.0
        
       tmmax=0.0
       ntot=0 
        do ii=1,nfr1
            do jj=1,nfr2
                smtx(ii,jj)=dpst(ii,jj)/glmax
            enddo
        enddo

       do ipp=1,3
        gap_open=gap(ipp)
        call fragdp (nfr1,nfr2,smtx,gap_open,ifrali,ntot) ! returns ntot 
       enddo

       if (slmax .lt. 0.3) then
        do ii=1,nfr1
         do jj=1,nfr2
            sspt(ii,jj)=sspta(ii,jj)
         enddo
        enddo
        slmax=slmaxn
       endif
       
        do ii=1,nfr1
            do jj=1,nfr2
                smtx(ii,jj)=0.5*(sspt(ii,jj)/slmax)+0.5*(dpst(ii,jj)/glmax)
            enddo
        enddo

       do ipp=1,3
        gap_open=gap(ipp)
        call fragdp (nfr1,nfr2,smtx,gap_open,ifrali,ntot) ! returns ntot 
       enddo

       call filter (ifrali,ntot,nfr2,nali,iali)
       call calbesttm (dx,dxs,nfr1,nfr2,lenfr,nali,iali,irt,mapr,tmmax)

       end

************************************************************************
*****          Get the alignment from the fragments and do heuristic
*****          iteration
************************************************************************

       subroutine calbesttm (dx,dxs,nfr1,nfr2,nlen,ntot,ifrali,irt,mapr,tmmax) 
       parameter (maxres=3000)
       parameter (maxfr=1200)
       dimension ifrali(100,maxfr),mapr(maxres),invmap(maxres)
       dimension lcmap(maxfr),invmapr(maxres)

       common /length/ nseq1,nseq2

       tmmax=0.0

        do iali=1,ntot
            do ii=1,nfr2
                lcmap(ii)=ifrali(iali,ii)
            enddo

            call fillinvmap (invmap)

            do ii=1,nfr2
                ifr1=lcmap(ii)
                ifr2=ii
               if (ifr1.gt.0) then
                istart1=(ifr1-1)*nlen+1
                istart2=(ifr2-1)*nlen+1
                iend1=istart1+nlen-1
                iend2=istart2+nlen-1
                if (ifr1.eq.nfr1) iend1=nseq1
                if (ifr2.eq.nfr2) iend2=nseq2
                if ((iend1-istart1).lt.(iend2-istart2)) then
                    iend2=(iend1-istart1)+istart2
                endif
                 itmp=istart1
                 do jj=istart2,iend2
                    invmap(jj)=itmp
                    itmp=itmp+1
                 enddo
               endif
            enddo
            
            tminp=tmmax
            
            if (irt .eq. 0) then
                call caltmsc (dx,dxs,invmap,tminp,invmapr,b1tmmax)
            elseif (irt.eq.1) then
                call caltmsc_fast (dx,dxs,invmap,tminp,invmapr,b1tmmax)
            endif

            if (b1tmmax.gt.tmmax) then
                do ii=1,nseq2
                    mapr(ii)=invmapr(ii)
                enddo
                tmmax=b1tmmax
            endif

        enddo ! for total no. of alignments

       return
       end

*********************************************************
**      Routine for slow calculation of best alignment
**      Scans many alignments!!
*********************************************************

       subroutine caltmsc (dx,dxs,map,tmch,mapr,tmscm)
       parameter (maxres=3000)
       dimension map(maxres),mapr(maxres),score(maxres,maxres)
       dimension scorea(maxres,maxres),gap(2),invmap_i(maxres)
       dimension invmapr(maxres),scoret(maxres,maxres)
       dimension ssp(maxres,maxres),ssp1(maxres,maxres),mapA(maxres,maxres)
    
       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)

       gap(1)=-1.0
       gap(2)=0.0
       
       tmscm=0.0
       call get_initial3 (dx,map,scoret,scorea)
       tminp_lc=tmch
       istart=1
       iend=2

       if (tmch .gt. 0.50) then
            iend=1
       endif

       do ii=istart,iend
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tminp_lc) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         endif
       enddo

       if (tmscm .gt. tminp_lc) tminp_lc=tmscm

       if (tminp_lc .gt. 0.50) then
            istart=2
       endif

       do ii=istart,iend
         gap_open=gap(ii)
         call dp (scoret,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo

       if (tmscm .gt. tminp_lc) tminp_lc=tmscm

       if (tminp_lc .lt. 0.50) then
       
       call get_score (dx,dxs,map,scorea,TM)
       do ii=2,2
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,6,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo
       
       endif

       return
       end

*********************************************************
**** General routine to fill the starting alignment !!
*********************************************************
       subroutine fillinvmap (invmapr)
       parameter (maxres=3000)
       dimension invmapr(maxres)
       common /length/ nseq1,nseq2

        do jj=1,nseq2
            invmapr(jj)=-1
        enddo

       return
       end

*********************************************************
**      Routine for fast calculation of best alignment
**      Scans less number of alignments!!
*********************************************************

       subroutine caltmsc_fast (dx,dxs,map,tmch,mapr,tmscm)
       parameter (maxres=3000)
       dimension map(maxres),mapr(maxres),score(maxres,maxres)
       dimension scorea(maxres,maxres),gap(2),invmap_i(maxres)
       dimension invmapr(maxres),scoret(maxres,maxres)
       dimension ssp(maxres,maxres),ssp1(maxres,maxres),mapA(maxres,maxres)

       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)

       gap(1)=-1.0
       gap(2)=0.0

       tmscm=0.0
       call get_initial3 (dx,map,scoret,scorea)
       tminp_lc=tmch

       do ii=1,1
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tminp_lc) then
            call make_iter(dx,dxs,score,14,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         endif
       enddo

       if (tminp_lc .lt. tmscm ) tminp_lc=tmscm
                    
       do ii=2,2
         gap_open=gap(ii)
         call dp (scoret,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,14,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo

       if (tminp_lc .lt. tmscm ) tminp_lc=tmscm

       if (tminp_lc .lt. 0.50) then
        call get_score (dx,dxs,map,scorea,TM)
            do ii=2,2
                gap_open=gap(ii)
                call dp (scorea,gap_open,invmap_i,dpout)
                call get_score (dx,dxs,invmap_i,score,TM)
                tTM=TM
                if (tTM .gt. tmscm) then
                    call make_iter(dx,dxs,score,14,invmap_i,tTM,invmapr)
                    call getbest (tTM,invmapr,tmscm,mapr)
                else
                    call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
                    call getbest (tTM,invmapr,tmscm,mapr)
                endif
            enddo
       endif

       return
       end
 
******************************************************************
