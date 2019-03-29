************************************************************************
*                          Fr-TM-align (Version 1.0)                ****
*                                                                   ****
*   Program to align two protein structures, with maximizing the    ****
*   TM-score. The program has been extensively modified from the    ****
*   previous versions.                                              ****
*                                                                   ****        
*                                                                   ****
*   Email bugs to: spandit3@mail.gatech.edu                         ****
*                                                                   ****
*  Bug fixes:                                                       ****
*  A minor bug (typo) fixed on line 84 (Jan 10, 2009)               ****
*                                                                   ****
************************************************************************  

      program FrTMalign

      parameter (maxres=3000)           ! no. of residues       
      parameter (maxlen=2*maxres)       ! for alignment
      parameter (npr=600)               ! no. of proteins to align
      
      character*100 pdb(2),outfile,outname,buffer,fnam,pdbn(npr)

      character*3 aanam1(-2:20),resa(maxres,npr),resn(maxres,0:1)
      character*1 aanam2(-2:20),seqa(maxres,npr),seqn(maxres,0:1)
      character aseq1(maxlen),aseq2(maxlen),aseq3(maxlen)
      character*60 list,pdbid

      dimension invmap0(maxres),ires(maxres,0:1)
      dimension xtm1(maxres),ytm1(maxres),ztm1(maxres)
      dimension xtm2(maxres),ytm2(maxres),ztm2(maxres)
      dimension m1(maxres),m2(maxres)
      dimension xyz(3,maxres,npr),length(npr)
     
      data aanam1 /'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &              'PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN',
     &              'ARG','HIS','PHE','TYR','TRP','CYX','MSE'/

      data aanam2 /'X','G','A','S','C','V','T','I','P','M','D','N',
     &              'L','K','E','Q','R','H','F','Y','W','C','M'/

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /pdbinfo/ ires1(maxres,npr),resa,seqa
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure
      common /n1n2/ n1(maxres),n2(maxres)
      common /dinfo/ d8
ccc   RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/

cccc Arguments

       narg=iargc()

       if (narg .eq. 0) call instru
       call getarg (1,fnam)
       if (fnam.eq.'?'.or.fnam.eq.'-h'.or.fnam.eq.'--'.or.
     &     fnam.eq.'--help'.or.fnam.eq.'-help') then
           call instru
       endif

       mout=-1      ! output file
       mfix=-1      ! Fixed length TM-score
       mmat=0       ! Transformatrix output
       no=1         ! Number of pdbs
       irun_type=0  ! TYPE of Algorithm used for Fragment alignment
                    ! 0 is slow version
                    ! 1 is fast version

       i=1

       do while (i.le.narg) 
        call getarg(i,fnam)
        
         if (fnam(1:2).eq.'-o') then
            mout=1
            call getarg((i+1),outname)
            i=i+1
         elseif (fnam(1:2).eq.'-l') then
            mfix=1
            call getarg(i+1,fnam)
            read (fnam,*)L_fix
            i=i+1
         elseif (fnam(1:2).eq.'-m') then 
            call getarg(i+1,fnam)
            read (fnam,'(i1)')mmat
            i=i+1
         elseif (fnam(1:2).eq.'-a') then
            call getarg(i+1,fnam)
            read (fnam,'(i1)')ilist
            if (ilist.eq.1) then
                call getarg(i+2,fnam)
                read(fnam,'(a)')list
             i=i+1
            endif
            i=i+1
          elseif (fnam(1:2).eq.'-f') then
            call getarg(i+1,fnam)
            read (fnam,'(i1)')irun_type
            i=i+1
         else
            call getarg(i,fnam)
            read(fnam,'(a)')pdb(no)
            no=no+1
         endif
        i=i+1

       enddo
        
      if (ilist.eq.1 ) then
         if (no.lt.2) call instru
       if (mmat.eq.1) then
       write(6,*)'Warning: Rotation matrix is not written when a list ',
     &   '         of protein is compared against one protein'
       write(6,*)''
       endif

       if (mout.eq.1) then
       write(6,*)'Warning: Superposed coordinates are not written when',
     &   '         a list of protein  is compared against one protein'
       write(6,*)''
       endif
      endif

       if (ilist.eq.0 .and. no.lt.3) call instru

******************************************************************
****    main program starts
******************************************************************
cc checking if a list of PDBs are to be aligned against one protein, if

       if (ilist.eq.1) then
         ii=60
         do while (list(ii:ii).eq.' ')
             ii=ii-1
         enddo

         open (unit=1,file=list(1:ii),status='unknown')
         nopdb=1
         do while(.true.)
          read(1,*,end=11)pdbid
          call readpdb(pdbid,nopdb,nlen,aanam1,aanam2,xyz)
          pdbn(nopdb)=pdbid
          length(nopdb)=nlen
          nopdb=nopdb+1
         enddo
         close (1)
   11    continue
         call readpdb (pdb(1),nopdb,nlen,aanam1,aanam2,xyz)
         length(nopdb)=nlen
         pdbn(nopdb)=pdb(1)
       else
           do i=1,2
            call readpdb(pdb(i),i,nlen,aanam1,aanam2,xyz)
            length(i)=nlen
           enddo
        nopdb=2
       endif

cc initialize the coordinates xa() for superposition

      do 1166 ipdb=1,nopdb-1     ! main cycle
    
       if (ilist.eq.1) then
          i=0
          do j=1,length(ipdb)
           do k=1,3
            xa(k,j,i)=xyz(k,j,ipdb)   
           enddo
            seqn(j,i)=seqa(j,ipdb)
            resn(j,i)=resa(j,ipdb)
            ires(j,i)=ires1(j,ipdb)
          enddo
          
          if (ipdb.eq.1) then ! read second protein only once !!
            i=1
            jj=nopdb
            do j=1,length(jj)
             do k=1,3
                xa(k,j,i)=xyz(k,j,jj)   
             enddo
                seqn(j,i)=seqa(j,jj)
                resn(j,i)=resa(j,jj)
                ires(j,i)=ires1(j,jj)
            enddo
            nseq2=length(jj)
            pdb(2)=pdbn(nopdb)
          endif
        nseq1=length(ipdb)
        pdb(1)=pdbn(ipdb)
       else

        do i=0,1
         ik=i+1
         do j=1,length(ik)
            do k=1,3
               xa(k,j,i)=xyz(k,j,ik)   
            enddo
               seqn(j,i)=seqa(j,ik)
               resn(j,i)=resa(j,ik)
               ires(j,i)=ires1(j,ik)
         enddo
        enddo
        nseq1=length(1)
        nseq2=length(2)
       endif
    
cc fixing d0 for search ONLY with d0 fixed for small protein
      d0_min=0.5
      aminlen=min(nseq1,nseq2)
      d8=1.5*aminlen**0.3+3.5      !remove pairs with dis>d8 during search 
      if(aminlen.gt.15) then
        d0=1.24*(aminlen-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
        d0=d0_min
      endif
      nseq=max(nseq1,nseq2)
       do i=1,nseq
        n1(i)=i
        n2(i)=i
       enddo
      d0_search=d0
        
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.3.0)d0_search=4.5

      call super_align (d0,d0_search,invmap0,irun_type)

cc resuperose to find residues d < d8
      if (d0_search.lt.4.5)d0_search=4.5
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
      isearch=2 
      call TMsearch (d0_input,d0_search,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm,isearch) !TM-score with dis<d8 only

cccc for Output TM-score set d0 
      d0_min=0.5                !for output
      anseq=nseq2               !length for defining final TMscore
      if(mfix.eq.1)anseq=L_fix  !input length
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.4.5)d0_search=4.5

cccc remove dis>d8 in normal TM-score calculation for final report
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(seqn(m1(i),0).eq.seqn(m2(i),1) )then
               n_eq=n_eq+1
            endif
         endif
      enddo
      seq_id=float(n_eq)/(n_al+0.00000001)

      n8_al=j
      d0_input=d0
      isearch=3
      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm,isearch) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8/anseq       !TM-score after cutoff
      

      if (ilist. eq . 1 ) then
      write(6,108)pdb(1)(1:5),pdb(2)(1:5),n8_al,rmsd8_al,TM8,seq_id
 108  format(2a7,'= Aligned length=',I4,', RMSD=',f6.2,
     &     ', TM-score=',f7.5,', ID=',f5.3)
      endif
 1166 continue




cccc General output
cccc Write rotation matrix

      if (mmat.eq.1 .and. ilist. eq. 0) then
       open (unit=8,file='trf.mat',status='unknown')
       L=0
        do i=1,n8_al
          k=m1(i)
          L=L+1
          r_1(1,L)=xa(1,k,0)
          r_1(2,L)=xa(2,k,0)
          r_1(3,L)=xa(3,k,0)
          r_2(1,L)=xtm1(i)
          r_2(2,L)=ytm1(i)
          r_2(3,L)=ztm1(i)
        enddo
        if(L.gt.3)then
           call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2

           write(8,*)'Rotation matrix to rotate chain 1 to chain 2'
           write(8,*)'i          t(i)         u(i,1)         u(i,2) ',
     &               '        u(i,3)'
           do i=1,3
              write(8,204)i,t(i),u(i,1),u(i,2),u(i,3)
           enddo
  204     format(I2,f18.10,f15.10,f15.10,f15.10)
        endif
      endif
      close (8)
      
cccc
cccc Write superposed coordinates in rasmol script

      if (mout.eq.1 .and. ilist .eq. 0) then
        open (unit=9,file=outname,status='unknown')

**  rasmol   script:
         write(9,900)'load inline'
         write(9,900)'select atomno<2000'
         write(9,900)'wireframe .45'
         write(9,900)'select none'
         write(9,900)'select atomno>2000'
         write(9,900)'wireframe .20'
         write(9,900)'color white'
 900     format(A)

         do i=1,n8_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &           (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
            if(dis2.le.5)then
               write(9,901)m1(i)
               write(9,900)'color red'
               write(9,901)2000+m2(i)
               write(9,900)'color red'
            endif
         enddo
 901     format('select atomno=',I4)
         write (9,900)'select all'
         write (9,900)'exit'
         write (9,102) pdb(1),nseq1
         write (9,103) pdb(2),nseq2
 102     format ('REMARK   Protein 1:',A10,' Length =',i4)
 103     format ('REMARK   Protein 2:',A10,' Length =',i4)
         write (9,104) n8_al,rmsd8_al,TM8
 104     format ('REMARK   Structural alignment summary: ',
     &     ' Aligned length=',i4,'; RMSD =',f6.2,'; TM-score=',f7.5)      

***   chain1:
         do i=1,n8_al
            write(9,1237)m1(i),resn(m1(i),0),ires(m1(i),0),
     &           xtm1(i),ytm1(i),ztm1(i)
         enddo
         write(9,1238)          !TER
         do i=2,n8_al
            write(9,1239)m1(i-1),m1(i) !connect atoms
         enddo
***   chain2:
         do i=1,n8_al
            write(9,1237)2000+m2(i),resn(m2(i),1),ires(m2(i),1),
     $           xtm2(i),ytm2(i),ztm2(i)
         enddo
         write(9,1238)
         do i=2,n8_al
            write(9,1239)2000+m2(i-1),2000+m2(i)
         enddo
      endif
      close(9)
 1237    format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238    format('TER')
 1239    format('CONECT',I5,I5)
   
cccc
cccc Write the structure based sequence alignment
      if ( ilist .eq . 0 ) then
        call header 
        write(6,105)pdb(1),nseq1
        write(6,106)pdb(2),nseq2,int(anseq)
        write(6,*)
        write(6,107)n8_al,rmsd8_al,TM8,seq_id

 105  format('Chain 1:',A10,'  Size=',I4)
 106  format('Chain 2:',A10,'  Size=',I4,
     &' (TM-score is normalized by ',I4,')')
 107  format('Aligned length=',I4,', RMSD=',f6.2,
     &     ', TM-score=',f7.5,', ID=',f5.3)
      write(6,*)

       ii=0
       i1_old=1
       i2_old=1
       do i=1,n8_al
         do j=i1_old,m1(i)-1
            ii=ii+1
            aseq1(ii)=seqn(j,0)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do j=i2_old,m2(i)-1
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seqn(j,1)
            aseq3(ii)=' '
         enddo

         ii=ii+1
         aseq1(ii)=seqn(m1(i),0)
         aseq2(ii)=seqn(m2(i),1)
         dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &     (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
         if(dis2.le.5)then
           aseq3(ii)=':'
         else
           aseq3(ii)='.'
         endif
         i1_old=m1(i)+1
         i2_old=m2(i)+1
       enddo

      do i=i1_old,nseq1
         ii=ii+1
         aseq1(ii)=seqn(i,0)
         aseq2(ii)='-'
         aseq3(ii)=' '
      enddo
      do i=i2_old,nseq2
         ii=ii+1
         aseq1(ii)='-'
         aseq2(ii)=seqn(i,1)
         aseq3(ii)=' '
      enddo
      write(6,50)
 50   format('(":" denotes the residue pairs of distance < 5.0 ',
     &     'Angstrom)')
      write(6,10)(aseq1(i),i=1,ii)
      write(6,10)(aseq3(i),i=1,ii)
      write(6,10)(aseq2(i),i=1,ii)
 10   format(10000A1)
      write(6,*)
      endif

      stop
      end









******************************************************************
****            Subroutines
******************************************************************
cccc making initial superposition

      subroutine super_align(dx,dxs,invmap0,ir_type)
      parameter (maxres=3000)
      dimension score(maxres,maxres),invmap(maxres),invmap0(maxres)
      dimension invmap_i(maxres)

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure
       
       call assignssp   ! secondary structure assignment

       call fragscan(dx,dxs,invmap0,ir_type,atm)

      return
      end

cccc for getting the best alignment
      subroutine getbest(aTM,invmapi,aTMmax,invmapr)
      parameter (maxres=3000)
      dimension invmapi(maxres),invmapr(maxres)
      common /length/ nseq1,nseq2

        if (aTM.gt.aTMmax) then
            aTMmax=aTM
                do j=1,nseq2
                    invmapr(j)=invmapi(j)
                enddo
        endif
      return
      end

cccc making secondary structure assignment
      subroutine assignssp
      parameter (maxres=3000)
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

       do i=1,nseq1
          isec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq1)then
             dis13=diszy(0,j1,j3)
             dis14=diszy(0,j1,j4)
             dis15=diszy(0,j1,j5)
             dis24=diszy(0,j2,j4)
             dis25=diszy(0,j2,j5)
             dis35=diszy(0,j3,j5)
             isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       do i=1,nseq2
          jsec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq2)then
             dis13=diszy(1,j1,j3)
             dis14=diszy(1,j1,j4)
             dis15=diszy(1,j1,j5)
             dis24=diszy(1,j2,j4)
             dis25=diszy(1,j2,j5)
             dis35=diszy(1,j3,j5)
             jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       call smooth               !smooth the assignment

      return
      end

cccc make secondary str
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

cccc smooth the secondary structure assignment
      subroutine smooth
      parameter (maxres=3000)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

***   smooth single -------------->
***   --x-- => -----
       do i=3,nseq1
          if(isec(i).eq.2.or.isec(i).eq.4)then
             j=isec(i)
             if(isec(i-2).ne.j)then
                if(isec(i-1).ne.j)then
                   if(isec(i+1).ne.j)then
                      if(isec(i+1).ne.j)then
                         isec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo
       do i=3,nseq2
          if(jsec(i).eq.2.or.jsec(i).eq.4)then
             j=jsec(i)
             if(jsec(i-2).ne.j)then
                if(jsec(i-1).ne.j)then
                   if(jsec(i+1).ne.j)then
                      if(jsec(i+1).ne.j)then
                         jsec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo

***   smooth double -------------->
***   --xx-- => ------
       do i=1,nseq1-5
          if(isec(i).ne.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    if(isec(i+3).eq.2)then
                        if(isec(i+4).ne.2)then
                            if(isec(i+5).ne.2)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(isec(i).ne.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    if(isec(i+3).eq.4)then
                        if(isec(i+4).ne.4)then
                            if(isec(i+5).ne.4)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo
       do i=1,nseq2-5
          if(jsec(i).ne.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    if(jsec(i+3).eq.2)then
                        if(jsec(i+4).ne.2)then
                            if(jsec(i+5).ne.2)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(jsec(i).ne.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    if(jsec(i+3).eq.4)then
                        if(jsec(i+4).ne.4)then
                            if(jsec(i+5).ne.4)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo

***   connect -------------->
***   x-x => xxx
       do i=1,nseq1-2
          if(isec(i).eq.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    isec(i+1)=2
                endif
            endif
          endif
 
          if(isec(i).eq.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    isec(i+1)=4
                endif
            endif
          endif
       enddo
       do i=1,nseq2-2
          if(jsec(i).eq.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    jsec(i+1)=2
                endif
            endif
          endif
 
          if(jsec(i).eq.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    jsec(i+1)=4
                endif
            endif
          endif
       enddo
 
      return
      end

      
cccc distance calculation

      function diszy(i,i1,i2)
       parameter (maxres=3000)
       common /coord/ xa(3,maxres,0:1)

        diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
       return
      end

cccc reading pdb file
       subroutine readpdb (filen,np,nlen,aanam1,aanam2,xyz) ! filen: file name, np: 0 or 1

        parameter (maxres=3000)           ! no. of residues       
        parameter (npr=600)

        character*100 buffer,fnam,filen
        character*3 aanam1(-2:20),resa(maxres,npr)
        character*1 aanam2(-2:20),seqa(maxres,npr)
        dimension xyz(3,maxres,npr)

        common /pdbinfo/ ires1(maxres,npr),resa,seqa
        
         nii=60
         do while (filen(nii:nii).eq.' ')
             nii=nii-1
         enddo

          open (unit=10,file=filen(1:nii),status='old',iostat=ios)

         if (ios.ne.0) then
          write(*,*)'Error in opening PDB file ',filen(1:nii),'. Exiting
     & now'
          stop
         endif
          
          i=0
          do while(.true.) 
            read(10,91,end=101) buffer
                if (i.gt.0) then
                  if(buffer(1:6).eq.'ENDMDL') then
                    write(*,*) 'Warning: PDB contains multiple models, 
     &               only first model will be used'
                  endif
                  if(buffer(1:3).eq.'TER'.or.buffer(1:3).eq.'END')goto 101
                endif

                if ( (buffer(1:4).eq.'ATOM').or.
     &            (buffer(1:6).eq.'HETATM'.and.buffer(17:20).eq.'MSE') ) then
                 if(buffer(13:16).eq.' CA '.or. buffer(13:16).eq.'  CA'
     &              .or.buffer(13:16).eq.'CA  ') then
                   if(buffer(17:17).eq.' '.or.buffer(17:17).eq.'A') then
                     i=i+1
                   read(buffer,90)fnam,resa(i,np),fnam,ires1(i,np),fnam,
     &                              xyz(1,i,np),xyz(2,i,np),xyz(3,i,np)
                     do j=-2,20
                       if (resa(i,np).eq.aanam1(j)) seqa(i,np)=aanam2(j)
                     enddo
                   endif
                 endif
                endif

          enddo
  101  continue
   90  format(a17,a3,a2,i4,a4,3f8.3)
   91  format(a100)
       close (10)
       if (i.eq.0) then
        write(*,*)'Error in reading PDB file ',filen(1:nii),'. Exiting
     & now'
        stop
       endif
       nlen=i

       !! If PDB is not in the present dir. then rename pdb file name!!
       if ( index(filen(1:nii),'/') .ne. 0) then
        io=0
        kk=1
        do while (kk .ne. 0)
            kk=index(filen(1+io:),'/')
            io=io+kk
        enddo
        fnam=filen(io+1:nii)
        filen=fnam
       endif


       return
       end

cccc Instruction for running of the program
       subroutine instru
       
       call header
       write(*,*) 'Instructions for running Fr-TM-align program:'
       write(*,*)   
       write(*,*) 'Required arguments:'
       write(*,*)   
       write(*,*) 'Two structures. Aligns structure.pdb to target.pdb'
       write(*,*) '(By default, TM-score is normalized by the length', 
     &  'of the target.pdb)'
       write(*,*) './TMalign structure.pdb target.pdb'
       write(*,*)   
       write(*,*) 'Optional arguments:'
       write(*,*)   
       write(*,*)'-f : <0 or 1 > Default is 0. This option is to choose'
       write(*,*)'     slow or fast algorithm for structural alignment'
       write(*,*)'     0 - is slow algorithm is used and '
       write(*,*)'     1 - fast algorithm is used for alignment'
       write(*,*) ''
       write(*,*) '-o : Output the superposed coordinates in a file'
       write(*,*) '    ./TMalign structure.pdb target.pdb -o TM.sup'
       write(*,*) '    (To view the superimposed structures by rasmol:',
     &  './rasmol -script TM.sup)'
       write(*,*)   
       write(*,*) '-l : TM-score to be normalized by the given length'
       write(*,*) '    ./TMalign structure.pdb target.pdb -L 100'
       write(*,*)   
       write(*,*) '-m : <0 or 1 > If -m = 1, stores the rotation'
       write(*,*) '     matrix in trf.mat file. Default value is 0.'

       write(*,*) '    ./TMalign structure.pdb target.pdb -m 1'
       write(*,*) 
       write(*,*) '-a : <0 or 1 > If -a = 0, aligns only two structures, pdbid read from'
       write(*,*) '     standard input. If -a = 1, program compares a list of proteins'
       write(*,*) '     against one target PDB'
       write(*,*) '    ./TMalign -a 1 list target.pdb '
       write(*,*) ' '
    
       stop
       end

      subroutine header
       write(6,*)
       write(6,*)'****************************************************',
     &     '***********'
       write(6,*)'**                         Fr-TM-TMalign            ',
     &     '         **'
       write(6,*)'** A protein structutal alignment program based on f',
     &     'ragment  **'
       write(6,*)'** assembly and TM-score                            ',
     &     '         **'
       write(6,*)'**                                                  ',
     &     '         **'
       write(6,*)'** Comments and bugs email to :                     ',
     &     '         **'
       write(6,*)'** spandit3@mail.gatech.edu                         ',
     &     '         **'
       write(6,*)'****************************************************',
     &     '***********'
       
       write(6,*)

       end

cccc 
