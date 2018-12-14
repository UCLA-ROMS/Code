      subroutine mpi_exchange8_tile (istr,iend,jstr,jend, A, nmaxA)
      implicit none
      include "mpif.h"
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4), parameter :: size_Z=16*(N+1), 
     &                       size_X=8*(N+1)*(Lm+4),
     &                                       size_E=8*(N+1)*(Mm+4)
      real(kind=8) sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &     rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &     sendW(size_E),                    sendE(size_E),
     &     recvW(size_E),                    recvE(size_E),
     &     sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &     rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)
      common /mess_buffers/ sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE,
     &                    sendW,recvW,                sendE,recvE,
     &                     sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE
C$OMP THREADPRIVATE(/mess_buffers/)
      integer(kind=4) inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, 
     &                                p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxA)
      integer(kind=4) i,j,k, kshft, indx, mess_count, ierr,
     &        req(16), comm(16), status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=(istr+iend-2*iwest)/(2*(iend-istr+1))
      jtg=(jstr+jend-2*jsouth)/(2*(jend-jstr+1))
      itg=8*itg
      jtg=8*jtg
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr.eq.iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+2, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend.eq.ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+1, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+4, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+3, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 6, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 5, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 8, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 7, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2).eq.0) then
          if (west_msg_exch.and.istr.eq.iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+1, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend.eq.ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+2, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2).eq.0) then
          if (south_msg_exch.and.jstr.eq.jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+3, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend.eq.jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+4, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &              p_SW, 5, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &              p_NE, 6, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &              p_SE, 7, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &              p_NW, 8, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i).gt.0) then
          mess_count=mess_count+1
          if (mess_count.lt.i) then
            comm(mess_count)=comm(i)
            req(mess_count)=req(i)
          endif
        endif
      enddo
      do while (mess_count.gt.0)
        call MPI_Waitany(mess_count, req, j, status, ierr)
        indx=comm(j)
        mess_count=mess_count-1
        do i=j,mess_count
          req(i)=req(i+1)
          comm(i)=comm(i+1)
        enddo
        if (indx.eq.1) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(iwest-2,j,k)=recvW(j-jmin       +kshft)
              A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.2) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(ieast+1,j,k)=recvE(j-jmin       +kshft)
              A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.3) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jsouth-2,k)=recvS(i-imin       +kshft)
              A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.4) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jnorth+1,k)=recvN(i-imin       +kshft)
              A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.5) then
          do k=1,nmaxA
            A(iwest-2,jsouth-2,k)=rv_SW(k        )
            A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
            A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
            A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
          enddo
        elseif (indx.eq.6) then
          do k=1,nmaxA
            A(ieast+1,jnorth+1,k)=rv_NE(k        )
            A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
            A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
            A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
          enddo
        elseif (indx.eq.7) then
          do k=1,nmaxA
            A(ieast+1,jsouth-2,k)=rv_SE(k        )
            A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
            A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
            A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
          enddo
        elseif (indx.eq.8) then
          do k=1,nmaxA
            A(iwest-2,jnorth+1,k)=rv_NW(k        )
            A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
            A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
            A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
          enddo
        endif
      enddo
      return
      end
      subroutine mpi_exchange8_2_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                                      B, nmaxB)
      implicit none
      include "mpif.h"
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4), parameter :: size_Z=16*(N+1), 
     &                       size_X=8*(N+1)*(Lm+4),
     &                                       size_E=8*(N+1)*(Mm+4)
      real(kind=8) sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &     rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &     sendW(size_E),                    sendE(size_E),
     &     recvW(size_E),                    recvE(size_E),
     &     sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &     rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)
      common /mess_buffers/ sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE,
     &                    sendW,recvW,                sendE,recvE,
     &                     sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE
C$OMP THREADPRIVATE(/mess_buffers/)
      integer(kind=4) inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, 
     &                                p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxB)
      integer(kind=4) i,j,k, kshft, indx, mess_count, ierr,
     &        req(16), comm(16), status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=(istr+iend-2*iwest)/(2*(iend-istr+1))
      jtg=(jstr+jend-2*jsouth)/(2*(jend-jstr+1))
      itg=8*itg
      jtg=8*jtg
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr.eq.iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+2, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend.eq.ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+1, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+4, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+3, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 6, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 5, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 8, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 7, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2).eq.0) then
          if (west_msg_exch.and.istr.eq.iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+1, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend.eq.ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+2, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2).eq.0) then
          if (south_msg_exch.and.jstr.eq.jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+3, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend.eq.jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+4, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &              p_SW, 5, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &              p_NE, 6, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &              p_SE, 7, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &              p_NW, 8, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i).gt.0) then
          mess_count=mess_count+1
          if (mess_count.lt.i) then
            comm(mess_count)=comm(i)
            req(mess_count)=req(i)
          endif
        endif
      enddo
      do while (mess_count.gt.0)
        call MPI_Waitany(mess_count, req, j, status, ierr)
        indx=comm(j)
        mess_count=mess_count-1
        do i=j,mess_count
          req(i)=req(i+1)
          comm(i)=comm(i+1)
        enddo
        if (indx.eq.1) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(iwest-2,j,k)=recvW(j-jmin       +kshft)
              A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(iwest-2,j,k)=recvW(j-jmin       +kshft)
              B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.2) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(ieast+1,j,k)=recvE(j-jmin       +kshft)
              A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(ieast+1,j,k)=recvE(j-jmin       +kshft)
              B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.3) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jsouth-2,k)=recvS(i-imin       +kshft)
              A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jsouth-2,k)=recvS(i-imin       +kshft)
              B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.4) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jnorth+1,k)=recvN(i-imin       +kshft)
              A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jnorth+1,k)=recvN(i-imin       +kshft)
              B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.5) then
          do k=1,nmaxA
            A(iwest-2,jsouth-2,k)=rv_SW(k        )
            A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
            A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
            A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
            B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
            B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
          enddo
        elseif (indx.eq.6) then
          do k=1,nmaxA
            A(ieast+1,jnorth+1,k)=rv_NE(k        )
            A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
            A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
            A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
            B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
            B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
          enddo
        elseif (indx.eq.7) then
          do k=1,nmaxA
            A(ieast+1,jsouth-2,k)=rv_SE(k        )
            A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
            A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
            A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
            B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
            B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
          enddo
        elseif (indx.eq.8) then
          do k=1,nmaxA
            A(iwest-2,jnorth+1,k)=rv_NW(k        )
            A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
            A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
            A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
            B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
            B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
          enddo
        endif
      enddo
      return
      end
      subroutine mpi_exchange8_3_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                            B, nmaxB, C, nmaxC)
      implicit none
      include "mpif.h"
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4), parameter :: size_Z=16*(N+1), 
     &                       size_X=8*(N+1)*(Lm+4),
     &                                       size_E=8*(N+1)*(Mm+4)
      real(kind=8) sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &     rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &     sendW(size_E),                    sendE(size_E),
     &     recvW(size_E),                    recvE(size_E),
     &     sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &     rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)
      common /mess_buffers/ sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE,
     &                    sendW,recvW,                sendE,recvE,
     &                     sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE
C$OMP THREADPRIVATE(/mess_buffers/)
      integer(kind=4) inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, 
     &                                p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxB)
      integer(kind=4) nmaxC
      real(kind=8) C(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxC)
      integer(kind=4) i,j,k, kshft, indx, mess_count, ierr,
     &        req(16), comm(16), status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      ksize=ksize+nmaxC
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=(istr+iend-2*iwest)/(2*(iend-istr+1))
      jtg=(jstr+jend-2*jsouth)/(2*(jend-jstr+1))
      itg=8*itg
      jtg=8*jtg
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr.eq.iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+2, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend.eq.ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+1, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+4, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+3, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 6, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 5, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 8, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 7, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2).eq.0) then
          if (west_msg_exch.and.istr.eq.iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=C(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=C(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+1, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend.eq.ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=C(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=C(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+2, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2).eq.0) then
          if (south_msg_exch.and.jstr.eq.jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=C(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=C(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+3, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend.eq.jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=C(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=C(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+4, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SW(k        +offset)=C(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxC+offset)=C(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxC+offset)=C(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxC+offset)=C(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &              p_SW, 5, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NE(k        +offset)=C(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxC+offset)=C(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxC+offset)=C(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxC+offset)=C(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &              p_NE, 6, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SE(k        +offset)=C(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxC+offset)=C(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxC+offset)=C(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxC+offset)=C(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &              p_SE, 7, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NW(k        +offset)=C(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxC+offset)=C(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxC+offset)=C(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxC+offset)=C(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &              p_NW, 8, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i).gt.0) then
          mess_count=mess_count+1
          if (mess_count.lt.i) then
            comm(mess_count)=comm(i)
            req(mess_count)=req(i)
          endif
        endif
      enddo
      do while (mess_count.gt.0)
        call MPI_Waitany(mess_count, req, j, status, ierr)
        indx=comm(j)
        mess_count=mess_count-1
        do i=j,mess_count
          req(i)=req(i+1)
          comm(i)=comm(i+1)
        enddo
        if (indx.eq.1) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(iwest-2,j,k)=recvW(j-jmin       +kshft)
              A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(iwest-2,j,k)=recvW(j-jmin       +kshft)
              B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxB
          do k=1,nmaxC
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              C(iwest-2,j,k)=recvW(j-jmin       +kshft)
              C(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.2) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(ieast+1,j,k)=recvE(j-jmin       +kshft)
              A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(ieast+1,j,k)=recvE(j-jmin       +kshft)
              B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxB
          do k=1,nmaxC
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              C(ieast+1,j,k)=recvE(j-jmin       +kshft)
              C(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.3) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jsouth-2,k)=recvS(i-imin       +kshft)
              A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jsouth-2,k)=recvS(i-imin       +kshft)
              B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxB
          do k=1,nmaxC
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              C(i,jsouth-2,k)=recvS(i-imin       +kshft)
              C(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.4) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jnorth+1,k)=recvN(i-imin       +kshft)
              A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jnorth+1,k)=recvN(i-imin       +kshft)
              B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxB
          do k=1,nmaxC
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              C(i,jnorth+1,k)=recvN(i-imin       +kshft)
              C(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.5) then
          do k=1,nmaxA
            A(iwest-2,jsouth-2,k)=rv_SW(k        )
            A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
            A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
            A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
            B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
            B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            C(iwest-1,jsouth-2,k)=rv_SW(k + nmaxC+offset)
            C(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxC+offset)
            C(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxC+offset)
          enddo
        elseif (indx.eq.6) then
          do k=1,nmaxA
            A(ieast+1,jnorth+1,k)=rv_NE(k        )
            A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
            A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
            A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
            B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
            B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            C(ieast+2,jnorth+1,k)=rv_NE(k + nmaxC+offset)
            C(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxC+offset)
            C(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxC+offset)
          enddo
        elseif (indx.eq.7) then
          do k=1,nmaxA
            A(ieast+1,jsouth-2,k)=rv_SE(k        )
            A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
            A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
            A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
            B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
            B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            C(ieast+2,jsouth-2,k)=rv_SE(k + nmaxC+offset)
            C(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxC+offset)
            C(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxC+offset)
          enddo
        elseif (indx.eq.8) then
          do k=1,nmaxA
            A(iwest-2,jnorth+1,k)=rv_NW(k        )
            A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
            A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
            A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
            B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
            B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            C(iwest-1,jnorth+1,k)=rv_NW(k + nmaxC+offset)
            C(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxC+offset)
            C(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxC+offset)
          enddo
        endif
      enddo
      return
      end
      subroutine mpi_exchange8_4_tile (istr,iend,jstr,jend, A, nmaxA,
     &                                  B, nmaxB, C, nmaxC, D, nmaxD)
      implicit none
      include "mpif.h"
      integer(kind=4), parameter ::
     &               LLm=435, MMm=660, N=60
      integer(kind=4), parameter ::
     &      NP_XI=8, NP_ETA=32, NSUB_X=1, NSUB_E=1
      integer(kind=4), parameter :: NNODES=NP_XI*NP_ETA,
     &    Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA
      integer(kind=4) ocean_grid_comm, mynode,  iSW_corn, jSW_corn,
     &                         iwest, ieast, jsouth, jnorth
      logical west_exchng,  east_exchng
      logical south_exchng, north_exchng
      common /mpi_comm_vars/  ocean_grid_comm, mynode,
     &     iSW_corn, jSW_corn, iwest, ieast, jsouth, jnorth
     &                , west_exchng,  east_exchng
     &                , south_exchng, north_exchng
      integer(kind=4), parameter :: padd_X=(Lm+2)/2-(Lm+1)/2,
     &                      padd_E=(Mm+2)/2-(Mm+1)/2
     &       , itemp=1
     &       , isalt=2
     &       , ntrc_salt=1
      integer(kind=4), parameter :: ntrc_pas=0
      integer(kind=4), parameter :: itrc_bio=itemp+ntrc_salt+ntrc_pas+1
      integer(kind=4), parameter :: iPO4=itrc_bio,iNO3=iPO4+1, 
     &                           iSIO3=iPO4+2,
     &     iNH4=iPO4+3,
     &     iFE=iPO4+4, iO2=iPO4+5, iDIC=iPO4+6,
     &     iALK=iPO4+7, iDOC=iPO4+8, iSPC=iPO4+9,
     &     iSPCHL=iPO4+10, iSPCACO3=iPO4+11, iDIATC=iPO4+12,
     &     iDIATCHL=iPO4+13, iZOOC=iPO4+14, iSPFE=iPO4+15,
     &     iDIATSI=iPO4+16, iDIATFE=iPO4+17, iDIAZC=iPO4+18,
     &     iDIAZCHL=iPO4+19, iDIAZFE=iPO4+20, iDON=iPO4+21,
     &     iDOFE=iPO4+22, iDOP=iPO4+23
      integer(kind=4), parameter :: iNO2 = iDOP + 1
      integer(kind=4), parameter :: iN2O = iNO2 + 1, iN2 = iN2O + 1
      integer(kind=4), parameter :: iSedOrgC = 1
      integer(kind=4), parameter :: iSedCaCO3 = iSedOrgC + 1
      integer(kind=4), parameter :: iSedSi = iSedCaCO3 + 1
      integer(kind=4), parameter :: NT_sed = 3
      integer(kind=4), parameter :: ntrc_bio = 24
     &     + 3
      integer(kind=4), parameter :: NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio
      integer(kind=4), parameter :: size_Z=16*(N+1), 
     &                       size_X=8*(N+1)*(Lm+4),
     &                                       size_E=8*(N+1)*(Mm+4)
      real(kind=8) sn_NW(size_Z),   sendN(size_X),   sn_NE(size_Z),
     &     rv_NW(size_Z),   recvN(size_X),   rv_NE(size_Z),
     &     sendW(size_E),                    sendE(size_E),
     &     recvW(size_E),                    recvE(size_E),
     &     sn_SW(size_Z),   sendS(size_X),   sn_SE(size_Z),
     &     rv_SW(size_Z),   recvS(size_X),   rv_SE(size_Z)
      common /mess_buffers/ sn_NW,rv_NW, sendN,recvN, sn_NE,rv_NE,
     &                    sendW,recvW,                sendE,recvE,
     &                     sn_SW,rv_SW, sendS,recvS, sn_SE,rv_SE
C$OMP THREADPRIVATE(/mess_buffers/)
      integer(kind=4) inode,jnode,  p_W, p_SW, p_S, p_SE, p_E, p_NE, 
     &                                p_N,
     &                                      p_NW,  exc_call_count
      logical west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      common /hidden_mpi_vars/ inode,jnode, p_W, p_SW, p_S, p_SE,
     &                      p_E, p_NE, p_N, p_NW,  exc_call_count,
     &        west_msg_exch,east_msg_exch, south_msg_exch,north_msg_exch
      save /hidden_mpi_vars/
      integer(kind=4) istr,iend,jstr,jend, nmaxA
      real(kind=8) A(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxA)
      integer(kind=4) nmaxB, offset
      real(kind=8) B(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxB)
      integer(kind=4) nmaxC
      real(kind=8) C(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxC)
      integer(kind=4) nmaxD
      real(kind=8) D(-1:Lm+2+padd_X,-1:Mm+2+padd_E,nmaxD)
      integer(kind=4) i,j,k, kshft, indx, mess_count, ierr,
     &        req(16), comm(16), status(MPI_STATUS_SIZE)
      integer(kind=4) ipass
      integer(kind=4) imin,imax,jmin,jmax, ishft,jshft,
     &          isize,jsize,ksize, itg,jtg
      if (istr.eq.iwest .and. .not.west_exchng) then
        imin=istr-1
      else
        imin=istr
      endif
      if (iend.eq.ieast .and. .not.east_exchng) then
        imax=iend+1
      else
        imax=iend
      endif
      if (jstr.eq.jsouth .and. .not.south_exchng) then
        jmin=jstr-1
      else
        jmin=jstr
      endif
      if (jend.eq.jnorth .and. .not.north_exchng) then
        jmax=jend+1
      else
        jmax=jend
      endif
      ishft=imax-imin+1
      jshft=jmax-jmin+1
      ksize=nmaxA
      ksize=ksize+nmaxB
      ksize=ksize+nmaxC
      ksize=ksize+nmaxD
      isize=2*ishft*ksize
      jsize=2*jshft*ksize
      ksize=4*ksize
      itg=(istr+iend-2*iwest)/(2*(iend-istr+1))
      jtg=(jstr+jend-2*jsouth)/(2*(jend-jstr+1))
      itg=8*itg
      jtg=8*jtg
      do i=1,16
        comm(i)=0
      enddo
      if (west_msg_exch.and.istr.eq.iwest) then
        call MPI_Irecv (recvW, jsize, MPI_DOUBLE_PRECISION,
     &          p_W, jtg+2, ocean_grid_comm, req(1), ierr)
        comm(1)=1
      endif
      if (east_msg_exch.and.iend.eq.ieast) then
        call MPI_Irecv (recvE, jsize, MPI_DOUBLE_PRECISION,
     &          p_E, jtg+1, ocean_grid_comm, req(2), ierr)
        comm(2)=2
      endif
      if (south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (recvS, isize, MPI_DOUBLE_PRECISION,
     &          p_S, itg+4, ocean_grid_comm, req(3), ierr)
        comm(3)=3
      endif
      if (north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (recvN, isize, MPI_DOUBLE_PRECISION,
     &          p_N, itg+3, ocean_grid_comm, req(4), ierr)
        comm(4)=4
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SW, ksize, MPI_DOUBLE_PRECISION,
     &             p_SW, 6, ocean_grid_comm, req(5), ierr)
        comm(5)=5
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NE, ksize, MPI_DOUBLE_PRECISION,
     &             p_NE, 5, ocean_grid_comm, req(6), ierr)
        comm(6)=6
      endif
      if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
        call MPI_Irecv (rv_SE, ksize, MPI_DOUBLE_PRECISION,
     &             p_SE, 8, ocean_grid_comm, req(7), ierr)
        comm(7)=7
      endif
      if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
        call MPI_Irecv (rv_NW, ksize, MPI_DOUBLE_PRECISION,
     &             p_NW, 7, ocean_grid_comm, req(8), ierr)
        comm(8)=8
      endif
      do ipass=0,1
        if (mod(inode+ipass,2).eq.0) then
          if (west_msg_exch.and.istr.eq.iwest) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=B(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=B(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=C(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=C(iwest+1,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=D(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=D(iwest+1,j,k)
              enddo
            enddo
            call MPI_Isend (sendW, jsize, MPI_DOUBLE_PRECISION,
     &              p_W, jtg+1, ocean_grid_comm, req(9), ierr)
            comm(9)=9
          endif
        else
          if (east_msg_exch.and.iend.eq.ieast) then
            do k=1,nmaxA
              kshft=2*jshft*(k-1) +1
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo
            offset=2*jshft*nmaxA +1
            do k=1,nmaxB
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=B(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=B(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxB
            do k=1,nmaxC
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=C(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=C(ieast  ,j,k)
              enddo
            enddo
            offset=offset + 2*jshft*nmaxC
            do k=1,nmaxD
              kshft=2*jshft*(k-1) +offset
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=D(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=D(ieast  ,j,k)
              enddo
            enddo
            call MPI_Isend (sendE, jsize, MPI_DOUBLE_PRECISION,
     &             p_E, jtg+2, ocean_grid_comm, req(10), ierr)
            comm(10)=10
          endif
        endif
      enddo
      do ipass=0,1
        if (mod(jnode+ipass,2).eq.0) then
          if (south_msg_exch.and.jstr.eq.jsouth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=B(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=B(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=C(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=C(i,jsouth+1,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendS(i-imin       +kshft)=D(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=D(i,jsouth+1,k)
              enddo
            enddo
            call MPI_Isend (sendS, isize, MPI_DOUBLE_PRECISION,
     &             p_S, itg+3, ocean_grid_comm, req(11), ierr)
            comm(11)=11
          endif
        else
          if (north_msg_exch.and.jend.eq.jnorth) then
            do k=1,nmaxA
              kshft=2*ishft*(k-1) +1
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo
            offset=2*ishft*nmaxA +1
            do k=1,nmaxB
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=B(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=B(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxB
            do k=1,nmaxC
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=C(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=C(i,jnorth  ,k)
              enddo
            enddo
            offset=offset + 2*ishft*nmaxC
            do k=1,nmaxD
              kshft=2*ishft*(k-1) +offset
              do i=imin,imax
                sendN(i-imin       +kshft)=D(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=D(i,jnorth  ,k)
              enddo
            enddo
            call MPI_Isend (sendN, isize, MPI_DOUBLE_PRECISION,
     &             p_N, itg+4, ocean_grid_comm, req(12), ierr)
            comm(12)=12
          endif
        endif
      enddo
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SW(k        )=A(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxA)=A(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxA)=A(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxA)=A(iwest+1,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SW(k        +offset)=B(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxB+offset)=B(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxB+offset)=B(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxB+offset)=B(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SW(k        +offset)=C(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxC+offset)=C(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxC+offset)=C(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxC+offset)=C(iwest+1,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SW(k        +offset)=D(iwest  ,jsouth  ,k)
            sn_SW(k + nmaxD+offset)=D(iwest+1,jsouth  ,k)
            sn_SW(k+2*nmaxD+offset)=D(iwest  ,jsouth+1,k)
            sn_SW(k+3*nmaxD+offset)=D(iwest+1,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SW, ksize, MPI_DOUBLE_PRECISION,
     &              p_SW, 5, ocean_grid_comm, req(13), ierr)
          comm(13)=13
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NE(k        )=A(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxA)=A(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxA)=A(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxA)=A(ieast  ,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NE(k        +offset)=B(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxB+offset)=B(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxB+offset)=B(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxB+offset)=B(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NE(k        +offset)=C(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxC+offset)=C(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxC+offset)=C(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxC+offset)=C(ieast  ,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NE(k        +offset)=D(ieast-1,jnorth-1,k)
            sn_NE(k + nmaxD+offset)=D(ieast  ,jnorth-1,k)
            sn_NE(k+2*nmaxD+offset)=D(ieast-1,jnorth  ,k)
            sn_NE(k+3*nmaxD+offset)=D(ieast  ,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NE, ksize, MPI_DOUBLE_PRECISION,
     &              p_NE, 6, ocean_grid_comm, req(14), ierr)
          comm(14)=14
        endif
        if (east_msg_exch.and.iend.eq.ieast .and. 
     &              south_msg_exch.and.jstr.eq.jsouth) then
          do k=1,nmaxA
            sn_SE(k        )=A(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxA)=A(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxA)=A(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxA)=A(ieast  ,jsouth+1,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_SE(k        +offset)=B(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxB+offset)=B(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxB+offset)=B(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxB+offset)=B(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_SE(k        +offset)=C(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxC+offset)=C(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxC+offset)=C(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxC+offset)=C(ieast  ,jsouth+1,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_SE(k        +offset)=D(ieast-1,jsouth  ,k)
            sn_SE(k + nmaxD+offset)=D(ieast  ,jsouth  ,k)
            sn_SE(k+2*nmaxD+offset)=D(ieast-1,jsouth+1,k)
            sn_SE(k+3*nmaxD+offset)=D(ieast  ,jsouth+1,k)
          enddo
          call MPI_Isend (sn_SE, ksize, MPI_DOUBLE_PRECISION,
     &              p_SE, 7, ocean_grid_comm, req(15), ierr)
          comm(15)=15
        endif
        if (west_msg_exch.and.istr.eq.iwest .and. 
     &              north_msg_exch.and.jend.eq.jnorth) then
          do k=1,nmaxA
            sn_NW(k        )=A(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxA)=A(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxA)=A(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxA)=A(iwest+1,jnorth  ,k)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            sn_NW(k        +offset)=B(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxB+offset)=B(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxB+offset)=B(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxB+offset)=B(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            sn_NW(k        +offset)=C(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxC+offset)=C(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxC+offset)=C(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxC+offset)=C(iwest+1,jnorth  ,k)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            sn_NW(k        +offset)=D(iwest  ,jnorth-1,k)
            sn_NW(k + nmaxD+offset)=D(iwest+1,jnorth-1,k)
            sn_NW(k+2*nmaxD+offset)=D(iwest  ,jnorth  ,k)
            sn_NW(k+3*nmaxD+offset)=D(iwest+1,jnorth  ,k)
          enddo
          call MPI_Isend (sn_NW, ksize, MPI_DOUBLE_PRECISION,
     &              p_NW, 8, ocean_grid_comm, req(16), ierr)
          comm(16)=16
        endif
      mess_count=0
      do i=1,16
        if (comm(i).gt.0) then
          mess_count=mess_count+1
          if (mess_count.lt.i) then
            comm(mess_count)=comm(i)
            req(mess_count)=req(i)
          endif
        endif
      enddo
      do while (mess_count.gt.0)
        call MPI_Waitany(mess_count, req, j, status, ierr)
        indx=comm(j)
        mess_count=mess_count-1
        do i=j,mess_count
          req(i)=req(i+1)
          comm(i)=comm(i+1)
        enddo
        if (indx.eq.1) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(iwest-2,j,k)=recvW(j-jmin       +kshft)
              A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(iwest-2,j,k)=recvW(j-jmin       +kshft)
              B(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxB
          do k=1,nmaxC
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              C(iwest-2,j,k)=recvW(j-jmin       +kshft)
              C(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxC
          do k=1,nmaxD
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              D(iwest-2,j,k)=recvW(j-jmin       +kshft)
              D(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.2) then
          do k=1,nmaxA
            kshft=2*jshft*(k-1) +1
            do j=jmin,jmax
              A(ieast+1,j,k)=recvE(j-jmin       +kshft)
              A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=2*jshft*nmaxA +1
          do k=1,nmaxB
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              B(ieast+1,j,k)=recvE(j-jmin       +kshft)
              B(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxB
          do k=1,nmaxC
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              C(ieast+1,j,k)=recvE(j-jmin       +kshft)
              C(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
          offset=offset + 2*jshft*nmaxC
          do k=1,nmaxD
            kshft=2*jshft*(k-1) +offset
            do j=jmin,jmax
              D(ieast+1,j,k)=recvE(j-jmin       +kshft)
              D(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
            enddo
          enddo
        elseif (indx.eq.3) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jsouth-2,k)=recvS(i-imin       +kshft)
              A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jsouth-2,k)=recvS(i-imin       +kshft)
              B(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxB
          do k=1,nmaxC
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              C(i,jsouth-2,k)=recvS(i-imin       +kshft)
              C(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxC
          do k=1,nmaxD
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              D(i,jsouth-2,k)=recvS(i-imin       +kshft)
              D(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.4) then
          do k=1,nmaxA
            kshft=2*ishft*(k-1) +1
            do i=imin,imax
              A(i,jnorth+1,k)=recvN(i-imin       +kshft)
              A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=2*ishft*nmaxA +1
          do k=1,nmaxB
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              B(i,jnorth+1,k)=recvN(i-imin       +kshft)
              B(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxB
          do k=1,nmaxC
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              C(i,jnorth+1,k)=recvN(i-imin       +kshft)
              C(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
          offset=offset + 2*ishft*nmaxC
          do k=1,nmaxD
            kshft=2*ishft*(k-1) +offset
            do i=imin,imax
              D(i,jnorth+1,k)=recvN(i-imin       +kshft)
              D(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
            enddo
          enddo
        elseif (indx.eq.5) then
          do k=1,nmaxA
            A(iwest-2,jsouth-2,k)=rv_SW(k        )
            A(iwest-1,jsouth-2,k)=rv_SW(k + nmaxA)
            A(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxA)
            A(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            B(iwest-1,jsouth-2,k)=rv_SW(k + nmaxB+offset)
            B(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxB+offset)
            B(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            C(iwest-1,jsouth-2,k)=rv_SW(k + nmaxC+offset)
            C(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxC+offset)
            C(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxC+offset)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            D(iwest-2,jsouth-2,k)=rv_SW(k        +offset)
            D(iwest-1,jsouth-2,k)=rv_SW(k + nmaxD+offset)
            D(iwest-2,jsouth-1,k)=rv_SW(k+2*nmaxD+offset)
            D(iwest-1,jsouth-1,k)=rv_SW(k+3*nmaxD+offset)
          enddo
        elseif (indx.eq.6) then
          do k=1,nmaxA
            A(ieast+1,jnorth+1,k)=rv_NE(k        )
            A(ieast+2,jnorth+1,k)=rv_NE(k + nmaxA)
            A(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxA)
            A(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            B(ieast+2,jnorth+1,k)=rv_NE(k + nmaxB+offset)
            B(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxB+offset)
            B(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            C(ieast+2,jnorth+1,k)=rv_NE(k + nmaxC+offset)
            C(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxC+offset)
            C(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxC+offset)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            D(ieast+1,jnorth+1,k)=rv_NE(k        +offset)
            D(ieast+2,jnorth+1,k)=rv_NE(k + nmaxD+offset)
            D(ieast+1,jnorth+2,k)=rv_NE(k+2*nmaxD+offset)
            D(ieast+2,jnorth+2,k)=rv_NE(k+3*nmaxD+offset)
          enddo
        elseif (indx.eq.7) then
          do k=1,nmaxA
            A(ieast+1,jsouth-2,k)=rv_SE(k        )
            A(ieast+2,jsouth-2,k)=rv_SE(k + nmaxA)
            A(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxA)
            A(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            B(ieast+2,jsouth-2,k)=rv_SE(k + nmaxB+offset)
            B(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxB+offset)
            B(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            C(ieast+2,jsouth-2,k)=rv_SE(k + nmaxC+offset)
            C(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxC+offset)
            C(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxC+offset)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            D(ieast+1,jsouth-2,k)=rv_SE(k        +offset)
            D(ieast+2,jsouth-2,k)=rv_SE(k + nmaxD+offset)
            D(ieast+1,jsouth-1,k)=rv_SE(k+2*nmaxD+offset)
            D(ieast+2,jsouth-1,k)=rv_SE(k+3*nmaxD+offset)
          enddo
        elseif (indx.eq.8) then
          do k=1,nmaxA
            A(iwest-2,jnorth+1,k)=rv_NW(k        )
            A(iwest-1,jnorth+1,k)=rv_NW(k + nmaxA)
            A(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxA)
            A(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxA)
          enddo
          offset=4*nmaxA
          do k=1,nmaxB
            B(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            B(iwest-1,jnorth+1,k)=rv_NW(k + nmaxB+offset)
            B(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxB+offset)
            B(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxB+offset)
          enddo
          offset=offset + 4*nmaxB
          do k=1,nmaxC
            C(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            C(iwest-1,jnorth+1,k)=rv_NW(k + nmaxC+offset)
            C(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxC+offset)
            C(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxC+offset)
          enddo
          offset=offset + 4*nmaxC
          do k=1,nmaxD
            D(iwest-2,jnorth+1,k)=rv_NW(k        +offset)
            D(iwest-1,jnorth+1,k)=rv_NW(k + nmaxD+offset)
            D(iwest-2,jnorth+2,k)=rv_NW(k+2*nmaxD+offset)
            D(iwest-1,jnorth+2,k)=rv_NW(k+3*nmaxD+offset)
          enddo
        endif
      enddo
      return
      end
