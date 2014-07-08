#include <misc.h>
#include <params.h>

#if (defined SPMD)
      subroutine mpi_move_left(sbuf, rbuf, leng)
      use mpishorthand, only: mpir8,mpicom
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid,       only: myid_y,npr_y,iam,myid_z
      use spmd_dyn,     only: comm_y,lreq
      implicit none
      include 'mpif.h'
      real(r8) sbuf(*), rbuf(*)
!
      integer leng, status(mpi_status_size), irq, ierr
!
       lreq=MPI_REQUEST_NULL
      if(myid_y.ne.0)      call mpi_isend(sbuf,leng,mpir8,myid_y-1,10,comm_y,lreq(1),ierr)
!      if(myid_y.ne.0)      call mpi_send(sbuf,leng,mpir8,myid_y-1,10,comm_y,ierr)
      if(myid_y.ne.(npr_y-1)) call mpi_irecv( rbuf,leng,mpir8,myid_y+1,10,comm_y,lreq(2),ierr)
!      if(myid_y.ne.(npr_y-1)) call mpi_recv( rbuf,leng,mpir8,myid_y+1,10,comm_y,status,ierr)
      return        
      end     

      subroutine mpi_move_right(sbuf, rbuf, leng)
      use mpishorthand, only: mpir8,mpicom
      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid,       only: myid_y,npr_y,iam,myid_z
      use spmd_dyn,     only: comm_y,rreq
      implicit none
      include 'mpif.h'
      real(r8) sbuf(*), rbuf(*)
!
      integer leng, status(mpi_status_size), ierr
!
      rreq=MPI_REQUEST_NULL
      if(myid_y.ne.0)      call mpi_irecv( rbuf,leng,mpir8,myid_y-1,10+myid_z,comm_y,rreq(2),ierr)
!      if(myid_y.ne.0)      call mpi_recv( rbuf,leng,mpir8,myid_y-1,101,comm_y,status,ierr)
      if(myid_y.ne.(npr_y-1)) call mpi_isend(sbuf,leng,mpir8,myid_y+1,10+myid_z,comm_y,rreq(1),ierr)
!      if(myid_y.ne.(npr_y-1)) call mpi_send(sbuf,leng,mpir8,myid_y+1,101,comm_y,ierr)
      return        
      end     
      subroutine mpi_w_a(req)
      use pmgrid,       only: myid_y,npr_y
      implicit none
      include 'mpif.h'
      integer  req(2)
      integer status(mpi_status_size,2),ierr
!      if(myid_y.ne.0.and.myid_y.ne.(npr_y-1))  call mpi_waitall(n,req,status,ierr)
       call mpi_waitall(2,req,status,ierr)
!      if(myid_y.ne.0) call mpi_wait(req(1),status,ierr)
!      if(myid_y.ne.(npr_y-1)) call mpi_wait(req(2),status,ierr)
      return
      end
#endif
