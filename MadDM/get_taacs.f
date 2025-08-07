c-------------------------------------------------------------------------c
      subroutine get_taacs(xmin, xmax, nsteps, logscale)
c-------------------------------------------------------------------------c
c                                                                         c
c  This subroutine prints out the thermally averaged annihilation cross   c
c  section over a range of x values defined by xmin and xmax.  nsteps     c
c  is the number of steps that are calculated and logscale chooses        c
c  whether the scale will be logrithmic or linear.  (logscale = 1 for     c
c  log, logscale = 0 for linear)                                          c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'

c input parameters
      double precision xmin, xmax
      integer nsteps, logscale

c parameters used in this subroutine only
      double precision x
      character(50) output_file, pname, output_file2
      integer i_proc, pname_index, i, j, k, m
      Integer process_relic_name_indexes(size(process_names)), process_index
      double precision matrix_taacs(nsteps, size(process_names))
      double precision x_values(nsteps)

c create array of x values
      do i=1, nsteps
        x_values(i) = (xmax-xmin)*dble(i-1)/dble(nsteps-1) + xmin
      enddo
c loop thorugh the dm + coannihilator particles
      i_proc=1
      do x1=1,ndmparticles
        do x2=x1,ndmparticles
c loop through the annihilation processes involving dm or coannihilator in the initial state
          do x3=1,ann_nprocesses(x1,x2)
c save the index of the annihilation processes in the process_relic_name_indexes array
            call get_process_index(x1,x2,x3,1,process_index) 
            process_relic_name_indexes(i_proc) = process_index
c calculate taacs and save in matrix_taacs, for annihilation processes
c (NOTE: for self-annihilation processes the cross-section includes a factor of 1/2, which should be cancelled)
            do k=1, size(x_values)
               matrix_taacs(k, i_proc) =  taacs_coupled(x_values(k),x1,x2,x3,1)
            enddo
            i_proc = i_proc + 1
          enddo
c loop through the DM -> DM processes
          do x4=1, dm2dm_nprocesses(x1, x2)
c save the index of the DM -> DM processes in the process_relic_name_indexes array
            call get_process_index(x1,x2,x4,2,process_index) 
              process_relic_name_indexes(i_proc) = process_index
              do k=1, size(x_values)
                matrix_taacs(k, i_proc) = taacs_coupled(x_values(k),x1,x2,x4,2)
              enddo
              i_proc = i_proc + 1
          enddo
        enddo
      enddo
        
c write the output
c open the output file
      output_file = 'output/taacs.csv'
      open(unit=42, file=output_file)
      write(42,'(A)') '### List of thermally averaged cross-sections from MadDM'
      write(42,'(A)') '### x = mDM/T'
      write(42,'(A)') '### Process index, Process name, Process PDGs'
      do j=1,size(process_names)
        i_proc = process_relic_name_indexes(j)
        pname = trim(process_names(i_proc))
        write(42,'(A2,I3,2(",",:,A))') '# ',i_proc,pname,trim(process_pdgs(i_proc))
      enddo
      
      write(42, '(A2,*(I3,:,","))') 'x,',(process_relic_name_indexes(j), j=1,size(process_names))
      do m=1, nsteps
        write(42, '(*(E11.5,:,","))') x_values(m), (matrix_taacs(m, j), j=1,size(process_names))
      enddo
      close(42)
c      output_file = 'output/processes_taacs.csv'
c      open(unit=44, file=output_file)
c      do m=1, size(process_names)
c        write(44, *) process_relic_name_indexes(m), ' , ', process_names(process_relic_name_indexes(m))
c      enddo
c      close(44)
      return
      end