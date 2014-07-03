
      subroutine HIST_HDR( hst_file_cnt, &
			   histout_cnt, &
                           histout_map, &
                           user_hst_names, &
                           hist_type, &
			   dyn_hst_fld_cnt, &
                           spcsym, &
                           spccnt, &
                           hetmap, &
                           usrmap, &
                           ptplen, &
                           filename, &
                           model  )
!-----------------------------------------------------------------------
!        ... Process all output history tape controls
!-----------------------------------------------------------------------

      use MASS_DIAGS
      use VAR_MOD, only : var_lim, hst_file_lim, class_prod_cnt, class_loss_cnt, clsmap
      use RXT_MOD, only : rxt_lim

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::      hst_file_cnt                        ! number of history files
      integer, intent(in)  ::      histout_cnt(20,2,hst_file_lim)      ! number of outputs in each catagory
      integer, intent(in)  ::      histout_map(1000,20,2,hst_file_lim) ! map of outputs
      integer, intent(in)  ::      dyn_hst_fld_cnt(2)
      integer, intent(in)  ::      spccnt(*)        ! number of symbols in each catagory
      integer, intent(in)  ::      hetmap(*)        ! wet dep map
      integer, intent(in)  ::      usrmap(*)        ! ext frc map
      integer, intent(out) ::      ptplen           ! total hist tape fields
      character(len=64), intent(in) :: hist_type         ! type of dyn hist tape ( short/long )
      character(len=16), intent(in)  :: user_hst_names(var_lim,4)
      character(len=16), intent(in)  :: spcsym(var_lim,*)     ! list of symbols
      character(len=*), intent(in)  :: filename              ! path/name of simulation data file
      character(len=*), intent(in)  :: model                 ! model name

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: inst = 1, avgr = 2

      integer  ::  i, j, m, n, typind, id_ox = 0
      integer  ::  file
      integer  ::  istat
      integer  ::  class, classno
      integer  ::  summ(2), sums(2)
      character(len=32) :: fld_name(4)
      character(len=16)  :: namtag(1000)
      character(len=4)  :: numa
      character(len=4)  :: oper_tag(2) = (/ 'inst', 'avrg' /)
      character(len=3)  :: num
      logical  ::  lexist
      
      CLOSE( 30 )
      OPEN( unit = 30, &
	    file = 'hist.h', &
	    status = 'replace', &
	    iostat = istat )
      if( istat /= 0 ) then
	 write(*,*) 'HIST_HDR: Failed to open hist.h'
	 stop
      end if
!-----------------------------------------------------------------------
!	... Locate index for ox
!-----------------------------------------------------------------------
      do i = 1,spccnt(6)
         if( spcsym(i,6) == 'OX' ) then
	    id_ox = i
	    exit
	 end if
      end do

      do typind = inst,avgr
         summ(typind) = SUM( histout_cnt(:2,typind,1) ) + SUM( histout_cnt(5:6,typind,1) ) &
                      + SUM( histout_cnt(8:11,typind,1) ) + histout_cnt(13,typind,1)
      end do
      write(30,'(''# define PMULTI '',i3)') summ(1) + summ(2)
      write(30,'(''# define PMULTA '',i3)') summ(2)
      sums(1) = histout_cnt(3,inst,1) + histout_cnt(4,inst,1) + histout_cnt(7,inst,1) &
				    + histout_cnt(12,inst,1)
      sums(2) = histout_cnt(3,avgr,1) + histout_cnt(4,avgr,1) + histout_cnt(12,avgr,1)
      write(30,'(''# define PSINGL '',i3)') sums(1) + sums(2)
      write(30,'(''# define PSINGLA '',i3)') sums(2)
      ptplen = summ(1) + sums(1) + summ(2) + sums(2)
      write(30,'(''# define PTPLEN '',i3)') ptplen
      i = SUM( histout_cnt(:4,inst,1) ) + SUM( histout_cnt(8:13,inst,1) ) &
        + SUM( histout_cnt(:4,avgr,1) ) + SUM( histout_cnt(8:13,avgr,1) )
      write(30,'(''# define HSTLEN '',i3)') i
      write(30,'(''# define HSTOFFSET '',i3)') SUM( histout_cnt(5:7,inst,1) )
      write(30,'(''# define HSTINSCNT '',i3)') SUM( histout_cnt(:17,inst,1) )
      write(30,'(''# define HSTINSCNTM '',i3)') &
              SUM( histout_cnt(:4,inst,1) ) &
            + SUM( histout_cnt(8:13,inst,1) )
      write(30,'(''# define HSTPHTCNT '',i3)') histout_cnt(8,inst,1)
      write(30,'(''# define HSTRXTCNT '',i3)') histout_cnt(9,inst,1)
      write(30,'(''# define HSTPHTCNTA '',i3)') histout_cnt(8,avgr,1)
      write(30,'(''# define HSTRXTCNTA '',i3)') histout_cnt(9,avgr,1)
      write(30,'(''# define HSTRXTIND '',i3)') SUM( histout_cnt(:4,inst,1) )
      write(30,'(''# define HSTUSRINDI '',i3)') SUM( histout_cnt(:4,inst,1) ) &
					     + SUM( histout_cnt(8:11,inst,1) )
      write(30,'(''# define HSTRXTINDA '',i3)') SUM( histout_cnt(:4,inst,1) ) &
                                            + SUM( histout_cnt(8:13,inst,1) ) &
                                            + SUM( histout_cnt(:4,avgr,1) )
      write(30,'(''# define HSTUSRINDA '',i3)') SUM( histout_cnt(:4,inst,1) ) &
					     + SUM( histout_cnt(8:13,inst,1) ) &
					     + SUM( histout_cnt(:11,avgr,1) )
      write(30,'(''# define HSTAVGCNT '',i3)') SUM( histout_cnt(:17,avgr,1) )
      write(30,'(''# define HSTXPTICNT '',i3)') histout_cnt(1,inst,1)
      write(30,'(''# define HSTXPTACNT '',i3)') histout_cnt(1,avgr,1)
      write(30,'(''# define HSTPCEICNT '',i3)') histout_cnt(2,inst,1)
      write(30,'(''# define HSTPCEACNT '',i3)') histout_cnt(2,avgr,1)
      write(30,'(''# define HSTSFLXICNT '',i3)') histout_cnt(3,inst,1)
      write(30,'(''# define HSTSFLXACNT '',i3)') histout_cnt(3,avgr,1)
      write(30,'(''# define HSTDVELICNT '',i3)') histout_cnt(4,inst,1)
      write(30,'(''# define HSTDVELACNT '',i3)') histout_cnt(4,avgr,1)
      write(30,'(''# define HSTWDEPICNT '',i3)') histout_cnt(10,inst,1)
      write(30,'(''# define HSTWDEPACNT '',i3)') histout_cnt(10,avgr,1)
      write(30,'(''# define HSTEXTICNT '',i3)') histout_cnt(11,inst,1)
      write(30,'(''# define HSTEXTACNT '',i3)') histout_cnt(11,avgr,1)
      write(30,'(''# define HSTUSRSICNT '',i3)') histout_cnt(12,inst,1)
      write(30,'(''# define HSTUSRMICNT '',i3)') histout_cnt(13,inst,1)
      write(30,'(''# define HSTUSRSACNT '',i3)') histout_cnt(12,avgr,1)
      write(30,'(''# define HSTUSRMACNT '',i3)') histout_cnt(13,avgr,1)
      write(30,'(''# define HSTPROD '',i3)') histout_cnt(14,inst,1)
      write(30,'(''# define HSTPRODA '',i3)') histout_cnt(14,avgr,1)
      write(30,'(''# define HSTLOSS '',i3)') histout_cnt(15,inst,1)
      write(30,'(''# define HSTLOSSA '',i3)') histout_cnt(15,avgr,1)
      write(30,'(''# define HSTPLCNT1 '',i3)') MAX( MAXVAL(class_prod_cnt(1,:)),MAXVAL(class_loss_cnt(1,:)) )
      write(30,'(''# define HSTPLCNT2 '',i3)') MAX( MAXVAL(class_prod_cnt(2,:)),MAXVAL(class_loss_cnt(2,:)) )
      write(30,'(''# define HSTPLCNT3 '',i3)') MAX( MAXVAL(class_prod_cnt(3,:)),MAXVAL(class_loss_cnt(3,:)) )
      write(30,'(''# define HSTPLCNT4 '',i3)') MAX( MAXVAL(class_prod_cnt(4,:)),MAXVAL(class_loss_cnt(4,:)) )
      write(30,'(''# define HSTPLCNT5 '',i3)') MAX( MAXVAL(class_prod_cnt(5,:)),MAXVAL(class_loss_cnt(5,:)) )
      write(30,'(''# define SW_MDI_CNT '',i3)') histout_cnt(16,inst,1)
      write(30,'(''# define SW_MDA_CNT '',i3)') histout_cnt(16,avgr,1)
      write(30,'(''# define SW_HINST '',i5)') SUM(histout_cnt(:4,inst,1)) &
                                              + SUM(histout_cnt(8:11,inst,1)) &
                                              + SUM(histout_cnt(14:17,inst,1))
      write(30,'(''# define SW_HTIMAV '',i5)') SUM(histout_cnt(:4,avgr,1)) &
                                               + SUM(histout_cnt(8:11,avgr,1)) &
                                               + SUM(histout_cnt(14:17,avgr,1))
      if( SUM(histout_cnt(14,:,1)+histout_cnt(15,:,1)) /= 0 ) then
         write(30,'(''# define HSTPRDLOSS'')')
      end if
      if( SUM( histout_cnt(17,:,1) ) /= 0 ) then
         write(30,'(''# define HSTMDIAGS'')')
      end if
      if( id_ox /= 0 ) then
         do i = 14,15
            if( ANY( histout_map(:histout_cnt(i,1,1),i,1,1) == id_ox ) .or. &
                ANY( histout_map(:histout_cnt(i,2,1),i,2,1) == id_ox ) ) then
	       if( i == 14 ) then
                  write(30,'(''# define OX_PROD '')')
	       else
                  write(30,'(''# define OX_LOSS '')')
	       end if
	    end if
         end do
      end if
      if( SUM( histout_cnt(12:13,inst,1) ) + SUM( histout_cnt(12:13,avgr,1) ) /= 0 ) then
         write(30,'(''# define HSTUSR'')')
      end if

!-----------------------------------------------------------------------
!	... Group history output flag
!-----------------------------------------------------------------------
      if( (histout_cnt(2,inst,1) + histout_cnt(2,avgr,1)) /= 0 ) then
	 write(30,'(''# define GRPHST'')')
      end if

      CLOSE(30)

      if( model == 'MOZART') then
!-----------------------------------------------------------------------
!        ... Write the history output tape information
!	     1. item count by category
!	     2. xported species
!	     3. "pce" species
!	     4. surface emissions
!	     5. deposition velocities
!	     6. washout rates
!	     7. "extraneous" forcing rates
!-----------------------------------------------------------------------
      OPEN( unit   = 30, &
            file   = TRIM(filename), &
            status = 'old', &
            position = 'append', &
	    iostat = istat )
      if( istat /= 0 ) then
         write(*,*) 'HIST_HDR: Failed to open '//TRIM(filename)
	 stop
      end if
      write(30,'(i4)') hst_file_cnt
      write(30,'(10i4)') histout_cnt(:,:,:hst_file_cnt)
      if( ptplen /= 0 ) then
file_loop : &
	 do file = 1,hst_file_cnt
output_type : &
	 do typind = 1,2
	    if( histout_cnt(1,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(1,typind,file)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(histout_map(j,1,typind,file),6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_VMR_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(1,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
		  end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(1,typind,file),1,typind,file)
	    end if

	    if( histout_cnt(2,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(2,typind,file)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(histout_map(j,2,typind,file),7) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_VMR_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(2,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
		  end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(2,typind,file),2,typind,file)
	    end if

	    if( histout_cnt(3,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(3,typind,file)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(histout_map(j,3,typind,file),6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_SRF_EMIS_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(3,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(3,typind,file),3,typind,file)
	    end if

	    if( histout_cnt(4,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(4,typind,file)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(histout_map(j,4,typind,file),6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_DEP_VEL_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(4,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(4,typind,file),4,typind,file)
	    end if

	    if( histout_cnt(8,typind,file) /= 0 ) then
	       do j = 1,histout_cnt(8,typind,file)
		  m = MOD( j-1,4 ) + 1
		  write(numa,'(i4)') histout_map(j,8,typind,file) + 1000
		  fld_name(m) = 'J_' // numa(2:4) // '_' // oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(8,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(8,typind,file),8,typind,file)
	    end if

	    if( histout_cnt(9,typind,file) /= 0 ) then
	       do j = 1,histout_cnt(9,typind,file)
		  m = MOD( j-1,4 ) + 1
		  write(numa,'(i4)') histout_map(j,9,typind,file) + 1000
		  fld_name(m) = 'R_' // numa(2:4) // '_' // oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(9,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(9,typind,file),9,typind,file)
	    end if

	    if( histout_cnt(10,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(10,typind,file)
		  m = MOD( j-1,4 ) + 1
		  i = hetmap(histout_map(j,10,typind,file))
	          write(fld_name(m),'(a)') TRIM( spcsym(i,6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_WET_DEP_RATE_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(10,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(10,typind,file),10,typind,file)
	    end if

	    if( histout_cnt(11,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(11,typind,file)
		  m = MOD( j-1,4 ) + 1
		  i = usrmap(histout_map(j,11,typind,file))
	          write(fld_name(m),'(a)') TRIM( spcsym(i,6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_EXT_FRC_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(11,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(11,typind,file),11,typind,file)
	    end if

	    if( histout_cnt(12,typind,file) /= 0 ) then
	       do j = 1,histout_cnt(12,typind,file)
		  m = MOD( j-1,4 ) + 1
	          fld_name(m) = TRIM( user_hst_names(j,2*(typind-1)+1) ) // '_' // oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(12,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	    end if

	    if( histout_cnt(13,typind,file) /= 0 ) then
	       do j = 1,histout_cnt(13,typind,file)
		  m = MOD( j-1,4 ) + 1
	          fld_name(m) = TRIM( user_hst_names(j,2*(typind-1)+2) ) // '_' // oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(13,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	    end if

	    if( histout_cnt(14,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(14,typind,file)
		  class   = histout_map(j,14,typind,file)/1000
		  classno = MOD( histout_map(j,14,typind,file),1000 )
		  n       = clsmap(classno,class,2)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(n,6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_PROD_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(14,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(14,typind,file),14,typind,file)
	    end if

	    if( histout_cnt(15,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(15,typind,file)
		  class   = histout_map(j,15,typind,file)/1000
		  classno = MOD( histout_map(j,15,typind,file),1000 )
		  n       = clsmap(classno,class,2)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(n,6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_LOSS_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(15,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(15,typind,file),15,typind,file)
	    end if

	    if( histout_cnt(16,typind,file) /= 0 ) then
	       do j = 1,histout_cnt(16,typind,file),8
	          fld_name(:) = ' '
		  do i = 1,4
		     select case (i)
			case( 1 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_ADV_'',a4)') oper_tag(typind)
		        case( 2 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_DPS_'',a4)') oper_tag(typind)
			case( 3 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_CNV_'',a4)') oper_tag(typind)
			case( 4 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_DIF_'',a4)') oper_tag(typind)
	                   write(30,'(4a32)') fld_name(:)
		     end select
	          end do
	          fld_name(:) = ' '
		  do i = 1,4
		     select case (i)
			case( 1 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_CHM_'',a4)') oper_tag(typind)
			case( 2 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_XFLX_'',a4)') oper_tag(typind)
			case( 3 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_YFLX_'',a4)') oper_tag(typind)
			case( 4 )
	                   write(fld_name(i),'(a)') TRIM( spcsym(histout_map(j,16,typind,file),6) )
	                   write(fld_name(i)(LEN_TRIM(fld_name(i))+1:),'(''_ZFLX_'',a4)') oper_tag(typind)
	                   write(30,'(4a32)') fld_name(:)
		     end select
	          end do
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(16,typind,file),16,typind,file)
	    end if

	    if( histout_cnt(17,typind,file) /= 0 ) then
	       fld_name(:) = ' '
	       do j = 1,histout_cnt(17,typind,file)
		  m = MOD( j-1,4 ) + 1
	          write(fld_name(m),'(a)') TRIM( spcsym(histout_map(j,17,typind,file),6) )
	          write(fld_name(m)(LEN_TRIM(fld_name(m))+1:),'(''_DRY_DEP_FLX_'',a4)') oper_tag(typind)
		  if( m == 4 .or. j == histout_cnt(17,typind,file) ) then
	             write(30,'(4a32)') fld_name(:m)
	             fld_name(:) = ' '
	          end if
	       end do
	       write(30,'(20i4)') histout_map(:histout_cnt(17,typind,file),17,typind,file)
	    end if
	 end do output_type
	 end do file_loop
      end if
!-----------------------------------------------------------------------
!        ... Write out the diagnostics
!-----------------------------------------------------------------------
      call MASS_DIAGS_SERIALIZE( 30 )

      CLOSE(30)

      end if

      end subroutine HIST_HDR
