module rxt_equations_mod
  use VAR_MOD, only : var_lim
  use RXT_MOD, only : rxt_lim, prd_lim, prd_limp1

  implicit none

  private

  public :: write_rxt_out_code

contains

  subroutine write_rxt_out_code ( &
       rxmcnt, &
       rxmap, &
       fixmap, &
       solsym, &
       fixsym, &
       prdcnt, &
       prdmap, &
       rxntot, &
       phtcnt, &
       outfile )

    use io,      only : temp_path

    implicit none

    integer, intent(in) :: rxmcnt(2)
    integer, intent(in) :: rxmap(rxt_lim,prd_lim+3,2)
    integer, intent(in) :: fixmap(var_lim,3,2)
    character(len=16), intent(in) :: solsym(var_lim)
    character(len=16), intent(in) :: fixsym(var_lim)
    integer, intent(in) :: phtcnt
    integer, intent(in) :: prdcnt
    integer, intent(in) :: rxntot
    integer, intent(in) :: prdmap(var_lim,prd_limp1)
    character(len=*), intent(in) :: outfile

    integer, parameter :: unitno = 33
    character(len=128) :: line
    character(len=64) :: mod_name
    logical ::  lexist
    integer :: pos

    inquire( file = trim( temp_path ) // outfile, exist = lexist )
    if( lexist ) then
       call system( 'rm ' // trim( temp_path ) // trim(outfile) )
    end if
    open( unit = unitno, file = trim( temp_path ) // trim(outfile) )

    pos = index(trim(outfile),'.F')
    mod_name = outfile(1:pos-1)

    line = ' '
    line(1:) = 'module '//trim(mod_name)
    write(unitno,100) trim(line)
    
    line = ' '
    line(3:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
    write(unitno,100) trim(line)

    line = ' '
    line(3:) = 'implicit none'
    write(unitno,100) trim(line)

    line = ' '
    line(3:) = 'private'
    write(unitno,100) trim(line)

    line = ' '
    line(3:) = 'public :: set_rates'
    write(unitno,100) trim(line)

    line = ' '
    line(1:) = 'contains'
    write(unitno,100) trim(line)

    line = ' '
    line(4:) = 'subroutine set_rates( rxt_rates, sol, ncol )'
    write(unitno,100) trim(line)

    line = ' '
    line(7:) = 'real(r8), intent(inout) :: rxt_rates(:,:,:)'
    write(unitno,100) trim(line)
    line = ' '
    line(7:) = 'real(r8), intent(in) :: sol(:,:,:)'
    write(unitno,100) trim(line)
    line = ' '
    line(7:) = 'integer, intent(in) :: ncol'
    write(unitno,100) trim(line)

    call write_rxt_equations ( &
         rxmcnt, &
         rxmap, &
         fixmap, &
         solsym, &
         fixsym, &
         prdcnt, &
         prdmap, &
         rxntot, &
         phtcnt, &
         unitno )

    line = ' '
    line(3:) = 'end subroutine set_rates'
    write(unitno,100) trim(line)

    line = ' '
    line(1:) = 'end module '//trim(mod_name)
    write(unitno,100) trim(line)

    close(unitno)

100 format(a)

  end subroutine write_rxt_out_code

  subroutine write_rxt_equations ( &
       rxmcnt, &
       rxmap, &
       fixmap, &
       solsym, &
       fixsym, &
       prdcnt, &
       prdmap, &
       rxntot, &
       phtcnt, unitno )

    implicit none

    integer, intent(in) :: rxmcnt(2)
    integer, intent(in) :: rxmap(rxt_lim,prd_lim+3,2)
    integer, intent(in) :: fixmap(var_lim,3,2)
    character(len=16), intent(in) :: solsym(var_lim)
    character(len=16), intent(in) :: fixsym(var_lim)
    integer, intent(in) :: phtcnt
    integer, intent(in) :: prdcnt
    integer, intent(in) :: rxntot
    integer, intent(in) :: prdmap(var_lim,prd_limp1)
    integer, intent(in) :: unitno

    character(len=80)  :: eq_piece
    character(len=80)  :: doc_piece
    character(len=6)   :: num
    character(len=120) :: eqline
    character(len=120) :: docline
    character(len=16)   :: symbol

    character(len=120) :: equations(rxntot)
    character(len=120) :: docs(rxntot)

    integer :: i,j,l
    integer :: rxno
    logical :: debug = .false.

    equations(:) =  ' '
    docs(:) =  ' '

    ! this is for case where all reactants are invariants
    do i = 1,prdcnt

       rxno = prdmap(i,1)
       write( num, '(i6)' ) rxno
       docline = 'rate_const'
!!$       eqline  = 'rxt_rates(:ncol,:,'//trim(num)//') = rxt_rates(:ncol,:,'//trim(num)//')'

       call get_fixed_reactants( fixmap, var_lim,  3, phtcnt, rxno, fixsym, doc_piece )

       docline = trim(docline)//trim(doc_piece)

!!$       equations(rxno) = trim(eqline)
       docs(rxno) = trim(docline)

    enddo

    do i = 1,2
       do j = 1,rxmcnt(i)

          rxno = rxmap(j,1,i)

          write(num, '(i6)' ) rxno
          docline = 'rate_const'
          eqline  = 'rxt_rates(:ncol,:,'//trim(num)//') = rxt_rates(:ncol,:,'//trim(num)//')'
          eq_piece = ' '
          doc_piece = ' '

          call get_fixed_reactants( fixmap, var_lim,  3, phtcnt, rxno, fixsym, doc_piece )

          do l = 2,i+1
             if( rxmap(j,l,i) == 0 ) then
                exit
             end if
             symbol = solsym(ABS(rxmap(j,l,i)))

             write(num,'(i6)') ABS(rxmap(j,l,i))
             eq_piece = trim(eq_piece)//'*sol(:ncol,:,' // trim(num) //')'
             doc_piece = trim(doc_piece)//'*' // trim(symbol)

          end do

          eqline = trim(eqline)//trim(eq_piece)
          docline = trim(docline)//trim(doc_piece)

          equations(rxno) = trim(eqline)
          docs(rxno) = trim(docline)

       enddo
    enddo

    do i = 1,rxntot
       write(unitno,'(a6,a120,a)') ' ',equations(i), ' ! '//trim(docs(i))
    enddo
    
    if (debug) then
       write(*,*) ' EQUATIONS : '
       do i = 1,rxntot
          write(*,'(i4,a120,a)') i, '  '//equations(i), ' ! '//trim(docs(i))
       enddo
    endif

  end subroutine write_rxt_equations

  subroutine get_fixed_reactants( &
       fixmap, &
       rowdim, &
       coldim, &
       phtcnt, &
       rxno, &
       fixsym, &
       doc_piece )

    use VAR_MOD, only : var_lim

    implicit none

    !-----------------------------------------------------------------------
    !     	... Dummy args
    !-----------------------------------------------------------------------
    integer, intent(in) ::  rowdim,        coldim,        phtcnt
    integer, intent(in) ::  fixmap(rowdim,coldim,2)
    integer, intent(in) :: rxno
    character(len=*), intent(in)  ::  fixsym(:)

    character(len=*), intent(out) ::  doc_piece

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  j, l, index
    character(len=16) :: symbol
    character(len=6) :: num

    integer  :: GET_INDEX
    integer  :: irx

    doc_piece = ' '

    irx = rxno
    if( rxno < phtcnt ) then
       irx = - rxno
    end if
    do j = 1,2
       index = GET_INDEX( fixmap(1,1,j), var_lim, 3, 1, irx )
       if( index /= 0 ) then
          do l = 2,3

             if( fixmap(index,l,j) == 0 ) then
                return
             end if

             symbol = fixsym(fixmap(index,l,j)) 
             write(num,'(i6)') fixmap(index,l,j)

             doc_piece = trim(doc_piece)//'*' // trim(symbol)

          end do
          exit
       end if
    end do

  end subroutine get_fixed_reactants

end module rxt_equations_mod
