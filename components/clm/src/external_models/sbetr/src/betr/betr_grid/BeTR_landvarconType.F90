module BeTR_landvarconType
  !DESCRIPTION
  !module for land type constants
implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_landvarcon_type
    integer,  public :: istsoil      !soil         landunit type (natural vegetation)
    integer,  public :: istcrop      !crop         landunit type
    integer,  public :: istice       !land ice     landunit type (glacier)
    integer,  public :: istice_mec   !land ice (multiple elevation classes) landunit type
    integer,  public :: istdlak      !deep lake    landunit type (now used for all lakes)
    integer,  public :: istwet       !wetland      landunit type (swamp, marsh, etc.)

    integer,  public :: isturb_MIN   !minimum urban type index
    integer,  public :: isturb_tbd   !urban tbd    landunit type
    integer,  public :: isturb_hd    !urban hd     landunit type
    integer,  public :: isturb_md    !urban md     landunit type
    integer,  public :: isturb_MAX   !maximum urban type index

    integer,  public :: max_lunit    !maximum value that lun%itype can have
  end type betr_landvarcon_type

  type(betr_landvarcon_type), public :: betr_landvarcon

end module BeTR_landvarconType
