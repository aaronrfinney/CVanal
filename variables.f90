module glob_var


!==================================================================
!
! Contains parameters and derived types
!
!==================================================================

  implicit none

!---Program parameters
 integer, parameter :: dp=selected_real_kind(15,307)
 integer, parameter :: limit=1000
 real(dp), parameter :: pi=3.14159265358979323846
 real(dp), parameter :: na=6.0221413e+23

 integer, parameter :: maxcoord=50   !Maximum atomic coordination

!---I/O channels
 integer, parameter :: un_log=300      !CVanal.log
 integer, parameter :: un_analinp=201  !CVanal.inp
 integer, parameter :: un_config=202   !CONFIG
 integer, parameter :: un_fld=203      !FIELD
 integer, parameter :: un_hist=204     !HISTORY
 integer, parameter :: un_xyz=302      !CVanal.xyz

!===Derived types===

!---History format
 type histfrm
    integer :: vrsn
    logical :: frm
 end type histfrm


!---Ca first shell
type firstshell
   integer :: ain
  integer :: nmol
  integer,dimension(10) :: imol
end type firstshell

!---Atom list data
 type atom_list
    integer :: ain          ! atom number
    integer :: mtyp         ! molecule type number
    integer :: imol         ! molecule number (running)
    integer :: iatp         ! atom_type number
    character(len=8):: atnm ! atom_type label
    real(dp) :: mass        ! weight of atom
    real(dp) :: chrg        ! charge of atom
    integer :: iclu         ! CV log
    integer :: aicl         ! atom type
 end type atom_list


!---Molecule list data
 type mol_list
    integer :: mtyp    ! molecule type number
    integer :: imtp    ! imtp'th molecule of this type
    integer :: n1      ! atom number of first atom in molecule
 end type mol_list


!---Molecule type data
 type moldefine
    integer :: nmtp                                     !# of molecular types
    integer :: nclmtp                                   !# of CV mol types
    integer,dimension(:),allocatable :: nmol            !# mols of each type
    integer,dimension(:),allocatable :: cntyp           !1 if a CV type
    integer,dimension(:),allocatable :: nmat            !# atoms for each type
    integer,dimension(:),allocatable :: nbnd            !# of bonds in each mol
    character(len=15),dimension(:),allocatable :: mname !name of each type
    character(len=15),dimension(:),allocatable :: cnnam !name of each CV type
    integer,dimension(:),allocatable :: nadd            !number of additive mols
    integer,dimension(:,:),allocatable :: iref          !ref atoms for -
                                                        ! axis_sense
    integer,dimension(:,:),allocatable :: imat          !atom index in each type
    integer,dimension(:,:,:),allocatable :: lbnd        !bonded atom-list -
                                                        ! for each type
 end type moldefine


!---Connection arrays
 type connection
    integer :: ncn                                  !# of connections
    integer, dimension(maxcoord) :: ain             !atom index at bond end
    integer, dimension(maxcoord) :: iatp            !atom type at bond end
    real(dp), dimension(maxcoord) :: dcn            !length of connection
 end type connection

end module glob_var





module variables
 use glob_var

 implicit none

!---Input variables
 integer :: nomit, nskip, nframes
 logical :: xyzcl
 real(dp) :: dcoord
 character(len=8), dimension(:), allocatable :: cnatnam
 type(histfrm) :: trjfm

 integer :: natms
 integer :: nuniq
 character(len=8), dimension(limit) :: auniq

 type(firstshell),dimension(:),allocatable :: cafs

 integer :: nstp, imcon
 real(dp) :: tstp
 real(dp), dimension(:,:), allocatable :: rxyz
 real(dp), dimension(9) :: cell,avcell

 integer :: numol

 type(atom_list),dimension(:),allocatable :: atlst
 type(atom_list),dimension(:),allocatable :: cvlst
 type(mol_list),dimension(:),allocatable :: mlst
 type(connection), dimension(:), allocatable :: conn
 type(moldefine) :: mdef

 integer :: nacn
 integer :: numconn
 real(dp), dimension(:,:), allocatable :: clxyz

end module variables
