c
c     #################################################################
c     ##                                                             ##
c     ##  module mdiengine  --  MDI engine interface                 ##
c     ##                                                             ##
c     #################################################################
c
c
c     "init_mdi" checks for the argument to the -mdi command line
c     option and if it is found, initializes the MDI Library
c
c
      module mdiengine
 1    use mdi , only : MDI_NAME_LENGTH, MDI_COMMAND_LENGTH
      integer :: mdi_comm = 0
      logical :: mdi_exit = .false.
      logical :: use_mdi = .false.
      logical :: mdi_ignore_nodes = .false.
c     flag whether MDI has changed the value of use_ewald
      logical :: mdi_set_ewald = .false.
      logical :: forces_need_update = .true.
      logical :: analyze_need_update = .false.
      logical :: mdi_cycle_analyze = .false.
      character(len=MDI_NAME_LENGTH) :: target_node = "@DEFAULT"
      character(len=MDI_NAME_LENGTH) :: current_node = " "
      character(len=MDI_NAME_LENGTH) :: mdi_initial_caller = " "
      real*8, pointer, dimension (:,:)   :: forces_ptr
c
c     variables for communicating dipole information
c
      integer mdi_nprobes;
      integer, allocatable :: mdi_probes(:)
      integer, allocatable :: mdi_probe_mask(:)
      real*8, allocatable :: mdi_fielde(:,:)
      real*8, allocatable :: mdi_dfield_pair(:,:,:)
      real*8, allocatable :: mdi_ufield_pair(:,:,:)
c
c     Only used by analyze.x
c     If it is necessary to exit the node before responding to a command,
c       this variable will contain the current command
c
      character(len=MDI_COMMAND_LENGTH) :: current_command = " "
      save
      contains
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine init_mdi  --  Initialize the MDI Library            ##
c     ##                                                             ##
c     #################################################################
c
      subroutine init_mdi
      use argue
      use iounit
ccccccccccccc      use efield
      use mpole
 1    use iso_c_binding
 2    use mdi , only : MDI_Init, MDI_Accept_Communicator,
     &    MDI_Register_node, MDI_Register_command
      implicit none
      logical found_mdi
      integer i
      integer ierr
      character*240 mdi_options
      character*240 string
      procedure(execute_command), pointer :: generic_command => null()
      type(c_ptr)                         :: class_obj
      generic_command => execute_command
      class_obj = c_null_ptr
c
c     check for the --mdi command line argument
c
      use_mdi = .false.
      found_mdi = .false.
      do i = 1, narg-1
         string = arg(i)
         call upcase (string)
         if (string(1:5) .eq. '--MDI') then
            mdi_options = arg(i+1)
            found_mdi = .true.
         end if
      end do
c
c     initialize the MDI Library
c
      if ( found_mdi ) then
        use_mdi = .true.
        call MDI_Init(mdi_options, ierr)
        if ( ierr .ne. 0 ) then
           write(iout,*)'INIT_MDI -- Could not initalize MDI'
           call fatal
        end if
c
c     set the execute_command callback function
c
c        CALL MDI_Set_Execute_Command_Func(generic_command, class_obj,
c     &                                  ierr)
c
c     accept an MDI communicator
c
        call MDI_Accept_Communicator(mdi_comm, ierr)
        if ( ierr .ne. 0 ) then
           write(iout,*)'INIT_MDI -- Could not accept MDI communicator'
           call fatal
        end if
c
c     register all MDI nodes and commands
c
        call MDI_Register_node("@DEFAULT", ierr)
        call MDI_Register_command("@DEFAULT", "EXIT", ierr)
        call MDI_Register_command("@DEFAULT", ">ACTIVE", ierr)
        call MDI_Register_command("@DEFAULT", "<CELL", ierr)
        call MDI_Register_command("@DEFAULT", "<CELL_DISPL", ierr)
        call MDI_Register_command("@DEFAULT", "<CHARGES", ierr)
        call MDI_Register_command("@DEFAULT", "<COORDS", ierr)
        call MDI_Register_command("@DEFAULT", ">COORDS", ierr)
        call MDI_Register_command("@DEFAULT", "<DIMENSIONS", ierr)
        call MDI_Register_command("@DEFAULT", "<ELEMENTS", ierr)
        call MDI_Register_command("@DEFAULT", ">EWALD", ierr)
        call MDI_Register_command("@DEFAULT", "<MASSES", ierr)
        call MDI_Register_command("@DEFAULT", ">MASSES", ierr)
        call MDI_Register_command("@DEFAULT", "<NATOMS", ierr)
        call MDI_Register_command("@DEFAULT", "<POLARITIES", ierr)
        call MDI_Register_command("@DEFAULT", ">POLARITIES", ierr)
        call MDI_Register_command("@DEFAULT", ">POLARIZE", ierr)
        call MDI_Register_command("@DEFAULT", "<POLEDIMS", ierr)
        call MDI_Register_command("@DEFAULT", "<TOTCHARGE", ierr)
        call MDI_Register_command("@DEFAULT", "<NPOLES", ierr)
        call MDI_Register_command("@DEFAULT", "<POLES", ierr)
        call MDI_Register_command("@DEFAULT", "<IPOLES", ierr)
        call MDI_Register_command("@DEFAULT", "<FIELD", ierr)
        call MDI_Register_command("@DEFAULT", "<DFIELD", ierr)
        call MDI_Register_command("@DEFAULT", "<UFIELD", ierr)
        call MDI_Register_command("@DEFAULT", "<RESIDUES", ierr)
        call MDI_Register_command("@DEFAULT", "<MOLECULES", ierr)
        call MDI_Register_command("@DEFAULT", "<MULTIPOLES", ierr)
        call MDI_Register_command("@DEFAULT", ">MULTIPOLES", ierr)
        call MDI_Register_command("@DEFAULT", ">NPROBES", ierr)
        call MDI_Register_command("@DEFAULT", ">PROBES", ierr)
        call MDI_Register_command("@DEFAULT", "<@", ierr)
        call MDI_Register_command("@DEFAULT", "@INIT_MD", ierr)
c
c       Should only be supported if this was called from analyze.x
c
        if ( TRIM(mdi_initial_caller) .eq. 'analyze' ) THEN
           call MDI_Register_command("@DEFAULT", "@", ierr)
           call MDI_Register_command("@DEFAULT", "<ENERGY", ierr)
           call MDI_Register_command("@DEFAULT", "<PE", ierr)
           call MDI_Register_command("@DEFAULT", "<FORCES", ierr)
        end if
c
        call MDI_Register_node("@INIT_MD", ierr)
        call MDI_Register_command("@INIT_MD", "EXIT", ierr)
        call MDI_Register_command("@INIT_MD", "<CELL", ierr)
        call MDI_Register_command("@INIT_MD", "<CELL_DISPL", ierr)
        call MDI_Register_command("@INIT_MD", "<CHARGES", ierr)
        call MDI_Register_command("@INIT_MD", "<COORDS", ierr)
        call MDI_Register_command("@INIT_MD", ">COORDS", ierr)
        call MDI_Register_command("@INIT_MD", "<DIMENSIONS", ierr)
        call MDI_Register_command("@INIT_MD", "<ELEMENTS", ierr)
        call MDI_Register_command("@INIT_MD", "<MASSES", ierr)
        call MDI_Register_command("@INIT_MD", ">MASSES", ierr)
        call MDI_Register_command("@INIT_MD", "<NATOMS", ierr)
        call MDI_Register_command("@INIT_MD", "<TOTCHARGE", ierr)
        call MDI_Register_command("@INIT_MD", "<NPOLES", ierr)
        call MDI_Register_command("@INIT_MD", "<POLARITIES", ierr)
        call MDI_Register_command("@INIT_MD", ">POLARITIES", ierr)
        call MDI_Register_command("@INIT_MD", ">POLARIZE", ierr)
        call MDI_Register_command("@INIT_MD", "<POLEDIMS", ierr)
        call MDI_Register_command("@INIT_MD", "<POLES", ierr)
        call MDI_Register_command("@INIT_MD", "<IPOLES", ierr)
        call MDI_Register_command("@INIT_MD", "<FIELD", ierr)
        call MDI_Register_command("@INIT_MD", "<DFIELD", ierr)
        call MDI_Register_command("@INIT_MD", "<UFIELD", ierr)
        call MDI_Register_command("@INIT_MD", "<RESIDUES", ierr)
        call MDI_Register_command("@INIT_MD", "<MOLECULES", ierr)
        call MDI_Register_command("@INIT_MD", "<MULTIPOLES", ierr)
        call MDI_Register_command("@INIT_MD", ">MULTIPOLES", ierr)
        call MDI_Register_command("@INIT_MD", ">NPROBES", ierr)
        call MDI_Register_command("@INIT_MD", ">PROBES", ierr)
        call MDI_Register_command("@INIT_MD", "<@", ierr)
        call MDI_Register_command("@INIT_MD", "@", ierr)
        call MDI_Register_command("@INIT_MD", "@FORCES", ierr)
        call MDI_Register_command("@INIT_MD", "@COORDS", ierr)
c
        call MDI_Register_node("@FORCES", ierr)
        call MDI_Register_command("@FORCES", "EXIT", ierr)
        call MDI_Register_command("@FORCES", "<CELL", ierr)
        call MDI_Register_command("@FORCES", "<CELL_DISPL", ierr)
        call MDI_Register_command("@FORCES", "<CHARGES", ierr)
        call MDI_Register_command("@FORCES", "<COORDS", ierr)
        call MDI_Register_command("@FORCES", ">COORDS", ierr)
        call MDI_Register_command("@FORCES", "<DIMENSIONS", ierr)
        call MDI_Register_command("@FORCES", "<ELEMENTS", ierr)
        call MDI_Register_command("@FORCES", "<ENERGY", ierr)
        call MDI_Register_command("@FORCES", "<FORCES", ierr)
        call MDI_Register_command("@FORCES", ">FORCES", ierr)
        call MDI_Register_command("@FORCES", ">+FORCES", ierr)
        call MDI_Register_command("@FORCES", "<KE", ierr)
        call MDI_Register_command("@FORCES", "<MASSES", ierr)
        call MDI_Register_command("@FORCES", ">MASSES", ierr)
        call MDI_Register_command("@FORCES", "<NATOMS", ierr)
        call MDI_Register_command("@FORCES", "<PE", ierr)
        call MDI_Register_command("@FORCES", "<TOTCHARGE", ierr)
        call MDI_Register_command("@FORCES", "<NPOLES", ierr)
        call MDI_Register_command("@FORCES", "<POLARITIES", ierr)
        call MDI_Register_command("@FORCES", ">POLARITIES", ierr)
        call MDI_Register_command("@FORCES", ">POLARIZE", ierr)
        call MDI_Register_command("@FORCES", "<POLEDIMS", ierr)
        call MDI_Register_command("@FORCES", "<POLES", ierr)
        call MDI_Register_command("@FORCES", "<IPOLES", ierr)
        call MDI_Register_command("@FORCES", "<FIELD", ierr)
        call MDI_Register_command("@FORCES", "<DFIELD", ierr)
        call MDI_Register_command("@FORCES", "<UFIELD", ierr)
        call MDI_Register_command("@FORCES", "<RESIDUES", ierr)
        call MDI_Register_command("@FORCES", "<MOLECULES", ierr)
        call MDI_Register_command("@FORCES", "<MULTIPOLES", ierr)
        call MDI_Register_command("@FORCES", ">MULTIPOLES", ierr)
        call MDI_Register_command("@FORCES", ">NPROBES", ierr)
        call MDI_Register_command("@FORCES", ">PROBES", ierr)
        call MDI_Register_command("@FORCES", "<VELOCITIES", ierr)
        call MDI_Register_command("@FORCES", ">VELOCITIES", ierr)
        call MDI_Register_command("@FORCES", "<@", ierr)
        call MDI_Register_command("@FORCES", "@", ierr)
        call MDI_Register_command("@FORCES", "@FORCES", ierr)
        call MDI_Register_command("@FORCES", "@COORDS", ierr)
c
        call MDI_Register_node("@COORDS", ierr)
        call MDI_Register_command("@COORDS", "EXIT", ierr)
        call MDI_Register_command("@COORDS", "<CELL", ierr)
        call MDI_Register_command("@COORDS", "<CELL_DISPL", ierr)
        call MDI_Register_command("@COORDS", "<CHARGES", ierr)
        call MDI_Register_command("@COORDS", "<COORDS", ierr)
        call MDI_Register_command("@COORDS", ">COORDS", ierr)
        call MDI_Register_command("@COORDS", "<DIMENSIONS", ierr)
        call MDI_Register_command("@COORDS", "<ELEMENTS", ierr)
        call MDI_Register_command("@COORDS", "<ENERGY", ierr)
        call MDI_Register_command("@COORDS", "<FORCES", ierr)
        call MDI_Register_command("@COORDS", "<KE", ierr)
        call MDI_Register_command("@COORDS", "<MASSES", ierr)
        call MDI_Register_command("@COORDS", ">MASSES", ierr)
        call MDI_Register_command("@COORDS", "<NATOMS", ierr)
        call MDI_Register_command("@COORDS", "<PE", ierr)
        call MDI_Register_command("@COORDS", "<TOTCHARGE", ierr)
        call MDI_Register_command("@COORDS", "<NPOLES", ierr)
        call MDI_Register_command("@COORDS", "<POLARITIES", ierr)
        call MDI_Register_command("@COORDS", ">POLARITIES", ierr)
        call MDI_Register_command("@COORDS", ">POLARIZE", ierr)
        call MDI_Register_command("@COORDS", "<POLEDIMS", ierr)
        call MDI_Register_command("@COORDS", "<POLES", ierr)
        call MDI_Register_command("@COORDS", "<IPOLES", ierr)
        call MDI_Register_command("@COORDS", "<FIELD", ierr)
        call MDI_Register_command("@COORDS", "<DFIELD", ierr)
        call MDI_Register_command("@COORDS", "<UFIELD", ierr)
        call MDI_Register_command("@COORDS", "<RESIDUES", ierr)
        call MDI_Register_command("@COORDS", "<MOLECULES", ierr)
        call MDI_Register_command("@COORDS", "<MULTIPOLES", ierr)
        call MDI_Register_command("@COORDS", ">MULTIPOLES", ierr)
        call MDI_Register_command("@COORDS", ">NPROBES", ierr)
        call MDI_Register_command("@COORDS", ">PROBES", ierr)
        call MDI_Register_command("@FORCES", "<VELOCITIES", ierr)
        call MDI_Register_command("@FORCES", ">VELOCITIES", ierr)
        call MDI_Register_command("@COORDS", "<@", ierr)
        call MDI_Register_command("@COORDS", "@", ierr)
        call MDI_Register_command("@COORDS", "@FORCES", ierr)
        call MDI_Register_command("@COORDS", "@COORDS", ierr)
c
      end if
      mdi_exit = .false.
c
c     zero nprobes
c
ccccccccccccc      nprobes = 0

      return
      end subroutine init_mdi

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exit_mdi  --  Exit       MDI Library            ##
c     ##                                                             ##
c     #################################################################
c

      subroutine exit_mdi
ccccccccccccc      use efield
      use_mdi = .false.
      mdi_exit = .true.
ccccccccccccc      if (allocated(fielde)) deallocate (fielde)
ccccccccccccc      if (allocated(dfield_pair)) deallocate (dfield_pair)
ccccccccccccc      if (allocated(ufield_pair)) deallocate (ufield_pair)
      nprobes = 0
      return
      end subroutine exit_mdi


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdi_set_forces  --  Set a pointer to the foces  ##
c     ##                                                             ##
c     #################################################################
c
      subroutine mdi_set_forces(forces)
      use atoms , only  : n
      real*8, target, intent(in)         :: forces(3,n)
      forces_ptr => forces
      end subroutine mdi_set_forces

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdi_listen  --  Listen for commands             ##
c     ##                                                             ##
c     #################################################################
c
      subroutine mdi_listen(node_name)
 1    use mdi , only : MDI_NAME_LENGTH, MDI_Recv_Command
      implicit none
      character(len=*), intent(in) :: node_name
      integer ierr
      character(len=:), allocatable :: message
c
c     allocate memory to recieve an MDI command
c
      ALLOCATE( character(MDI_NAME_LENGTH) :: message )
      WRITE(6,*)'MDI at node: ',node_name
c
c     ignore this node if the mdi_ignore_nodes flag is set
c
      if ( mdi_ignore_nodes ) then
         return
      end if
c
c     set the current node
c
      current_node = node_name
c
c     do not listen if an "EXIT" command has been received
c
      if ( mdi_exit ) then
         return
      end if
c
c     check if this is the target node
c
      if ( target_node .ne. " " ) then
         if (target_node .eq. "@" .or. target_node .eq. node_name) then
            target_node = " "
         else
            return
         end if
      end if
      WRITE(6,*)'MDI listening at node: ',node_name
c
c     If using the analyze code, every time mdi_listen is called,
c     the energy, forces, etc. have been updated.
c
      if ( TRIM(mdi_initial_caller) .eq. 'analyze' ) then
         analyze_need_update = .false.
         if ( current_command .ne. " " ) then
c
c           respond to this command
c
            call execute_command(current_command, mdi_comm, ierr)
            mdi_cycle_analyze = .false.
            current_command = " "
         end if
      end if
c
c     listen for commands from the driver
c
      response_loop: do
c
c       receive a new command from the driver
c
        call MDI_Recv_Command(message, mdi_comm, ierr)
c
c       respond to this command
c
        call execute_command(message, mdi_comm, ierr)
c
c       check if the engine should stop listening for commands
c
        if ( target_node .ne. " " ) exit
        if ( mdi_exit ) exit
      end do response_loop
      return
      end subroutine mdi_listen
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdi_set_steps  --  Set the number of steps      ##
c     ##                                                             ##
c     #################################################################
c
      subroutine mdi_set_steps(istep, nstep)
      implicit none
      integer, intent(in) :: istep
      integer, intent(inout) :: nstep
      if ( mdi_exit ) then
c
c     if the EXIT command has been received, stop iterating
c
         nstep = istep
      else
c
c     otherwise, ensure that the md will continue iterating
c
         if ( nstep .le. istep ) then
            nstep = istep + 1
         end if
      end if
      return
      end subroutine mdi_set_steps
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine execute_command  --  Generic command response   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine execute_command(command, comm, ierr)
 1    use mdi , only : MDI_Send, MDI_CHAR, MDI_NAME_LENGTH,
     & MDI_Check_command_exists, MDI_COMM_NULL
      use iounit
      implicit none
      character(len=*), intent(in) :: command
      integer, intent(in)          :: comm
      integer, intent(out)         :: ierr
      integer                      :: cflag
      ierr = 0
c
c     confirm that this command is supported at this node
c
      call MDI_Check_command_exists(current_node, command,
     &                              MDI_COMM_NULL, cflag, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)
     &      'EXECUTE_COMMAND -- MDI_Check_command_exists failed'
         call fatal
      end if
      if ( cflag .ne. 1 ) then
         write(iout,*)'EXECUTE_COMMAND -- Unsupported command: ',
     &                command
         call fatal
      end if

c
c     respond to the command
c
      select case( TRIM(command) )
      case( "EXIT" )
        call exit_mdi
      case( ">ACTIVE" )
         call recv_active(comm)
      case( "<CELL" )
         call send_cell(comm)
      case( "<CELL_DISPL" )
         call send_celldispl(comm)
      case( "<CHARGES" )
         call send_charges(comm)
      case( "<COORDS" )
         call send_coords(comm)
      case( ">COORDS" )
         call recv_coords(comm)
      case( "<DIMENSIONS" )
         call send_dimensions(comm)
      case( "<ELEMENTS" )
         call send_elements(comm)
      case( "<ENERGY" )
         call send_energy(comm)
      case( ">EWALD" )
         call recv_ewald(comm)
      case( "<FORCES" )
         call send_forces(comm)
      case( ">FORCES" )
         call recv_forces(comm)
      case( ">+FORCES" )
         call add_forces(comm)
      case( "<KE" )
         call send_ke(comm)
      case( "<MASSES" )
         call send_masses(comm)
      case( ">MASSES" )
         call recv_masses(comm)
      case( "<NATOMS" )
         call send_natoms(comm)
      case( "<NPOLES" )
         call send_npoles(comm)
      case( "<PE" )
         call send_pe(comm)
      case( "<POLARITIES" )
         call send_polarities(comm)
      case( ">POLARITIES" )
         call recv_polarities(comm)
      case( ">POLARIZE" )
         call recv_polarize(comm)
      case( "<POLEDIMS" )
         call send_poledims(comm)
      case( "<TOTCHARGE" )
         call send_totcharge(comm)
      case( "<POLES" )
         call send_poles(comm)
      case( "<IPOLES" )
         call send_pole_indices(comm)
      case( "<FIELD" )
         call send_field(comm)
      case( "<DFIELD" )
         call send_dfield_components(comm)
      case( "<UFIELD" )
         call send_ufield_components(comm)
      case( "<RESIDUES" )
         call send_residues(comm)
      case( "<MOLECULES" )
         call send_molecules(comm)
      case( ">MULTIPOLES" )
         call recv_multipoles(comm)
      case( "<MULTIPOLES" )
         call send_poles(comm)
      case( ">NPROBES" )
         call recv_nprobes(comm)
      case( ">PROBES" )
         call recv_probes(comm)
      case( "<VELOCITIES" )
         call send_velocities(comm)
      case( ">VELOCITIES" )
         call recv_velocities(comm)
      case( "<@" )
         call MDI_Send(current_node, MDI_NAME_LENGTH, MDI_CHAR, comm,
     &                 ierr)
         if ( ierr .ne. 0 ) then
            write(iout,*)'EXECUTE_COMMAND -- MDI_Send failed'
            call fatal
         end if
      case( "@" )
         target_node = "@"
      case( "@INIT_MD" )
         target_node = "@INIT_MD"
      case( "@FORCES" )
         target_node = "@FORCES"
      case( "@COORDS" )
         target_node = "@COORDS"
      case default
        write(iout,*)'EXECUTE_COMMAND -- Command name not recognized: ',
     &                command
        call fatal
      end select
      return
      end subroutine execute_command
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_active  --  Respond to ">ACTIVE"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_active(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Recv
      use usage , only  : nuse, use, iuse
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, i, j
      integer, allocatable         :: mdiactive(:)

      allocate( mdiactive(n) )
c
c     receive the list of active atoms
c
      call MDI_Recv(mdiactive, n, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_ACTIVE -- MDI_Recv failed'
         call fatal
      end if
c
c     set the active atoms
c
ccccccccccccc      if (allocated(iuse))  deallocate (iuse)
ccccccccccccc      if (allocated(use))  deallocate (use)
      allocate (iuse(n))
      allocate (use(0:n))
      nuse = 0
      use(0) = .false.
      do i = 1, n
         if ( mdiactive(i) .eq. 0 ) then
           use(i) = .false.
         else
           use(i) = .true.
           nuse = nuse + 1
         end if
      end do
      j = 0
      do i = 1, n
         if (use(i)) then
            j = j + 1
            iuse(j) = i
         end if
      end do

      deallocate( mdiactive )
      return
      end subroutine recv_active
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_celldispl  --  Respond to "<CELL_DISPL"    ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_celldispl(comm)
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      real*8                       :: cell_displ(3)

c
c     construct the cell_displ array
c
      cell_displ = 0.0d0
c
c     send the cell_displ
c
      call MDI_Send(cell_displ, 3, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_CELLDISPL -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_celldispl
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_cell  --  Respond to "<CELL"               ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_cell(comm)
      use iounit , only : iout
      use boxes , only : lvec
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      real*8                       :: conv, cell(9)
c
c     get the conversion factor from angstrom to a.u.
c
      call MDI_Conversion_Factor("angstrom", "atomic_unit_of_length",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_CELL -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     construct the cell array
c
      cell = RESHAPE( lvec, SHAPE(cell) )
      cell = cell * conv
c
c     send the cell vector
c
      call MDI_Send(cell, 9, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_CELL -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_cell
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_charges  --  Respond to "<CHARGES"         ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_charges(comm)
      use atoms , only  : n
      use charge , only  : nion, iion, pchg
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send
c
c      use mpole
c
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: charges(:)

      allocate( charges(n) )
c
c     construct the charges array
c
      do iatom=1, n
         charges(iatom) = 0.0d0
      end do
      do iatom=1, nion
        charges(iion(iatom)) = charges(iion(iatom)) + pchg(iatom)
      end do
c
c     send the charges
c
      call MDI_Send(charges, n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_CHARGES -- MDI_Send failed'
         call fatal
      end if
      deallocate( charges )
      return
      end subroutine send_charges
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_coords  --  Respond to "<COORDS"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_coords(comm)
      use atoms , only  : n, x, y, z
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: coords(:)
      real*8                       :: conv

      allocate( coords(3*n) )
c
c     get the conversion factor from angstrom to a.u.
c
      call MDI_Conversion_Factor("angstrom", "atomic_unit_of_length",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_COORDS -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     construct the coordinates array
c
      do iatom=1, n
        coords( 3*(iatom-1) + 1 ) = x(iatom) * conv
        coords( 3*(iatom-1) + 2 ) = y(iatom) * conv
        coords( 3*(iatom-1) + 3 ) = z(iatom) * conv
      end do
c
c     send the coordinates
c
      call MDI_Send(coords, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_COORDS -- MDI_Send failed'
         call fatal
      end if
      deallocate( coords )
      return
      end subroutine send_coords
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_coords  --  Respond to ">COORDS"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_coords(comm)
      use atoms , only  : n, x, y, z
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: coords(:)
      real*8                       :: conv

      allocate( coords(3*n) )
c
c     get the conversion factor from a.u. to angstrom
c
      call MDI_Conversion_Factor("atomic_unit_of_length", "angstrom",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_COORDS -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     receive the coordinates
c
      call MDI_Recv(coords, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_COORDS -- MDI_Recv failed'
         call fatal
      end if
c
c     replace the system coords with the received coords
c
      do iatom=1, n
        x(iatom) = coords( 3*(iatom-1) + 1 ) * conv
        y(iatom) = coords( 3*(iatom-1) + 2 ) * conv
        z(iatom) = coords( 3*(iatom-1) + 3 ) * conv
      end do
c
c     the forces are now out-of-date
c
      forces_need_update = .true.
      analyze_need_update = .true.
      deallocate( coords )
      return
      end subroutine recv_coords
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_velocities  --  Respond to "<VELOCITIES"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_velocities(comm)
      use atoms , only  : n
      use moldyn , only : v
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdi_velocities(:)
      real*8                       :: conv_d, conv_t, conv

      allocate( mdi_velocities(3*n) )
c
c     get the conversion factor from angstrom to a.u.
c
      call MDI_Conversion_Factor("angstrom", "atomic_unit_of_length",
     &                           conv_d, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_VELOCITIES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the conversion factor from seconds to a.u.
c
      call MDI_Conversion_Factor("second", "atomic_unit_of_time",
     &                           conv_t, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_VELOCITIES -- MDI_Conversion_Factor failed'
         call fatal
      end if
      conv = conv_d / conv_t
c
c     construct the coordinates array
c
      do iatom=1, n
        mdi_velocities( 3*(iatom-1) + 1 ) = v(1, iatom) * conv
        mdi_velocities( 3*(iatom-1) + 2 ) = v(2, iatom) * conv
        mdi_velocities( 3*(iatom-1) + 3 ) = v(3, iatom) * conv
      end do
c
c     send the coordinates
c
      call MDI_Send(mdi_velocities, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_COORDS -- MDI_Send failed'
         call fatal
      end if
      deallocate( mdi_velocities )
      return
      end subroutine send_velocities
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_velocities  --  Respond to ">VELOCITIES"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_velocities(comm)
      use atoms , only  : n
      use moldyn , only : v
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdi_velocities(:)
      real*8                       :: conv_d, conv_t, conv

      allocate( mdi_velocities(3*n) )
c
c     get the conversion factor from angstrom to a.u.
c
      call MDI_Conversion_Factor("angstrom", "atomic_unit_of_length",
     &                           conv_d, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_VELOCITIES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the conversion factor from seconds to a.u.
c
      call MDI_Conversion_Factor("second", "atomic_unit_of_time",
     &                           conv_t, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_VELOCITIES -- MDI_Conversion_Factor failed'
         call fatal
      end if
      conv = conv_t / conv_d
c
c     receive the coordinates
c
      call MDI_Recv(mdi_velocities, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_VELOCITIES -- MDI_Recv failed'
         call fatal
      end if
c
c     replace the system velocities with the received velocities
c
      do iatom=1, n
        v(1, iatom) = mdi_velocities( 3*(iatom-1) + 1 ) * conv
        v(2, iatom) = mdi_velocities( 3*(iatom-1) + 2 ) * conv
        v(3, iatom) = mdi_velocities( 3*(iatom-1) + 3 ) * conv
      end do
      deallocate( mdi_velocities )
      return
      end subroutine recv_velocities
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_dimensions  --  Respond to "<DIMENSIONS"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_dimensions(comm)
      use boxes , only : orthogonal, monoclinic, triclinic, octahedron
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      integer                      :: dimensions(3)

c
c     construct the dimensions array
c
      dimensions = 1
      if ( orthogonal .or. monoclinic .or.
     &     triclinic .or. octahedron ) THEN
         dimensions = 2
      END IF
c
c     send the cell_displ
c
      call MDI_Send(dimensions, 3, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_DIMENSIONS -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_dimensions
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_elements  --  Respond to "<ELEMENTS"       ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_elements(comm)
      use atoms , only  : n
      use atmtyp , only  : atomic
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      integer, allocatable          :: mdielem(:)

      allocate( mdielem(n) )
c
c     construct the elements array
c
      do iatom=1, n
         mdielem(iatom) = atomic(iatom)
      end do
c
c     send the elements
c
      call MDI_Send(mdielem, n, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_ELEMENTS -- MDI_Send failed'
         call fatal
      end if
      deallocate( mdielem )
      return
      end subroutine send_elements
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_energy  --  Respond to "<ENERGY"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_energy(comm)
      use energi , only  : esum
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      real*8                       :: etotal, conv
c
c     variables for getting the kinetic energy
c
      real*8 eksum, temperature
      real*8 ekin(3,3)
c
c     check if this node needs to temporarily exit
c
      if ( TRIM(mdi_initial_caller) .eq. 'analyze' ) then
         if ( analyze_need_update ) then
            mdi_cycle_analyze = .true.
            current_command = '<ENERGY'
            target_node = '@DEFAULT'
            return
         end if
      end if
c
c     get the conversion factor from kilocalorie_per_mol to a.u.
c
      call MDI_Conversion_Factor("kilocalorie_per_mol",
     &                           "atomic_unit_of_energy",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_ENERGY -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the kinetic energy
c     there is no kinetic energy when using analyze.x
c
      eksum = 0.0d0
      if ( TRIM(mdi_initial_caller) .ne. 'analyze' ) THEN
         call kinetic(eksum, ekin, temperature)
      end if
c
c     get the total energy
c
      etotal = esum + eksum
c
c     convert the energy into atomic units
c
      etotal = etotal * conv
c
c     send the charges
c
      call MDI_Send(etotal, 1, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_ENERGY -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_energy
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_ewald  --  Respond to ">EWALD"             ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_ewald(comm)
ccccccccccccc      use iounit , only : iout
ccccccccccccc 1    use mdi , only    : MDI_INT, MDI_Recv
ccccccccccccc      use limits , only : use_ewald
ccccccccccccc      implicit none
      integer, intent(in)          :: comm
ccccccccccccc      integer                      :: ierr
ccccccccccccc      integer                      :: ewald_flag
cccccccccccccc
cccccccccccccc     receive the ewald flag
cccccccccccccc
ccccccccccccc      call MDI_Recv(ewald_flag, 1, MDI_INT, comm, ierr)
ccccccccccccc      if ( ierr .ne. 0 ) then
ccccccccccccc         write(iout,*)'SEND_EWALD -- MDI_Send failed'
ccccccccccccc         call fatal
ccccccccccccc      end if
cccccccccccccc
cccccccccccccc     it is not currently possible to turn on ewald if it was never in the keyfile
cccccccccccccc
ccccccccccccc      if ( .not. mdi_set_ewald ) then
ccccccccccccc         if ( .not. use_ewald ) then
ccccccccccccc            write(iout,*)'MDI ERROR: EWALD KEYWORD MISSING FROM KEYFILE'
ccccccccccccc            call fatal
ccccccccccccc         end if
ccccccccccccc      end if
cccccccccccccc
cccccccccccccc     set the correct value of use_ewald
cccccccccccccc
ccccccccccccc      if ( ewald_flag .eq. 0 ) then
ccccccccccccc        use_ewald = .false.
ccccccccccccc      else
ccccccccccccc        use_ewald = .true.
ccccccccccccc      end if
ccccccccccccc      mdi_set_ewald = .true.
cccccccccccccc
cccccccccccccc     redo cutoff initialization
cccccccccccccc
ccccccccccccc      call cutoffs
      return
      end subroutine recv_ewald
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_forces  --  Respond to "<FORCES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_forces(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdiforces(:)
      real*8                       :: lenconv, econv, conv

c
c     check if this node needs to temporarily exit
c
      if ( TRIM(mdi_initial_caller) .eq. 'analyze' ) then
         if ( analyze_need_update ) then
            mdi_cycle_analyze = .true.
            current_command = '<FORCES'
            target_node = '@DEFAULT'
            return
         end if
      end if
c
c     allocate an array for the forces
c
      allocate( mdiforces(3*n) )
c
c     get the conversion factor from angstrom to a.u.
c
      call MDI_Conversion_Factor("angstrom",
     &                           "atomic_unit_of_length",
     &                           lenconv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the conversion factor from kilocalorie_per_mol to a.u.
c
      call MDI_Conversion_Factor("kilocalorie_per_mol",
     &                           "atomic_unit_of_energy",
     &                           econv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
      conv = econv / lenconv
c
c     construct the coordinates array
c
      do iatom=1, n
        mdiforces( 3*(iatom-1) + 1 ) = forces_ptr(1,iatom) * conv
        mdiforces( 3*(iatom-1) + 2 ) = forces_ptr(2,iatom) * conv
        mdiforces( 3*(iatom-1) + 3 ) = forces_ptr(3,iatom) * conv
      end do
c
c     send the coordinates
c
      call MDI_Send(mdiforces, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_FORCES -- MDI_Send failed'
         call fatal
      end if
      deallocate( mdiforces )
      return
      end subroutine send_forces
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_forces  --  Respond to ">FORCES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_forces(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdiforces(:)
      real*8                       :: lenconv, econv, conv

      allocate( mdiforces(3*n) )
c
c     get the conversion factor from a.u. to angstrom
c
      call MDI_Conversion_Factor("atomic_unit_of_length", "angstrom",
     &                           lenconv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the conversion factor from a.u. to kilocalorie_per_mol
c
      call MDI_Conversion_Factor("atomic_unit_of_energy",
     &                           "kilocalorie_per_mol",
     &                           econv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
      conv = econv / lenconv
c
c     receive the coordinates
c
      call MDI_Recv(mdiforces, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Recv failed'
         call fatal
      end if
c
c     replace the system forces with the received forces
c
      do iatom=1, n
        forces_ptr(1,iatom) = mdiforces( 3*(iatom-1) + 1 ) * conv
        forces_ptr(2,iatom) = mdiforces( 3*(iatom-1) + 2 ) * conv
        forces_ptr(3,iatom) = mdiforces( 3*(iatom-1) + 3 ) * conv
      end do
      deallocate( mdiforces )
      return
      end subroutine recv_forces
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine add_forces  --  Respond to ">+FORCES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine add_forces(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdiforces(:)
      real*8                       :: lenconv, econv, conv

      allocate( mdiforces(3*n) )
c
c     get the conversion factor from a.u. to angstrom
c
      call MDI_Conversion_Factor("atomic_unit_of_length", "angstrom",
     &                           lenconv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the conversion factor from a.u. to kilocalorie_per_mol
c
      call MDI_Conversion_Factor("atomic_unit_of_energy",
     &                           "kilocalorie_per_mol",
     &                           econv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
      conv = econv / lenconv
c
c     receive the coordinates
c
      call MDI_Recv(mdiforces, 3*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Recv failed'
         call fatal
      end if
c
c     replace the system forces with the received forces
c
      do iatom=1, n
         forces_ptr(1,iatom) = forces_ptr(1,iatom) +
     &        mdiforces( 3*(iatom-1) + 1 ) * conv
         forces_ptr(2,iatom) = forces_ptr(2,iatom) +
     &        mdiforces( 3*(iatom-1) + 2 ) * conv
         forces_ptr(3,iatom) = forces_ptr(3,iatom) +
     &        mdiforces( 3*(iatom-1) + 3 ) * conv
      end do
      deallocate( mdiforces )
      return
      end subroutine add_forces
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_ke  --  Respond to "<KE"                   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_ke(comm)
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      real*8                       :: etotal, conv
c
c     variables for getting the kinetic energy
c
      real*8 eksum, temperature
      real*8 ekin(3,3)
c
c     get the conversion factor from kilocalorie_per_mol to a.u.
c
      call MDI_Conversion_Factor("kilocalorie_per_mol",
     &                           "atomic_unit_of_energy",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     get the kinetic energy
c
      call kinetic(eksum, ekin, temperature)
c
c     convert the energy into atomic units
c
      etotal = eksum * conv
c
c     send the charges
c
      call MDI_Send(etotal, 1, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_KE -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_ke
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_masses  --  Respond to "<MASSES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_masses(comm)
      use atoms , only  : n
      use atmtyp , only  : mass
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdimass(:)

      allocate( mdimass(n) )
c
c     construct the mass array
c
      do iatom=1, n
         mdimass(iatom) = mass(iatom)
      end do
c
c     send the mass
c
      call MDI_Send(mdimass, n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_MASSES -- MDI_Send failed'
         call fatal
      end if
      deallocate( mdimass )
      return
      end subroutine send_masses
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_masses  --  Respond to ">MASSES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_masses(comm)
      use atoms , only  : n
      use atmtyp , only  : mass
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8, allocatable          :: mdimass(:)

      allocate( mdimass(n) )
c
c     receive the mass
c
      call MDI_Recv(mdimass, n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_MASSES -- MDI_Recv failed'
         call fatal
      end if
c
c     construct the mass array
c
      do iatom=1, n
         mass(iatom) = mdimass(iatom)
      end do
      deallocate( mdimass )
      return
      end subroutine recv_masses
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_natoms  --  Respond to "<NATOMS"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_natoms(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
c
c     send the number of atoms
c
      call MDI_Send(n, 1, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_NATOMS -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_natoms
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_npoles  --  Respond to "<NPOLES"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_npoles(comm)
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Send
      use mpole , only  : npole
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
c
c     send the number of multipole sites
c
      call MDI_Send(npole, 1, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_NPOLES -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_npoles

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_nprobes --  Respond to ">NPROBES"          ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_nprobes(comm)
ccccccccccccc         use iounit , only : iout
ccccccccccccc1        use mdi , only    : MDI_DOUBLE, MDI_INT, MDI_Recv
ccccccccccccc         use efield , only : nprobes, probes, probe_mask
ccccccccccccc         use mpole , only  : npole
ccccccccccccc         implicit none
         integer, intent(in)          :: comm
ccccccccccccc         integer                      :: ierr, iprobe
ccccccccccccc
cccccccccccccc
cccccccccccccc     receive the number of probes
cccccccccccccc
ccccccccccccc         call MDI_Recv(nprobes, 1, MDI_INT, comm, ierr)
ccccccccccccc         if ( ierr .ne. 0 ) then
ccccccccccccc            write(iout,*)'RECV_ -- MDI_Recv failed'
ccccccccccccc            call fatal
ccccccccccccc         end if
cccccccccccccc
cccccccccccccc     allocate all arrays associated with the probes
cccccccccccccc
ccccccccccccc        if ( allocated( probes ) ) deallocate (probes)
ccccccccccccc        if ( allocated( probe_mask ) ) deallocate (probe_mask)
ccccccccccccc        allocate (probes(nprobes))
ccccccccccccc        allocate (probe_mask(npole))
ccccccccccccc
      return
      end subroutine recv_nprobes

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_pe  --  Respond to "<PE"                   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_pe(comm)
      use energi , only  : esum
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_factor
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      real*8                       :: etotal, conv
c
c     variables for getting the kinetic energy
c
      real*8 eksum, temperature
      real*8 ekin(3,3)
c
c     check if this node needs to temporarily exit
c
      if ( TRIM(mdi_initial_caller) .eq. 'analyze' ) then
         if ( analyze_need_update ) then
            mdi_cycle_analyze = .true.
            current_command = '<PE'
            target_node = '@DEFAULT'
            return
         end if
      end if
c
c     get the conversion factor from kilocalorie_per_mol to a.u.
c
      call MDI_Conversion_Factor("kilocalorie_per_mol",
     &                           "atomic_unit_of_energy",
     &                           conv, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_FORCES -- MDI_Conversion_Factor failed'
         call fatal
      end if
c
c     convert the energy into atomic units
c
      etotal = esum * conv
c
c     send the charges
c
      call MDI_Send(etotal, 1, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_ENERGY -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_pe
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_totcharge  --  Respond to "<TOTCHARGE"     ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_totcharge(comm)
      use atoms , only  : n
      use charge , only  : nion, pchg
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iatom
      real*8                       :: totcharge

c
c     construct the charges array
c
      totcharge = 0.0d0
      do iatom=1, nion
        totcharge = totcharge + pchg(iatom)
      end do
c
c     send the charges
c
      call MDI_Send(totcharge, 1, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_TOTCHARGE -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_totcharge
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_probes --  Respond to ">PROBES"            ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_probes(comm)
ccccccccccccc        use iounit , only : iout
ccccccccccccc1       use mdi , only    : MDI_DOUBLE, MDI_INT, MDI_Recv
ccccccccccccc        use efield , only : nprobes, probes, probe_mask
ccccccccccccc        use mpole , only : npole
ccccccccccccc        implicit none
        integer, intent(in)          :: comm
ccccccccccccc        integer                      :: ierr, iprobe, i, count
ccccccccccccc
ccccccccccccc        if ( nprobes .eq. 0 ) then
ccccccccccccc           write(iout,*)'Must receive >NPROBES before >PROBES'
ccccccccccccc           call fatal
ccccccccccccc        end if
ccccccccccccc        if ( .not. allocated( probes ) ) then
ccccccccccccc           write(iout,*)'>NPROBES: probes not allocated'
ccccccccccccc           call fatal
ccccccccccccc        end if
ccccccccccccc        if ( .not. allocated( probe_mask ) ) then
ccccccccccccc           write(iout,*)'>NPROBES: probe_mask not allocated'
ccccccccccccc           call fatal
ccccccccccccc        end if
cccccccccccccc
cccccccccccccc     receive the probes
cccccccccccccc
ccccccccccccc       call MDI_Recv(probes, nprobes, MDI_INT, comm, ierr)
ccccccccccccc       if ( ierr .ne. 0 ) then
ccccccccccccc          write(iout,*)'RECV_ -- MDI_Recv failed'
ccccccccccccc          call fatal
ccccccccccccc       end if
ccccccccccccc
cccccccccccccc
cccccccccccccc     Make probe mask
cccccccccccccc
ccccccccccccc      do i=1, npole
ccccccccccccc        probe_mask(i) = 0
ccccccccccccc      end do
ccccccccccccc
ccccccccccccc      count = 1
ccccccccccccc      do i=1, nprobes
ccccccccccccc        probe_mask(probes(i)) = count
ccccccccccccc        count = count + 1
ccccccccccccc      end do
ccccccccccccc
cccccccccccccc      write (*,*) probe_mask
ccccccccccccc
      return
      end subroutine recv_probes


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_polarize  --  Respond to ">POLARIZE"       ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_polarize(comm)
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Recv
      use potent , only  : use_polar
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr
      integer                      :: polarizability
c
c     recv the polarizability flag
c
      call MDI_Recv(polarizability, 1, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_POLARIZE -- MDI_Recv failed'
         call fatal
      end if
c
c     set the flag that controls use of the atomic dipole polarization PE term
c
      if ( polarizability .eq. 0 ) then
         use_polar = .false.
      else
         use_polar = .true.
      end if
      forces_need_update = .true.
      analyze_need_update = .true.
      return
      end subroutine recv_polarize


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_poledims  --  Respond to "<POLEDIMS"           ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_poledims(comm)
      use iounit , only : iout
 1    use mdi , only    : MDI_INT, MDI_Send
      use mpole , only  : maxpole
      implicit none
      integer, intent(in)          :: comm
      integer                      :: poledims, ierr
c
c     determine the dimensionality of the multipoles
c
      if ( maxpole .eq. 1 ) then
         poledims = 1
      else if ( maxpole .eq. 4 ) then
         poledims = 2
      else if ( maxpole .eq. 13 ) then
         poledims = 3
      else
         write(iout,*)'SEND_POLEDIMS -- Invalid maxpole'
         call fatal
      end if
c
c     send the dimensionality of the multipoles
c
      call MDI_Send(poledims, 1, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_POLEDIMS -- MDI_Send failed'
         call fatal
      end if
      return
      end subroutine send_poledims

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_poles  --  Respond to "<POLES"             ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_poles(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_Factor
      use mpole , only  : maxpole, npole, pole, ipole
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iipole, iatom, icomp
      integer                      :: stride
      real*8, allocatable          :: poles_buf(:)
      real*8                       :: conv
c
c     determine the stride of the multipole data
c
      if ( maxpole .eq. 1 ) then
         stride = 1
      else if ( maxpole .eq. 4 ) then
         stride = 4
      else if ( maxpole .eq. 13 ) then
         stride = 9
      end if
c
c     prepare the poles buffer
c
      allocate( poles_buf(stride*n) )
      poles_buf = 0.0d0
      do iipole=1, npole
         iatom = ipole(iipole)
c
c        monopole / charge
c
         poles_buf(stride*(iatom-1) + 1) = pole(1, iipole)
         if ( stride .gt. 1 ) then
c
c           dipole dx, dy, and dz terms
c
            poles_buf(stride*(iatom-1) + 2) = pole(2, iipole)
            poles_buf(stride*(iatom-1) + 3) = pole(3, iipole)
            poles_buf(stride*(iatom-1) + 4) = pole(4, iipole)
         end if
         if ( stride .gt. 4 ) then
c
c           quadrupole qxx, qxy, qxz, qyy, qyz terms
c           Note: the qzz term is not necessary for traceless quadrupoles
c              e.g., qxx + qyy + qzz = 0
c
            poles_buf(stride*(iatom-1) + 5) = pole(5, iipole)
            poles_buf(stride*(iatom-1) + 6) = pole(6, iipole)
            poles_buf(stride*(iatom-1) + 7) = pole(7, iipole)
            poles_buf(stride*(iatom-1) + 8) = pole(9, iipole)
            poles_buf(stride*(iatom-1) + 9) = pole(10, iipole)
         end if
      end do
c
c     send the poles
c
      call MDI_Send(poles_buf, stride*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_POLES -- MDI_Send failed'
         call fatal
      end if
      deallocate( poles_buf )
      return
      end subroutine send_poles


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_multipoles  --  Respond to ">MULTIPOLES"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_multipoles(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      use mpole , only  : maxpole, npole, pole, ipole
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iipole, stride, iatom
      real*8, allocatable          :: poles_buf(:)
      real*8                       :: conv
c
c     determine the stride of the multipole data
c
      if ( maxpole .eq. 1 ) then
         stride = 1
      else if ( maxpole .eq. 4 ) then
         stride = 4
      else if ( maxpole .eq. 13 ) then
         stride = 9
      end if
c
c     prepare the poles buffer
c
      allocate( poles_buf(stride*n) )
c
c     send the poles
c
      call MDI_Recv(poles_buf, stride*n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_MULTIPOLES -- MDI_Recv failed'
         call fatal
      end if
c
c     set rpole
c
      do iipole=1, npole
         iatom = ipole(iipole)
c
c        monopole / charge
c
         pole(1, iipole) = poles_buf(stride*(iatom-1) + 1)
         if ( stride .gt. 1 ) then
c
c           dipole dx, dy, and dz terms
c
            pole(2, iipole) = poles_buf(stride*(iatom-1) + 2)
            pole(3, iipole) = poles_buf(stride*(iatom-1) + 3)
            pole(4, iipole) = poles_buf(stride*(iatom-1) + 4)
         end if
         if ( stride .gt. 4 ) then
c
c           quadrupole qxx, qxy, qxz, qyy, qyz terms
c           Note: the qzz term is not necessary for traceless quadrupoles
c              e.g., qxx + qyy + qzz = 0
c
            pole(5, iipole) = poles_buf(stride*(iatom-1) + 5)
            pole(6, iipole) = poles_buf(stride*(iatom-1) + 6)
            pole(7, iipole) = poles_buf(stride*(iatom-1) + 7)
            pole(9, iipole) = poles_buf(stride*(iatom-1) + 8)
            pole(10, iipole) = poles_buf(stride*(iatom-1) + 9)
c
c           add the symmetry terms
c
            pole(8, iipole) = poles_buf(stride*(iatom-1) + 6)
            pole(11, iipole) = poles_buf(stride*(iatom-1) + 7)
            pole(12, iipole) = poles_buf(stride*(iatom-1) + 9)
c
c           add the qzz term from the traceless property
c
            pole(13,iipole) = -1.0*(pole(5,iipole) + pole(9,iipole))
         end if
      end do
      deallocate( poles_buf )
      forces_need_update = .true.
      analyze_need_update = .true.
      return
      end subroutine recv_multipoles


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_polarities  --  Respond to "<POLARITIES"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_polarities(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Send, MDI_Conversion_Factor
      use mpole , only  : maxpole, npole, ipole
      use polar , only  : polarity
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iipole, iatom
      real*8, allocatable          :: poles_buf(:)
      real*8                       :: conv
c
c     prepare the poles buffer
c
      allocate( poles_buf(n) )
      poles_buf = 0.0d0
c
c     get the polarities
c
      do iipole=1, npole
         iatom = ipole(iipole)
         poles_buf(iatom) = polarity(iipole)
      end do
c
c     send the poles
c
      call MDI_Send(poles_buf, n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_POLARITIES -- MDI_Send failed'
         call fatal
      end if
      deallocate( poles_buf )
      return
      end subroutine send_polarities


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine recv_polarities  --  Respond to ">POLARITIES"   ##
c     ##                                                             ##
c     #################################################################
c
      subroutine recv_polarities(comm)
      use atoms , only  : n
      use iounit , only : iout
 1    use mdi , only    : MDI_DOUBLE, MDI_Recv, MDI_Conversion_Factor
      use mpole , only  : maxpole, npole, ipole
      use polar , only  : polarity
      implicit none
      integer, intent(in)          :: comm
      integer                      :: ierr, iipole, iatom
      real*8, allocatable          :: poles_buf(:)
      real*8                       :: conv
c
c     prepare the poles buffer
c
      allocate( poles_buf(n) )
c
c     send the poles
c
      call MDI_Recv(poles_buf, n, MDI_DOUBLE, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'RECV_POLARITIES -- MDI_Recv failed'
         call fatal
      end if
c
c     set polarity
c
      do iipole=1, npole
         iatom = ipole(iipole)
         polarity(iipole) = poles_buf(iatom)
      end do
      deallocate( poles_buf )
      forces_need_update = .true.
      analyze_need_update = .true.
      return
      end subroutine recv_polarities


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_residues  --  Respond to "<RESIDUES"       ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_residues(comm)
ccccccccccccc      use iounit , only : iout
ccccccccccccc 1    use mdi , only    : MDI_INT, MDI_Send, MDI_Conversion_Factor
ccccccccccccc      use pdb , only  : resnum
ccccccccccccc      use atoms , only : n
ccccccccccccc      implicit none
      integer, intent(in)          :: comm
ccccccccccccc      integer                      :: ierr, iatom
ccccccccccccc      integer, allocatable         :: res_buf(:)
ccccccccccccc
cccccccccccccc
cccccccccccccc     if residues are not used in the simulation, resnum will not be allocated.
cccccccccccccc
ccccccccccccc      allocate( res_buf(n) )
ccccccccccccc      if (.not. allocated (resnum ) ) then
ccccccccccccc         allocate( resnum(n) )
ccccccccccccc      end if
ccccccccccccc      resnum = 0.0d0
ccccccccccccc
cccccccccccccc
cccccccccccccc     prepare the residue buffer
cccccccccccccc
ccccccccccccc      do iatom=1, n
ccccccccccccc          res_buf(iatom) = resnum(iatom)
ccccccccccccc      end do
cccccccccccccc
cccccccccccccc     send the residues
cccccccccccccc
ccccccccccccc      call MDI_Send(res_buf, n, MDI_INT, comm, ierr)
ccccccccccccc      if ( ierr .ne. 0 ) then
ccccccccccccc         write(iout,*)'SEND_RESIDUES -- MDI_Send failed'
ccccccccccccc         call fatal
ccccccccccccc      end if
ccccccccccccc      deallocate( res_buf )
      return
      end subroutine send_residues

c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_molecules  --  Respond to "<MOLECULES"     ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_molecules(comm)
        use iounit , only : iout
   1    use mdi , only    : MDI_INT, MDI_Send, MDI_Conversion_Factor
        use molcul , only  : molcule
        use atoms , only : n
        implicit none
        integer, intent(in)          :: comm
        integer                      :: ierr, iatom
        integer, allocatable         :: mol_buf(:)

c
c     prepare the residue buffer
c
      allocate( mol_buf(n) )
      do iatom=1, n
          mol_buf(iatom) = molcule(iatom)
      end do
c
c     send the residues
c
      call MDI_Send(mol_buf, n, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_MOLECULES -- MDI_Send failed'
         call fatal
      end if
      deallocate( mol_buf )
      return
      end subroutine send_molecules


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_pole_indices  --  Respond to "<IPOLES"     ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_pole_indices(comm)
        use iounit , only : iout
   1    use mdi , only    : MDI_INT, MDI_Send, MDI_Conversion_Factor
        use mpole , only  : ipole, npole
        use atoms , only : n
        implicit none
        integer, intent(in)          :: comm
        integer                      :: ierr, polei, pole_ind
        integer, allocatable         :: pole_buf(:)

        allocate( pole_buf(n) )
        pole_buf = 0
c
c     prepare the pole index buffer
c
      do polei=1, npole
          pole_ind = ipole(polei)
          pole_buf(pole_ind) = polei
      end do
c
c     send the pole indices
c
      call MDI_Send(pole_buf, n, MDI_INT, comm, ierr)
      if ( ierr .ne. 0 ) then
         write(iout,*)'SEND_IPOLE -- MDI_Send failed'
         call fatal
      end if
      deallocate( pole_buf )
      return
      end subroutine send_pole_indices


c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_field  --  Respond to "<FIELD"             ##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_field(comm)
ccccccccccccc      use atoms , only   : n
ccccccccccccc      use charge , only  : nion, iion, pchg
ccccccccccccc      use efield , only  : fielde
ccccccccccccc      use iounit , only  : iout
ccccccccccccc1     use mdi , only     : MDI_DOUBLE, MDI_Send
ccccccccccccc      use mpole , only   : npole
ccccccccccccc      use uprior , only  : use_pred
ccccccccccccc
ccccccccccccc      implicit none
      integer, intent(in)          :: comm
ccccccccccccc      integer                      :: ierr, ipole, dim
ccccccccccccc      real*8, allocatable          :: field(:)
ccccccccccccc      real*8                       :: epot
ccccccccccccc      real*8, allocatable          :: derivs(:,:)
ccccccccccccc      logical                      :: use_pred_original
cccccccccccccc
cccccccccccccc     the @DEFAULT node must calculate the latest UFIELD
cccccccccccccc
ccccccccccccc      allocate( field(3*npole) )
ccccccccccccc      if ( current_node .eq. "@DEFAULT" .and. forces_need_update ) then
cccccccccccccc
cccccccccccccc     turn off prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred_original = use_pred
ccccccccccccc         use_pred = .false.
cccccccccccccc
cccccccccccccc     allocate array to hold gradients
cccccccccccccc
ccccccccccccc         allocate( derivs(3,n) )
cccccccccccccc
cccccccccccccc     calculate the gradients
cccccccccccccc
ccccccccccccc         mdi_ignore_nodes = .true.
ccccccccccccc         call gradient (epot,derivs)
ccccccccccccc         mdi_ignore_nodes = .false.
ccccccccccccc         forces_need_update = .false.
cccccccccccccc
cccccccccccccc     reset prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred = use_pred_original
cccccccccccccc
cccccccccccccc     deallocate the gradients
cccccccccccccc
ccccccccccccc         deallocate( derivs )
ccccccccccccc      end if
cccccccccccccc
cccccccccccccc     construct the field array
cccccccccccccc
ccccccccccccc      do ipole=1, npole
ccccccccccccc         do dim=1, 3
ccccccccccccc            field(3*(ipole-1) + dim) = fielde(dim, ipole)
ccccccccccccc         end do
ccccccccccccc      end do
cccccccccccccc
cccccccccccccc     send the field
cccccccccccccc
ccccccccccccc      call MDI_Send(field, 3*npole, MDI_DOUBLE, comm, ierr)
ccccccccccccc      if ( ierr .ne. 0 ) then
ccccccccccccc         write(iout,*)'SEND_FIELD -- MDI_Send failed'
ccccccccccccc         call fatal
ccccccccccccc      end if
ccccccccccccc      deallocate( field )
      return
      end subroutine send_field
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_dfield_components  --  Respond to "<DFIELD"##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_dfield_components(comm)
ccccccccccccc      use atoms , only   : n
ccccccccccccc      use efield , only  : nprobes, dfield_pair, fielde
ccccccccccccc      use iounit , only  : iout
ccccccccccccc1     use mdi , only     : MDI_DOUBLE, MDI_Send
ccccccccccccc      use mpole , only   : npole
ccccccccccccc      use uprior , only  : use_pred
ccccccccccccc
ccccccccccccc      implicit none
      integer, intent(in)          :: comm
ccccccccccccc      integer                      :: ierr, i, j, dim
ccccccccccccc      real*8, allocatable          :: field(:)
ccccccccccccc      real*8                       :: epot
ccccccccccccc      real*8, allocatable          :: derivs(:,:)
ccccccccccccc      logical                      :: use_pred_original
cccccccccccccc
cccccccccccccc     the @DEFAULT node must calculate the latest DFIELD
cccccccccccccc
ccccccccccccc      allocate( field(3*nprobes*npole) )
ccccccccccccc      if ( current_node .eq. "@DEFAULT" .and. forces_need_update ) then
cccccccccccccc
cccccccccccccc     turn off prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred_original = use_pred
ccccccccccccc         use_pred = .false.
cccccccccccccc
cccccccccccccc     allocate array to hold gradients
cccccccccccccc
ccccccccccccc         allocate( derivs(3,n) )
cccccccccccccc
cccccccccccccc     calculate the gradients
cccccccccccccc
ccccccccccccc         mdi_ignore_nodes = .true.
ccccccccccccc         call gradient (epot,derivs)
ccccccccccccc         forces_need_update = .false.
ccccccccccccc         mdi_ignore_nodes = .false.
cccccccccccccc
cccccccccccccc     reset prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred = use_pred_original
cccccccccccccc
cccccccccccccc     deallocate the gradients
cccccccccccccc
ccccccccccccc         deallocate( derivs )
ccccccccccccc      end if
cccccccccccccc
cccccccccccccc     construct the field array
cccccccccccccc
ccccccccccccc      do i=1, nprobes
ccccccccccccc        do j=1, npole
ccccccccccccc          do dim=1, 3
ccccccccccccc             field(3*npole*(i-1)+3*(j-1)+dim) = dfield_pair(dim, j, i)
ccccccccccccc          end do
ccccccccccccc        end do
ccccccccccccc      end do
cccccccccccccc
cccccccccccccc     send the field
cccccccccccccc
ccccccccccccc      call MDI_Send(field, 3*npole*nprobes, MDI_DOUBLE, comm, ierr)
ccccccccccccc      if ( ierr .ne. 0 ) then
ccccccccccccc         write(iout,*)'SEND_CHARGES -- MDI_Send failed'
ccccccccccccc         call fatal
ccccccccccccc      end if
ccccccccccccc      deallocate( field )
      return
      end subroutine send_dfield_components
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine send_ufield_components  --  Respond to "<UFIELD"##
c     ##                                                             ##
c     #################################################################
c
      subroutine send_ufield_components(comm)
ccccccccccccc      use atoms , only   : n
ccccccccccccc      use efield , only  : nprobes, ufield_pair, fielde
ccccccccccccc      use iounit , only  : iout
ccccccccccccc1     use mdi , only     : MDI_DOUBLE, MDI_Send
ccccccccccccc      use mpole , only   : npole
ccccccccccccc      use uprior , only  : use_pred
ccccccccccccc
ccccccccccccc      implicit none
      integer, intent(in)          :: comm
ccccccccccccc      integer                      :: ierr, i, j, dim
ccccccccccccc      real*8, allocatable          :: field(:)
ccccccccccccc      real*8                       :: epot
ccccccccccccc      real*8, allocatable          :: derivs(:,:)
ccccccccccccc      logical                      :: use_pred_original
cccccccccccccc
cccccccccccccc     the @DEFAULT node must calculate the latest UFIELD
cccccccccccccc
ccccccccccccc      allocate( field(3*nprobes*npole) )
ccccccccccccc      if ( current_node .eq. "@DEFAULT" .and. forces_need_update ) then
cccccccccccccc
cccccccccccccc     turn off prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred_original = use_pred
ccccccccccccc         use_pred = .false.
cccccccccccccc
cccccccccccccc     allocate array to hold gradients
cccccccccccccc
ccccccccccccc         allocate( derivs(3,n) )
cccccccccccccc
cccccccccccccc     calculate the gradients
cccccccccccccc
ccccccccccccc         mdi_ignore_nodes = .true.
ccccccccccccc         call gradient (epot,derivs)
ccccccccccccc         forces_need_update = .false.
ccccccccccccc         mdi_ignore_nodes = .false.
cccccccccccccc
cccccccccccccc     reset prediction of induced dipoles
cccccccccccccc
ccccccccccccc         use_pred = use_pred_original
cccccccccccccc
cccccccccccccc     deallocate the gradients
cccccccccccccc
ccccccccccccc         deallocate( derivs )
ccccccccccccc      end if
cccccccccccccc
cccccccccccccc     construct the field array
cccccccccccccc
ccccccccccccc      do i=1, nprobes
ccccccccccccc        do j=1, npole
ccccccccccccc          do dim=1, 3
ccccccccccccc             field(3*npole*(i-1)+3*(j-1)+dim) = ufield_pair(dim, j, i)
ccccccccccccc          end do
ccccccccccccc        end do
ccccccccccccc      end do
cccccccccccccc
cccccccccccccc     send the field
cccccccccccccc
ccccccccccccc      call MDI_Send(field, 3*npole*nprobes, MDI_DOUBLE, comm, ierr)
ccccccccccccc      if ( ierr .ne. 0 ) then
ccccccccccccc         write(iout,*)'SEND_CHARGES -- MDI_Send failed'
ccccccccccccc         call fatal
ccccccccccccc      end if
ccccccccccccc      deallocate( field )
      return
      end subroutine send_ufield_components

      end module mdiengine
