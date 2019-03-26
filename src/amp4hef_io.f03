module amp4hef_io
  implicit none

  character(128),protected :: path_tbldir=&
  !(path_tbldir!)

contains


  function get_path( fileName ,option ) result(rslt)
  character(*),intent(in) :: fileName
  character(*),intent(in),optional :: option
  character(144) :: rslt,optVal
  optVal = 'file'
  if (present(option)) optVal = option
  select case (fileName)
  case ('amp4hef.tbl')
    select case (trim(optVal))
    case ('dir'     ) ;rslt = trim(path_tbldir)
    case ('file'    ) ;rslt = trim(path_tbldir)//'amp4hef.tbl'
    case ('pathName') ;rslt = 'path_tbldir'
    end select
  end select
  end function

  subroutine check_file( fileName )
  character(*),intent(in) :: fileName
  logical :: file_exists
  inquire( file=trim(get_path(fileName)), exist=file_exists )
  if (file_exists) return
  write(*,'(A)') "ERROR in amp4hef:"
  write(*,'(A)') "  the file "//trim(get_path(fileName))
  write(*,'(A)') "  does not seem to exist. You can set the path to the directory"
  write(*,'(A)') "  where this file is with"
  write(*,'(A)') "    call set_path('"//trim(get_path(fileName,'pathName'))//"',yourPath)"
  write(*,'(A)') "  where for example  yourPath='/home/user/project/files/'"
  write(*,'(A)') "  It must include the final delimiter (the final slash in this example)!"
  write(*,*)
  stop
  end subroutine

  subroutine set_path(pathName,pathVal)
  character(*),intent(in) :: pathName,pathVal
  select case(pathName)
  case ('path_tbldir') ;path_tbldir = pathVal
  end select
  end subroutine


  subroutine split_line( Nwords,iBgn,iEnd ,line )
  character(*),intent(in) :: line
  integer,intent(out) :: Nwords,iBgn(len(line)),iEnd(len(line))
  integer :: ii,charLen
  charLen = len(line)
  Nwords = 0
  ii = 0
  do
    do ;ii=ii+1 ;if(ii.gt.charLen)exit 
      if(line(ii:ii).ne.' ')exit
    enddo
    if (ii.gt.charLen) exit
    Nwords = Nwords+1
    iBgn(Nwords) = ii
    do ;ii=ii+1 ;if(ii.gt.charLen)exit 
      if(line(ii:ii).eq.' ')exit
    enddo
    iEnd(Nwords) = ii-1
    if (ii.gt.charLen) exit
  enddo
  end subroutine


end module


