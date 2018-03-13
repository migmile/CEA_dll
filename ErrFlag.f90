
module CommonData
	character(len=256) ::cmData_dir_name
    character(len=256) ::cmLog_dir_name
    integer(kind=4) :: ErrorFlag
    
    contains
    subroutine WrtError(text)
    character(len=*) text
      open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
      if (ErrorFlag/=0) then 
           write(20,*)"___________________"
           write(20,*) ErrorFlag
      endif 
      write(20,*) text

      close(20)
    end subroutine WrtError
    
end Module CommonData
