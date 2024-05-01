 program NOCalculator

    !   This program computes the MO difference between two matfiles and returns a
    !   matfile that is used to compute Natural Orbitals 

    !****x* Main/matDiffNOCalc
    !   A. M. Kinyua, 2023
    !*  NAME
    !*      Natural Orbitals Calculator 
        
        use mqc_gaussian
        use iso_fortran_env, only: int32, int64, real64

    !
    ! Variables ----------------------------------------------------
        implicit none 
        type(mqc_gaussian_unformatted_matrix_file)::fileInfo
        character(len=:),allocatable::command,fileName,help_path,density_difference
        character(len=256),dimension(:),allocatable::fileList
        integer(kind=int64)::iOut=6,iPrint=1,numFile=0,i,j,k,l,stat_num,iUnit,nBasis,nElec
        type(mqc_matrix),dimension(:),allocatable::mo_list
        type(mqc_scf_integral)::temp_int1,temp_int2,mo_diff,overlap
        type(mqc_gaussian_unformatted_matrix_file)::temp_file
        type(mqc_wavefunction)::common_wave
        type(mqc_matrix)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
        
    ! USAGE --------------------------------------------------------
    !*      NOCalc.exe [-f <input_file_with_matrix_lists>] 
    !*
    !*    OPTIONS
    !*
    !
    !     Print program information.
    !
           write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
               '  Natural Orbitals Calculator',NEW_LINE('a'), &
               ' ',repeat('*',73),NEW_LINE('a'), &
               NEW_LINE('a'),repeat(' ',30),'Version 23.1.0',NEW_LINE('a'),NEW_LINE('a'),&
               'L. M. Thompson &  A. M. Kinyua, Louisville KY, 2022.',NEW_LINE('a')


    !
    !     Parse input Options ---------------------------------------
    !   

        j=1
        do i=1,command_argument_count()
            if(i.ne.j) cycle
            call mqc_get_command_argument(i,command)
           if(command.eq.'-f') then
    !
    !*      -f matrix_file                   Input matrix file with initial set of molecular orbitals.
    !*                                       The first line contains the number of matrix files in the
    !*                                       input, and then on each line is a separate matrix file.
    !*
               call mqc_get_command_argument(i+1,fileName)
               j = i+2
            elseif(command.eq.'--densitydifference') then
               call mqc_get_command_argument(i+1,command)
               density_difference = command
               j = i+2
            endIf
        endDo


    !
    !*  Parse input file ---------------------------------------------
    !
        if(.not.allocated(fileName)) call mqc_error('No input file provided', iOut)
        open(newunit=iUnit, file=fileName,status='old',iostat=stat_num)
        if(stat_num/=0) call mqc_error_a('Error opening file',iOut,'fileName',fileName)
        read(unit=iUnit,fmt='(i20)',iostat=stat_num) numFile
        if(stat_num/=0) call mqc_error('Error reading file number',iOut)
        allocate(fileList(numFile))
        do i = 1, numFile
            read(unit=iUnit,fmt='(A)',iostat=stat_num) fileList(i)
            if((stat_num<0).and.(i<=numFile)) call mqc_error('File EOF reached early',iOut)
        endDo
        close(unit=iUnit)
    !
    !******
    !  RUN NATORBITALS FOR VERY INSTANCE; OTHER THINGS TOO, IN A LOOP
    !

            call temp_file%getArray('SHELL TO ATOM MAP',sh2AtMp,fileName=fileList(1))
            call temp_file%getArray('SHELL TYPES',shlTyp,fileName=fileList(1))
            call temp_file%getArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh,fileName=fileList(1))
            call temp_file%getArray('PRIMITIVE EXPONENTS',prmExp,fileName=fileList(1))
            call temp_file%getArray('CONTRACTION COEFFICIENTS',conCoef,fileName=fileList(1))
            call temp_file%getArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo,fileName=fileList(1))
            call temp_file%getArray('COORDINATES OF EACH SHELL',shCoor,fileName=fileList(1))
    ! 
            call temp_file%getESTObj('overlap',est_integral=overlap,fileName=fileList(1))
            call temp_file%getESTObj('wavefunction',common_wave,fileName=fileList(1))
            nBasis = common_wave%nbasis%ival()
            nElec = common_wave%nelectrons%ival()
        !allocate(mo_list(numFile))
        !do i=1, numFile
            call temp_file%getESTObj('density',est_integral=temp_int1,fileName=fileList(1))
        if(.not.allocated(mo_list)) allocate(mo_list(1))
        do i=1,numFile
          !  if(allocated(mo_list(1))) deallocate(mo_list)
            call temp_file%getESTObj('density',est_integral=temp_int2,fileName=fileList(i))
            if(density_difference.eq.'yes') then
                mo_diff = temp_int1-temp_int2
                mo_list(1) = mo_diff%getBlock('full')
                call mo_list(1)%print(6,'MO density DIFFERENCE at step '//trim(num2char(i)))
            elseif(density_difference.eq.'no') then
                mo_diff = temp_int2
                mo_list(1) = mo_diff%getBlock('full')
                call mo_list(1)%print(6,'MO density at step '//trim(num2char(i)))
            endIf
! computing natural orbitals 
       call naturalOrbital(nBasis,nElec,mo_list,overlap,sortedNatValsOut=common_wave%mo_energies,&
             Density_in_in=.true.)
      call mo_list(1)%print(6,'Natural orbitals at step '//trim(num2char(i)))
      call common_wave%mo_energies%print(6,'Occupation numbers at step '//trim(num2char(i)))
! writing output matrix file
      call writeOutputFiles(temp_file,'NatOrbs-wholedensitymatrix-'//trim(num2char(i))//'-',mo_list,nBasis,&
             1.0e-14,sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor,common_wave)
   endDo
!  
! 
! 
 contains 

   subroutine naturalOrbital(nBasis,nElec,mo_list,overlapMatrix,activeSpace, &
       sortedNatValsOut,occ_thresh_in,virt_thresh_in,NOAvg_in,Density_in_in)
!
! Build the natural orbitals and overwrite incoming MOs. Optionally, return string 
! containing active space information based on fractional occupation of each set of MOs.
! The threshold for specifying upper and lower limits can be specified using the optional
! arguments occ_thresh_in and virt_thresh_in. If NOAvg set, return natural orbitals of 
! average density from all input MOs and return only a single set of orbitals. If 
! Density_in set, mo_list contains densities and not MOs.
!
   implicit none
   integer(kind=int64),intent(in)::nBasis,nElec
   type(mqc_matrix),dimension(:),allocatable,intent(inOut)::mo_list
   type(mqc_scf_integral),intent(in)::overlapMatrix
   character(len=*),optional,intent(inOut)::activeSpace
   type(mqc_scf_eigenvalues),optional,intent(inOut)::sortedNatValsOut
   real(kind=real64),optional,intent(in)::occ_thresh_in,virt_thresh_in
   logical,optional,intent(in)::NOAvg_in,Density_in_in

   ! Processing variables
   type(mqc_matrix)::densityMatrix,totalDensity,natOrbs,shalf,sminushalf,&
     tmpMO,IDmatrix
   type(mqc_vector)::natVals,sortedNatVals
   integer(kind=int64)::numFile,i,j,counter,actOrbs,inactOrbs,actElecs,&
     nAlpha,nBeta,nBasUse
   real(kind=real64)::occ_thresh,virt_thresh
   logical::NOAvg,Density_in


   nAlpha = nElec/2 
   nBeta = nElec/2 
   numFile = size(mo_list)
!*****
   if(present(occ_thresh_in)) then
     occ_thresh = occ_thresh_in
   else
     occ_thresh = 1.98
   endIf
   if(present(virt_thresh_in)) then
     virt_thresh = virt_thresh_in
   else
     virt_thresh = 0.02
   endIf
   if(present(NOAvg_in)) then
     NOAvg = NOAvg_in
   else
     NOAvg = .false.
   endIf
   if(present(Density_in_in)) then
     Density_in = Density_in_in
   else
     Density_in = .false.
   endIf
   shalf = overlapMatrix%getBlock('alpha-alpha')
   call shalf%power(0.5)
   sminushalf = overlapMatrix%getBlock('alpha-alpha')
   call sminushalf%power(-0.5)
   activeSpace = '['
   call densityMatrix%init(nBasis*2,nBasis*2)
   do i = 1, numFile
     if(NOAvg) then
       if(density_in) then
         densityMatrix = densityMatrix + mo_list(i)
       else
         densityMatrix = densityMatrix + matmul(moBuilder(mo_list(i),nAlpha,nBeta,nBasis), &
           dagger(moBuilder(mo_list(i),nAlpha,nBeta,nBasis)))
       endIf
       if(i.ne.numFile) then
         cycle
       else
         densityMatrix = densityMatrix/numFile
       endIf
     else
       if(density_in) then
         densityMatrix = mo_list(i)
       else
         densityMatrix = matmul(moBuilder(mo_list(i),nAlpha,nBeta,nBasis), & 
           dagger(moBuilder(mo_list(i),nAlpha,nBeta,nBasis)))
       endIf
     endIf
     totalDensity = densityMatrix%mat([1,nBasis],[1,nBasis]) + &
       densityMatrix%mat([nBasis+1,nBasis*2],[nBasis+1,nBasis*2]) + &
       densityMatrix%mat([nBasis+1,nBasis*2],[1,nBasis]) + &
       densityMatrix%mat([1,nBasis],[nBasis+1,nBasis*2])
     totalDensity = matmul(matmul(shalf,totalDensity),shalf)
     call totalDensity%diag(natVals,natOrbs)
     natOrbs = matmul(sminushalf,natOrbs)
     call natvals%print(6,'LMTLMT natvals')

    !   Natural values and corresponding orbitals are returned from eigensystem
    !   with lowest value first.  Desired order is to have the first orbitals
    !   correspond to inactive orbitals.

     call mo_list(i)%init(2*nBasis,2*nBasis)
     call sortedNatVals%init(size(natVals))

! *****
!    call sortedNatVals%print(6,'sorted natural orbitals ')
! *****
     counter = 1
     do j = nBasis, 1, -1
       call mo_list(i)%mput(natOrbs%mat([1,nBasis],[j,j]),[1,nBasis],[counter,counter])
       call mo_list(i)%mput(natOrbs%mat([1,nBasis],[j,j]),[nBasis+1,nBasis*2],[counter+nBasis,counter+nBasis])
       call sortedNatVals%put(natVals%at(j),counter)
       counter = counter + 1
     end do
! turn the below on to pick between real or imaginary natural orbitals
     !mo_list(i) = real(mo_list(i))
     !mo_list(i) = aimag(mo_list(i))
     IDmatrix = matmul(matmul(dagger(mo_list(i)),overlapMatrix),mo_list(i))
     call IDmatrix%print(6,'Identity matrix '//trim(num2char(i)))
     call sortedNatVals%print(6,'sorted natural orbitals ')
     !flush(6)
     if(present(sortedNatValsOut)) then
          call sortedNatVals%print(6,'LMTLMT ONs before')
         call mqc_eigenvalues_allocate(sortedNatValsOut,'mo energies','general',sortedNatVals/2.0,sortedNatVals/2.0)
     endIf
     !
     inactOrbs = 0
     actOrbs = 0
     actElecs = 0
     do j = 1, nBasis
         if(natVals%at(j).ge.occ_thresh) then
             inactOrbs = inactOrbs + 1
         elseif((natVals%at(j).lt.occ_thresh).and.(natVals%at(j).gt.virt_thresh)) then
             actOrbs = actOrbs + 1
         endif
     end do

     actElecs = nElec - (inactOrbs * 2)

     activeSpace = trim(activeSpace)//trim(num2char(actOrbs))
     activeSpace = trim(activeSpace)//','// trim(num2char(actElecs))

     if(i.lt.numFile) then
         activeSpace = trim(activeSpace)//':'
     else
         activeSpace = trim(activeSpace)//']'
     endif

    enddo

    if(NOAvg) then
        tmpMO = mo_list(numFile)
        deallocate(mo_list) 
        allocate(mo_list(1))
        mo_list(1) = tmpMO
    endIf
end subroutine
!
!
    
   function moBuilder(moMatrix,nAElec,nBElec,nBasFunc,occAlpha,occBeta)
!
! Return occupied orbitals from the set of spinor-basis MOs. Either the
! first n orbitals can be extracted, or a specific orbital, of each spin.
!
   implicit none
   type(mqc_matrix)::moBuilder,moMatrix
   integer(kind=int64)::nAElec,nBElec,nElectrons,nBasFunc
   integer(kind=int64),dimension(:),optional::occAlpha,occBeta


   nElectrons = nAElec + nBElec
   call moBuilder%init(nBasFunc*2,nElectrons)
 
 !****
    call moBuilder%print(6,'Initialized moBuildr')
    call moMatrix%print(6,'MO Matrix 1')
 !****

   if(present(occAlpha)) then
     call moBuilder%mput(moMatrix%mat([1,nBasFunc*2],occAlpha),[1,nBasFunc*2],[1,nAElec])
   else  
     call moBuilder%mput(moMatrix%mat([1,nBasFunc*2],[1,nAElec]),[1,nBasFunc*2],[1,nAElec])
   endif

  !****
   call moBuilder%print(6,'Second moBuildr')
   call moMatrix%print(6,'MO Matrix 2')
  !****

   if(nBElec.gt.0) then
     if(present(occBeta)) then
       call moBuilder%mput(moMatrix%mat([1,nBasFunc*2],occBeta), & 
         [1,nBasFunc*2],[nAElec+1,nElectrons])
     else
       call moBuilder%mput(moMatrix%mat([1,nBasFunc*2],[nBasFunc+1,nBasFunc+nBElec]), &
         [1,nBasFunc*2],[nAElec+1,nElectrons])
     endif
   endif
  return

 end function moBuilder

subroutine writeOutputFiles(fileInfo,outFileName,mo_list,nBasis,THRESH_ZERO,sh2AtMp,shlTyp,&
     nPrmSh,prmExp,conCoef,conCoTwo,shCoor,wavefunction,ERI)

 implicit none

 ! input/output variables
  type(mqc_gaussian_unformatted_matrix_file),intent(inOut)::fileInfo
  character(len=*),intent(in)::outFileName
  type(mqc_matrix),dimension(:),allocatable,intent(in)::mo_list
  integer(kind=int64),intent(in)::nBasis
  real(kind=real64),intent(in)::THRESH_ZERO
  type(mqc_matrix),intent(in)::sh2AtMp,shlTyp,nPrmSh,prmExp,conCoef,conCoTwo,shCoor
  type(mqc_wavefunction),intent(in)::wavefunction
  type(mqc_twoERIs),dimension(:),optional,allocatable,intent(in)::ERI

 ! scratch variables
  integer::i
  character(len=80)::newFileName,spinSymStr
  type(mqc_scf_integral)::outMat
  type(mqc_matrix)::amat,bMat,abMat,baMat
  type(mqc_scf_integral)::density

   do i = 1, size(mo_list)
     newFileName = outFileName
     call build_string_add_int(i,newFileName,20)
     newFileName = trim(newFileName) // '.mat'
     call outMat%init(nBasis*2,nBasis*2)
     aMat = mo_list(i)%mat([1,nBasis],[1,nBasis])
     bMat = mo_list(i)%mat([nBasis+1,nBasis*2],[nBasis+1,nBasis*2])
     abMat = mo_list(i)%mat([nBasis+1,nBasis*2],[1,nBasis])
     baMat = mo_list(i)%mat([1,nBasis],[nBasis+1,nBasis*2])
     fileInfo%icgu = 111
     if(MQC_Matrix_Norm((aMat-bMat)).lt.THRESH_ZERO.and.abMat%norm().lt.THRESH_ZERO.and. &
       baMat%norm().lt.THRESH_ZERO) then
       spinSymStr = 'space'
       call mqc_integral_allocate(outMat,'',spinSymStr,aMat)
       if(MQC_Matrix_HaveComplex(mo_list(i))) fileInfo%icgu = fileInfo%icgu + 10
     elseIf(abMat%norm().lt.THRESH_ZERO.and.baMat%norm().lt.THRESH_ZERO) then
       spinSymStr = 'spin'
       call mqc_integral_allocate(outMat,'',spinSymStr,aMat,bMat)
       fileInfo%icgu = fileinfo%icgu + 1
       if(MQC_Matrix_HaveComplex(mo_list(i))) fileInfo%icgu = fileInfo%icgu + 10
     else
       spinSymStr = 'general'
       call mqc_integral_allocate(outMat,'',spinSymStr,aMat,bMat,abMat,baMat)
       fileInfo%icgu = fileInfo%icgu + 100
       fileInfo%icgu = fileInfo%icgu + 10
     endIf
     
     call mo_list(i)%print(6,'MO density '//trim(num2char(i)))
    
     if(fileInfo%isOpen()) call Close_MatF(fileInfo%UnitNumber)

     call fileInfo%create(newFileName)
     call fileInfo%writeArray('SHELL TO ATOM MAP',sh2AtMp)
     call fileInfo%writeArray('SHELL TYPES',shlTyp)
     call fileInfo%writeArray('NUMBER OF PRIMITIVES PER SHELL',nPrmSh)
     call fileInfo%writeArray('PRIMITIVE EXPONENTS',prmExp)
     call fileInfo%writeArray('CONTRACTION COEFFICIENTS',conCoef)
     call fileInfo%writeArray('P(S=P) CONTRACTION COEFFICIENTS',conCoTwo)
     call fileInfo%writeArray('COORDINATES OF EACH SHELL',shCoor)
     call fileInfo%writeESTObj('mo coefficients',est_integral=outMat,override=spinSymStr)
     !call fileInfo%writeESTObj('density',est_integral=outMat,override=spinSymStr)
     call wavefunction%mo_energies%print(6,'LMTLMT ONs')
     call fileInfo%writeESTObj('mo energies',est_eigenvalues=wavefunction%mo_energies,override=spinSymStr)
     call fileInfo%writeESTObj('overlap',est_integral=wavefunction%overlap_matrix,override=spinSymStr)
     call fileInfo%writeESTObj('core hamiltonian',est_integral=wavefunction%core_hamiltonian,override=spinSymStr)
     if(i.eq.1.and.present(ERI)) then
       if(ERI(1)%type().eq.'regular') then
         call fileInfo%write2ERIs('regular',ERI)
       elseIf(ERI(1)%type().eq.'raffenetti1') then
         call fileInfo%write2ERIs('raffenetti',ERI)
       else
         call mqc_error_a('Output two electron integrals need to be either regular or raffenetti type',6,&
           'ERI(1)%type()',ERI(1)%type())
       endIf
     endIf
     call Close_MatF(fileInfo%UnitNumber)
    write(6,'(A,A)') 'Writing matrix file ',newFileName
   enddo
  return
 end subroutine
!
end program NOCalculator






























