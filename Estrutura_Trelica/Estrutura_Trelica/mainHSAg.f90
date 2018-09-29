program HSAg_teste
    use HSAg
    use funcobj
    use ocl
    use HSAr
    
    implicit none
    
    !---------------------------------------------
    !variáveis do programa
    !---------------------------------------------

  
    !---Variáveis usadas na execução do HSA
    type(HSA_par):: param
    type(HSA_sai),allocatable:: sHSA(:)
    integer:: ex   !Contador de execuções do HSA
    integer:: NEA  !Número de execuções do HSA
    real(8):: ffop !tolerância para definir ótimo prático: (fob-fot)/fot < ffop
    character(len=500):: tipo_problema !definição do tipo de problema a ser executado
    real(8),allocatable:: Rot(:), Rop(:)    !Confiabilidade aparente para o ótimo e ótimo prático, resp.
    integer,allocatable:: Eot(:), Eop(:)    !Núm. de avaliações da função objetivo relativo a Rot e Rop, resp.
    
    !---Variáveis auxiliares
    integer:: i, ii, iii !contadores
    integer:: ncla
    character(len=500):: nomearquivo, texto
    
    !---Variáveis definindo as unidades de entrada e saída de dados    
    integer,parameter:: uni_dad = 501, uni_print=502, uni_parametros = 503, uni_sai = 504, uni_Rot = 505, uni_Rop = 506
    
    
    
    !---------------------------------------------
    !definição das namelist para entrada dedados
    !---------------------------------------------
    namelist /parametros_HSA/ param
    namelist /execHSA/ tipo_problema, NEA, ffop
    
    

    !---------------------------------------------
    !Início das operações
    !---------------------------------------------

    
    !-----------------
    !Obtendo nome do arquivo (argumentos ou entrada no terminal)
    ncla = COMMAND_ARGUMENT_COUNT ()
    if(ncla/=0) then
       call GET_COMMAND_ARGUMENT(1,nomearquivo,ii,iii)
    else
        write(*,'(A)',advance='no')"Nome dos arquivos:"
        read(*,*) nomearquivo
    endif
    !-----------------
    
       
    !-----------------
    !abrindo arquivo e lendo os parâmetros do HSA
    open( unit=uni_dad , file=trim(nomearquivo)//".dad" , action='read',iostat=iii)
    if(iii/=0) then
        write(*,*)"ERR>> problema na abertura do arquivo:" , trim(nomearquivo)//".dad"
    else
        write(*,*)"arquibo aberto: ", trim(nomearquivo)//".dad"
    endif

    rewind(uni_dad)
    read(uni_dad,nml=execHSA)
    if(NEA<1) stop "ERR>> NEA<1"
    
    rewind(uni_dad)
    read(uni_dad,nml=parametros_HSA)
    
    !close(uni_dad)
    !-----------------   
    
    
    !-----------------
    !abrindo arquivos de saída de informações e escrevendo parâmetros do HSA usados na otimização
    open( unit=uni_print , file=trim(nomearquivo)//"_saida_geral.txt" , action='WRITE',iostat=iii)
    open( unit=uni_parametros , file=trim(nomearquivo)//"_parHSA.txt" , action='WRITE',iostat=iii)
    param%uni_imp_par = uni_parametros
    write(uni_print,nml=parametros_HSA)
    open( unit=uni_sai , file=trim(nomearquivo)//".sHS" , action='WRITE',iostat=iii)
    open(unit=uni_Rot,file=trim(nomearquivo)//"_R.plt", action='WRITE',iostat=iii)
    open(unit=uni_Rop,file=trim(nomearquivo)//"_Rop.plt", action='WRITE',iostat=iii)
    !-----------------
    
    !-----------------
    !Alocando as variáveis
    allocate( sHSA(NEA), Rot(NEA), Eot(NEA), Rop(NEA), Eop(NEA))
    !-----------------
    
    
    !-----------------
    !Executando as otimizações
    
    selectcase( trim(tipo_problema) )
    case('ocl1')    
        
        !Iniciando as informações da rotina de avaliação da função objetivo para o problema ocl1
        call iniciar_ocl1(uni_dad)

        do ex = 1 , NEA
            call iniciar_HSA_sai(param,sHSA(ex))
            call otimiza_HSA(f_ocl1,param,sHSA(ex))
            call impr_ocl1(param%nVD,sHSA(ex)%ot%vd,texto)
            write(*,'(i3.3," fob=",ES14.6,2x,A)')ex,sHSA(ex)%ot%fob,trim(texto(param%nVD*2:))
            write(uni_sai,'(A,i2.2," fob= ",ES14.6,2x,A)') "ocl1 -> ex=",ex,sHSA(ex)%ot%fob,trim(texto)
        enddo            
        
        call avaliar_HSAr(NEA,ffop,param,sHSA,Eot,Rot,Eop,Rop)
        call imprimir_HSAr(uni_Rot,Eot,Rot,uni_Rop,Eop,Rop)    
        

        if(NEA==1) then
               open(unit=200,file=trim(nomearquivo)//"_hAF.plt", action='WRITE',iostat=iii)
               call imprimir_HSA_sai_hAF_plt(param,sHSA(1),200)
        endif
        
            
            !do i = 1 , param%maxHMS
            !    call impr_ocl1(param%nVD,sHSA%HM(i)%vd,texto)
            !    write(*,'(A,i2.2,2x,A," fob= ",ES14.6)') "ocl1 -> HM=",i,trim(texto),sHSA%HM(i)%fob
            !    write(uni_sai,'(A,i2,2x,A," fob= ",ES14.6)') "ocl1 -> HM=",i,trim(texto),sHSA%HM(i)%fob
            !enddo
            !write(uni_sai,*)"nAF_ot=",sHSA%naf_ot
        
    case default
        
        call iniciar_HSA_sai(param,sHSA(1))
    
        call otimiza_HSA(fob1,param,sHSA(1))
        
    end select
    
    !write(uni_print,*) "------------------------- saida ---------------------------"
    
    
    !call imprimir_HSA_sai(param,sHSA,uni_sai)
    
    !write(uni_sai,*)"ot_fob=",sHSA(1)%ot%fob
    !write(uni_sai,*)"ot_vd=",sHSA%ot%vd
    !write(uni_sai,*)"ot_vc=",sHSA%ot%vc
    !
    
    !close(uni_print)
    close(uni_parametros)
    


    !-----------------
    
    
    
    
    
    end program HSAg_teste
    
    
    
      
