 program Programa_Trelica

    


    use Analise_trelica_trapezoidal
    use HSAg
    
    implicit none
    
    integer :: nd                                      ! numero de variaveis diescretas
    integer :: nc                                      ! numero de variaveis continuas
    integer, allocatable :: Vd(:)                      ! valor das variaveis discretas
    real(8), allocatable :: Vc(:)                      ! valor das variaveis continuas
    real(8) :: fob                                     ! função objetivo
    character(14) :: TCPLOT = "VARIABLES ="
    
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
    integer:: ii, iii !contadores
    integer:: ncla
    character(len=500):: nomearquivo, texto
    
    !---Variáveis definindo as unidades de entrada e saída de dados    
    integer,parameter:: uni_dad = 501, uni_print=502, uni_parametros = 503, uni_sai = 504, uni_Rot = 505, uni_Rop = 506, uni_result = 507 
    
    
    
    !---------------------------------------------
    !definição das namelist para entrada dedados
    !---------------------------------------------
    namelist /parametros_HSA/ param
    namelist /execHSA/ tipo_problema, NEA, ffop
    ! -----------------------------------------------------------------------------------------------------------------------------
  
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
    open( unit=uni_dad , file=trim(nomearquivo) , action='read',iostat=iii)
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
    case('trelica1')    
        
        !-------------DADOS DE ENTRADA - TRELIÇA TRAPEZOIDAL---------------!
    !Abertura de arquivo e leitura de dados---------------------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S1)
    close (uni_dad_vento)
    
    ! Chamada para cálculo do fator topográfico S1 --------------------------------------------------------------------------------
    call Fator_Topografico_S1 (caso, theta, z, d, S1)
      
    !Leitura dos dados do problema referentes ao fator S2 -------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S2)
    close (uni_dad_vento)
    
    ! Chamada para cálculo do fator topográfico S2 --------------------------------------------------------------------------------
    call Transf_dados_S2 (tab_S2)
    call Fator_Topografico_S2(cat, cls, z, tab_S2, S2)
    
    !Leitura dos dados do problema referentes a pressão dinâmica ------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = pressao_dinamica)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes ao coeficiente de pressão -----------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Cofc_pressao)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes aos dados geométricos de entrada -----------------------------------------------------------
    open (unit = uni_dad_vento, file = 'Vento_Entrada.txt', action='read')        
        read (uni_dad_vento, nml = Dados_geometricos)
    close (unit = uni_dad_vento)
    
    !Armazenamento dos dados do problema referentes aos vários perfis de terças  -----------------------------------------------------------
    open (unit = uni_dad_tercas, file = 'dados_tercas.txt', action='read')        
        call lista_tercas(uni_dad_tercas, terca)
    close (unit = uni_dad_tercas)
    
    !Armazenamento das propriedades geométricas referentes aos vários perfis 2L da AISC  -----------------------------------------------------------    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do


        
        do ex = 1 , NEA
                 open( unit=uni_validacao , file="pp_validacao"//trim(str(ex))//".txt" , action='write')
                 write(uni_validacao, '(A, 2A7)') TCPLOT, 'nAF', 'pp'
            
            call iniciar_HSA_sai(param,sHSA(ex))
            call otimiza_HSA(funcao_objetivo,param,sHSA(ex))
            !call impr_ocl1(param%nVD,sHSA(ex)%ot%vd,texto)
            write(*,'(i3.3," fob=",ES14.6,2x,A)')ex,sHSA(ex)%ot%fob,trim(texto(param%nVD*2:))
            write(uni_sai,'(A,i2.2," fob= ",ES14.6,2x,A)') "ocl1 -> ex=",ex,sHSA(ex)%ot%fob,trim(texto)
            open( unit=uni_result , file="resultados_NEA_"//trim(str(ex))//".txt" , action='write',iostat=iii)
                write(uni_result, '(A, 7A7, 2A12)') TCPLOT, 'it', 'Vd(1)', 'Vd(2)', 'Vd(3)', 'Vd(4)', 'fob', 'pp', 'AF_ot','fob_ot'
                do i = 1, param%maxAF
                    write(uni_result, '(I20 , I6 , I8 , I6, I6, 2G15.6, I6, G15.6)') i, sHSA(ex)%HAF(i)%Vd(1:4), sHSA(ex)%HAF(i)%fob, peso_total(i), sHSA(ex)%OT_HAF(i), sHSA(ex)%HAF(sHSA(ex)%OT_HAF(i))%fob
                end do
            close(uni_result)
            close(uni_validacao)
            num_AF = 0.d0
            
            if(ex == 1) then
                inclinacao_diagonais = "/\"
            else if(ex == 2) then
                inclinacao_diagonais = "\/"
            end if
            
        enddo            
        
       ! call avaliar_HSAr(NEA,ffop,param,sHSA,Eot,Rot,Eop,Rop)
       ! call imprimir_HSAr(uni_Rot,Eot,Rot,uni_Rop,Eop,Rop)    
        

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
        
        !call iniciar_HSA_sai(param,sHSA(1))
    
        !call otimiza_HSA(fob1,param,sHSA(1))
        
    end select
    
    !write(uni_print,*) "------------------------- saida ---------------------------"
    
    
    !call imprimir_HSA_sai(param,sHSA,uni_sai)
    
    !write(uni_sai,*)"ot_fob=",sHSA(1)%ot%fob
    !write(uni_sai,*)"ot_vd=",sHSA%ot%vd
    !write(uni_sai,*)"ot_vc=",sHSA%ot%vc
    !
    
    !close(uni_print)
    close(uni_parametros)
    
    end program Programa_Trelica   