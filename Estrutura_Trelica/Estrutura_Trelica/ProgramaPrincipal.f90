 program Programa_Trelica

    


    use Analise_trelica_trapezoidal
    use HSAg
    
    implicit none
    
    integer :: nd                                      ! numero de variaveis diescretas
    integer :: nc                                      ! numero de variaveis continuas
    integer, allocatable :: Vd(:)                      ! valor das variaveis discretas
    real(8), allocatable :: Vc(:)                      ! valor das variaveis continuas
    real(8) :: fob                                     ! fun��o objetivo
    real(8) :: fob_min
    character(14) :: TCPLOT = "VARIABLES ="
    real(8) :: media_it_fob = 0.d0
    real(8) :: variancia_it_fob = 0.d0
    real(8) :: desvio_padrao_it_fob = 0.d0
    real(8) :: termo_quadrado= 0.d0
    
    !---Vari�veis usadas na execu��o do HSA
    type(HSA_par):: param
    type(HSA_sai),allocatable:: sHSA(:)
    integer:: ex   !Contador de execu��es do HSA
    integer:: NEA  !N�mero de execu��es do HSA
    real(8):: ffop !toler�ncia para definir �timo pr�tico: (fob-fot)/fot < ffop
    character(len=500):: tipo_problema !defini��o do tipo de problema a ser executado
    real(8),allocatable:: Rot(:), Rop(:)    !Confiabilidade aparente para o �timo e �timo pr�tico, resp.
    integer,allocatable:: Eot(:), Eop(:)    !N�m. de avalia��es da fun��o objetivo relativo a Rot e Rop, resp.
    
    !---Vari�veis auxiliares
    integer:: ii, iii, mm, kk, jj !contadores
    integer:: ncla
    character(len=500):: nomearquivo, texto
    
    !---Vari�veis definindo as unidades de entrada e sa�da de dados    
    integer,parameter:: uni_dad = 501, uni_print=502, uni_parametros = 503, uni_sai = 504, uni_Rot = 505, uni_Rop = 506, uni_result = 507, uni_estatistica = 510 
    
    
    
    !---------------------------------------------
    !defini��o das namelist para entrada dedados
    !---------------------------------------------
    namelist /parametros_HSA/ param
    namelist /execHSA/ tipo_problema, NEA, ffop
    ! -----------------------------------------------------------------------------------------------------------------------------
  
    !---------------------------------------------
    !In�cio das opera��es
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
    !abrindo arquivo e lendo os par�metros do HSA
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
    !abrindo arquivos de sa�da de informa��es e escrevendo par�metros do HSA usados na otimiza��o
    open( unit=uni_print , file=trim(nomearquivo)//"_saida_geral.txt" , action='WRITE',iostat=iii)
    open( unit=uni_parametros , file=trim(nomearquivo)//"_parHSA.txt" , action='WRITE',iostat=iii)
    param%uni_imp_par = uni_parametros
    write(uni_print,nml=parametros_HSA)
    open( unit=uni_sai , file=trim(nomearquivo)//".sHS" , action='WRITE',iostat=iii)
    open(unit=uni_Rot,file=trim(nomearquivo)//"_R.plt", action='WRITE',iostat=iii)
    open(unit=uni_Rop,file=trim(nomearquivo)//"_Rop.plt", action='WRITE',iostat=iii)
    !-----------------
    
    !-----------------
    !Alocando as vari�veis
    allocate( sHSA(NEA), Rot(NEA), Eot(NEA), Rop(NEA), Eop(NEA))
    !-----------------
    
    
    
    
    !-----------------
    !Executando as otimiza��es
    
    selectcase( trim(tipo_problema) )
    case('trelica1')    
        
        !-------------DADOS DE ENTRADA - TRELI�A TRAPEZOIDAL---------------!
    !Abertura de arquivo e leitura de dados---------------------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S1)
    close (uni_dad_vento)
    
    ! Chamada para c�lculo do fator topogr�fico S1 --------------------------------------------------------------------------------
    call Fator_Topografico_S1 (caso, theta, z, d, S1)
      
    !Leitura dos dados do problema referentes ao fator S2 -------------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Fator_Top_S2)
    close (uni_dad_vento)
    
    ! Chamada para c�lculo do fator topogr�fico S2 --------------------------------------------------------------------------------
    call Transf_dados_S2 (tab_S2)
    call Fator_Topografico_S2(cat, cls, z, tab_S2, S2)
    
    !Leitura dos dados do problema referentes a press�o din�mica ------------------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = pressao_dinamica)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes ao coeficiente de press�o -----------------------------------------------------------
    open (unit=uni_dad_vento , file="Vento_Entrada.txt" , action='read')
        read(uni_dad_vento, nml = Cofc_pressao)
    close (uni_dad_vento)
    
    !Leitura dos dados do problema referentes aos dados geom�tricos de entrada -----------------------------------------------------------
    open (unit = uni_dad_vento, file = 'Vento_Entrada.txt', action='read')        
        read (uni_dad_vento, nml = Dados_geometricos)
    close (unit = uni_dad_vento)
    
    !Leitura dos dados do problema referentes aos coeficientes de pondera��o e ao caso de carga de vento (0 ou 1) -------------------------
    open (unit = uni_dad_vento, file = 'Vento_Entrada.txt', action='read')        
        read (uni_dad_vento, nml = Coeficientes_Ponderacao)
    close (unit = uni_dad_vento)
    
    !Armazenamento dos dados do problema referentes aos v�rios perfis de ter�as  -----------------------------------------------------------
    open (unit = uni_dad_tercas, file = 'dados_tercas.txt', action='read')        
        call lista_tercas(uni_dad_tercas, terca)
    close (unit = uni_dad_tercas)
    
    !Armazenamento das propriedades geom�tricas referentes aos v�rios perfis 2L da AISC  -----------------------------------------------------------    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do

    select case(verif)
        case('HSA')

        do ex = 1 , NEA

            call iniciar_HSA_sai(param,sHSA(ex))
            call otimiza_HSA(funcao_objetivo,param,sHSA(ex))
            !call impr_ocl1(param%nVD,sHSA(ex)%ot%vd,texto)
            write(*,'(i3.3," fob=",ES14.6,2x,A)')ex,sHSA(ex)%ot%fob,trim(texto(param%nVD*2:))
            write(uni_sai,'(A,i2.2," fob= ",ES14.6,2x,A)') "ocl1 -> ex=",ex,sHSA(ex)%ot%fob,trim(texto)
            open( unit=uni_result , file="resultados_NEA_"//trim(str(ex))//".txt" , action='write',iostat=iii)
                write(uni_result, '(A, 7A7, 2A12)') TCPLOT, 'it', 'Vd(1)', 'Vd(2)', 'Vd(3)', 'Vd(4)', 'fob', 'pp', 'AF_ot','fob_ot'
                do i = 1, sHSA(ex)%IT
                    j = sHSA(ex)%OT_HAF(i)
                    write(uni_result, '(I20 , I6 , I8 , I6, I6, 2G15.6, I6, G15.6)') i, sHSA(ex)%HAF(j)%Vd(1:4), sHSA(ex)%HAF(j)%fob, peso_total(j), j, sHSA(ex)%HAF(j)%fob
                end do
                media_it_fob = media_it_fob + sHSA(ex)%IT_OT
                termo_quadrado = termo_quadrado + sHSA(ex)%IT_OT**2
                if(ex == NEA) then
                    media_it_fob = media_it_fob/NEA
                    variancia_it_fob = 1.d0/(NEA-1.d0)*(termo_quadrado - NEA*(media_it_fob)**2)
                    desvio_padrao_it_fob = sqrt(variancia_it_fob)
                    open(unit = uni_estatistica, file="dados_estatisticos.txt", action='write', iostat=iii)
                        write(uni_estatistica, '(A, A8, A12, A15)') TCPLOT, 'm�dia', 'vari�ncia','desvio-padr�o'
                        write(uni_estatistica, '(G20.6,G16.7,G16.7 )') media_it_fob, variancia_it_fob, desvio_padrao_it_fob
                    close(unit=uni_estatistica)
                    media_it_fob = 0.d0
                    termo_quadrado = 0.d0
                    desvio_padrao_it_fob = 0.d0
                end if
          end do

	  case('minimo_abs')

            allocate(Vd(4))
            write(uni_result, '(A, 7A7, A12)') TCPLOT, 'Vd(1)', 'Vd(2)', 'Vd(3)', 'Vd(4)', 'fob', 'pp', 'fob_min'    
            nd = 4
            nc = 0
            do ii = 1, 61
                do jj = 1, 61
                    do kk = 1, 61
                        do mm = 1, 61
                            Vd(1) = ii
                            Vd(2) = jj
                            Vd(3) = kk
                            Vd(4)=  mm
                            call funcao_objetivo(nd, nc, Vd, Vc, fob)
                                
                                
                            if (Vd(4) == 1 .AND. Vd(3) == 1 .AND. Vd(2) == 1 .AND. Vd(1) == 1) then
                                fob_min = fob
                            else if (fob < fob_min .AND. fob == peso_total_abs) then
                                fob_min = fob
                                end if
                                    
                            write(uni_result, '(I20, 3I6, 3G15.6)') Vd(1), Vd(2), Vd(3), Vd(4), fob, peso_total_abs, fob_min
                        end do
                    end do
                end do
            end do
            deallocate(Vd)
            
       end select         
            close(uni_result)
            num_AF = 0.d0
            
            !if(ex == 1) then
            !    inclinacao_diagonais = "/\"
            !else if(ex == 2) then
            !    inclinacao_diagonais = "\/"
            !end if
            
 !          
        

                    
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